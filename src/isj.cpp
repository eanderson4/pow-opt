#include "isj.h"

rgrid *  isj::solveModel( isolve * is){

  rgrid * rg = new rgrid();
  
  time_t tstart;
  tstart = clock();
    
  int Nl = getGrid()->numBranches();
  int Ng = getGrid()->numGens();
  int Nm = _Cm.n_cols;
  
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  int n=0;  

  if (cplex.solve()){
    bool systemfail=true;
    while(systemfail){
      n++; if (n>100) throw itlimit;
      cout<<"Iteration: "<<n<<endl;
      cout<<"Checking Risk Constraint"<<endl;
      IloNumArray xsolve(getEnv(),Ng);
      IloNumArray betasolve(getEnv(),Ng);
      IloNumArray ysolve(getEnv(),Nl);
      IloNumArray sdsolve(getEnv(),Nl);

      cplex.getValues(xsolve,getG());
      cplex.getValues(ysolve,getF());
      cplex.getValues(betasolve,_beta);
      cplex.getValues(sdsolve,_sd);
      double genCost = getIcost().getCost(getGrid(),xsolve);
      vec x=_gc->convert(xsolve);      
      vec y=_gc->convert(ysolve);      
      vec beta=_gc->convert(betasolve);      
      vec sd=_gc->convert(sdsolve);      
      vec ones(Nm,1,fill::ones);

      cout<<"Cost: "<<genCost<<endl;
      x.t().print("x: ");
      y.t().print("y: ");
      beta.t().print("beta: ");
      sd.t().print("sd: ");
      
      vec pi = getA()*getCg()*beta;
      //      pi.t().print("pi: ");
      mat term = getA()*(getCg()*beta*ones.t() - getCm());    
      mat SIGy = term*_SIG*term.t();
      //      SIGy.print("sigy: ");
      vec z=_gc->risk(y,SIGy.diag(),_L,_p,_pc);
      z.t().print("z: ");
      //      cout<<"Risk: "<<sum(z)<<endl;

      vec sdtest(Nl);
      for(int e=0;e<Nl;e++){
	sdtest(e) = sqrt( pi(e)*pi(e)*_sig_delta - 2 * pi(e)*_sig(e) + _sigger(e,e));
      }
      sdtest.print("sdtest: ");
      

      systemfail = postCC(y,z,SIGy.diag(),&cplex);

      if(systemfail) cplex.solve();

    }
    //Risk constraint satisfied, record solution
    float total= float(clock() - tstart) / CLOCKS_PER_SEC;  
    cout<<"\n - Solve info -"<<endl;
    cout<<"MODEL solved in "<<total<<endl;
    rg->getSolveInfo(&cplex,total);
    getBaseResults(&cplex, rg);
    if(isLoadShed()) getIshed().getLoadShed(&cplex, rg);
    
    cout<<"STATUS: "<<rg->getStatus()<<endl;
    cout<<"OBJECTIVE: "<<cplex.getObjValue()<<"\n"<<endl;
    double genCost = getIcost().getCost(getGrid(),rg->getG());
    rg->setGenCost(genCost);
    
    cout<<"\n\nRisk constraint satisfied\n\n"<<endl;
    cout<<"Iterations: "<<n<<endl;
    cout<<"Time: "<<total<<endl;
    _addCut.t().print("cuts: ");
  }

  else{
    cerr<<"Not solved"<<endl;
    cerr<<cplex.getStatus()<<endl;
  }
  cplex.end();
  
  return rg;

}      

bool isj::postCC(vec y, vec z,vec SIGy,IloCplex * cplex, int iteration){
  stringstream ss;
  //define tolerance for line risk > 0
  double tol = pow(10,-5);
  int Nl = getGrid()->numBranches();
  int Ng = getGrid()->numGens();
  double account=0;
  
  //  IloNumArray yplus(getEnv(),Nl);
  //  cplex->getValues(fp,getYplus());
  
  ranvar rv;
  double Ueps = rv.ginv(_eps,_L,_p,_pc);
  double r = sum(z);
  cout<<"Risk: "<<r<<endl;
  if(r<=_eps+tol) return false;
  else{
    cout<<"CUTTING ----------"<<endl;
    for(int i=0; i<Nl; i++){
      if (z(i)>0){
       	account += z(i);
	if(_addCut(i)==0){
	  _riskConstraint.setExpr( _riskConstraint.getExpr() + _z[i]);
	  _yup[i].setExpr( _yplus[i] - getF()[i] );
	  _ydown[i].setExpr( _yplus[i] + getF()[i] );
	  double U = getGrid()->getBranch(i).getRateA();
	  _yplus[i].setBounds(0,U*Ueps);	  
	  getModel()->add(_yup[i]);
	  getModel()->add(_ydown[i]);
	  ss.str("");
	  ss<<"yplus"<<i<<"[0,"<<U*Ueps<<"]";
	  _yplus[i].setName( ss.str().c_str() );
	  
	  _pibeta[i].setExpr(_pi[i]);
	  for(int j=0;j<Ng;j++){
	    _pibeta[i].setExpr( _pibeta[i].getExpr() - _A(i,j)*_beta[j]);
	  }
	  getModel()->add(_pibeta[i]);
	  ss.str("");
	  ss<<"pibeta"<<i<<"[0,inf]";
	  _pibeta[i].setName( ss.str().c_str() );

	  cout<<_sig(i)<<endl;
	  _sdfe[i].setExpr( -_pi[i]*_pi[i]*_sig_delta + 2*_pi[i]*_sig(i) - _sigger(i,i) + _sd[i]*_sd[i] );
	  cout<<"sdfe"<<i<<": "<<_sdfe[i]<<endl;
	  getModel()->add(_sdfe);
	  ss.str("");
	  ss<<"sdfe"<<i<<"[0,inf]";
	  _sdfe[i].setName( ss.str().c_str() );
	  
	}
	//Add cuts for each line with positive risk
	double y_i = abs(y(i));
	cout<<"z_"<<i<<": "<<z(i)<<", y_"<<i<<": "<<y_i<<endl;
	double U=getGrid()->getBranch(i).getRateA();
	double dmu=rv.deriveMu(_L,_p,_pc,y_i/U,sqrt(SIGy(i))/U);
	double dsigma=rv.deriveSigma(_L,_p,_pc,y_i/U,sqrt(SIGy(i))/U);
	cout<<"dmu: "<<dmu<<", dsigma: "<<dsigma<<endl;
	//	cout<<dz<<" "<<dz/U<<endl;
	//	cout<<"z_"<<i<<" >= "<<dz<<"(y_"<<i<<" - "<<y_i<<")/"<<U<<" + "<<z(i)<<endl;
	cout<<"z_"<<i<<" >= "<<dmu/U<<" y_"<<i<<" + "<<z(i)-dmu*y_i/U<<endl;

	IloRange cut(getEnv(),-IloInfinity,0);
	cut.setExpr( dmu/U*(_yplus[i] - y_i) + dsigma*(_sd[i] - sqrt(SIGy(i))/U) + z(i) - _z[i]);
	cout<<cut<<endl;
	getModel()->add(cut);
	_addCut(i)=_addCut(i)+1;
	
	//Output cut information
	/*	if(i==28){
	  double a=z(i)-dz*y_i/U;
	  double b=dz/U;
	  double u=1/U;
	  double v=b;
	  //	  cerr<<iteration<<"\t"<<y_i/U<<"\t"<<z(i)<<"\t"<<z(i)-dz*y_i/U<<"\t"<<dz/U<<"\t"<<u<<"\t"<<v<<endl;
	  }*/

	if(z(i)==1){
	  y.t().print("y: ");
	  throw;
	}
	
      }
    }
    cout<<"Accounted: "<<account<<endl;
    return true;
  }     

  return true;
}

void isj::lineLimitStatus(bool status){
  int Nl = getGrid()->numBranches();
  
  ranvar rv;
  double Ueps = rv.ginv(_eps,_L,_p,_pc);

  for(int i=0;i<Nl;i++){
    double U = getGrid()->getBranch(i).getRateA();
    if(!status){
      getF()[i].setBounds(-U*Ueps,U*Ueps);
    }
    else{
      double U = getGrid()->getBranch(i).getRateA();
      getF()[i].setBounds(-U,U);
    }
  }


}


void isj::setup(){
  cout<<"Setup slack joint chance constraint"<<endl;
  stringstream ss;

  ranvar rv;
  double Ueps = rv.ginv(_eps,_L,_p,_pc);

  IloEnv env = getEnv();
  int Nl = getGrid()->numBranches();
  int Ng = getGrid()->numGens();
  int Nm = _Cm.n_cols;
  _addCut = vec(Nl,fill::zeros);
    _A = getA();
  //Calculate branch variance terms --------------------

  cout<<"Build Static Variance Terms: "<<endl;

  mat Ag = getA()*getCg();
  mat Ak = getA()*getCm();
  vec ones2(Nm,1,fill::ones);
  _sig_delta = accu(_SIG);
  _sig = Ak*_SIG*ones2;
  _sigger = Ak*_SIG*Ak.t();
  

    
  //Build CPLEX model ----------------------------------
  cout<<"Build CPLEX model for base: "<<endl;

  _z = IloNumVarArray(env,Nl,0,_eps);
  _yplus = IloNumVarArray(env,Nl,0,IloInfinity);
  _sd = IloNumVarArray(env,Nl,0,IloInfinity);
  _beta = IloNumVarArray(env,Ng,0,1);
  _pi = IloNumVarArray(env,Nl,-IloInfinity,IloInfinity);

  _yup = IloRangeArray(env, Nl,0,IloInfinity);
  _ydown = IloRangeArray(env,Nl,0,IloInfinity);
  _ydown = IloRangeArray(env,Nl,0,IloInfinity);
  _pibeta = IloRangeArray(env,Nl,0,0);
  _sdfe = IloRangeArray(env,Nl,0,IloInfinity);

  _riskConstraint = IloRange(env,0,_eps,"riskconstraint");
  _riskConstraint.setExpr( IloSum(_z) );
  _betaSum = IloRange(env,1,1,"betasum");
  _betaSum.setExpr( IloSum(_beta) );

  for(int i=0;i<Nl;i++){
    double U = getGrid()->getBranch(i).getRateA();
    getF()[i].setBounds(-U*Ueps,U*Ueps);

    ss.str("");
    ss<<"pi"<<i<<"[-inf,inf]";
    _pi[i].setName( ss.str().c_str() );
    ss.str("");
    ss<<"sd"<<i<<"[0,inf]";
    _sd[i].setName( ss.str().c_str() );
    ss.str("");
    ss<<"z"<<i<<"[0,eps]";
    _z[i].setName( ss.str().c_str() );
  }
  for(int j=0;j<Ng;j++){
    ss.str("");
    ss<<"bj"<<j<<"[0,1]";
    _pi[j].setName( ss.str().c_str() );
  }
  
  addCost(_beta,_sig_delta);


  getModel()->add(_sd);

  getModel()->add(_z);
  getModel()->add(_beta);
  getModel()->add(_riskConstraint);
  getModel()->add(_betaSum);



}

  /*
  for(int i=0;i<Nl;i++){
    _yup[i].setExpr( _yplus[i] - getF()[i] );
    _ydown[i].setExpr( _yplus[i] + getF()[i] );
    double U = getGrid()->getBranch(i).getRateA();
    getF()[i].setBounds(-U*Ueps,U*Ueps);
    _yplus[i].setBounds(0,U*Ueps);
    ss.str("");
    ss<<"yplus"<<i<<"[0,"<<U*Ueps<<"]";
    _yplus[i].setName( ss.str().c_str() );
    ss.str("");
    ss<<"z"<<i<<"[0,"<<U*Ueps<<"]";
    _z[i].setName( ss.str().c_str() );
    ss.str("");
    ss<<"fup"<<i;
    _yup[i].setName( ss.str().c_str() );
    ss.str("");
    ss<<"fdown"<<i;
    _ydown[i].setName( ss.str().c_str() );
    }*/
  //  _riskConstraint.setExpr( IloSum(_z) );
  
  //  getModel()->add(_z);
  //  getModel()->add(_yplus);
  //  getModel()->add(_yup);
  //  getModel()->add(_ydown);
