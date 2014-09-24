#include "isj.h"

rgrid *  isj::solveModel( isolve * is){

  rgrid * rg = new rgrid();

  double tol = pow(10,-6);
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
      IloNumArray zsolve(getEnv(),Nl);

      cplex.getValues(xsolve,getG());
      cplex.getValues(ysolve,getF());
      cplex.getValues(betasolve,_beta);
      cplex.getValues(sdsolve,_sd);
      cplex.getValues(zsolve,_z);
      double genCost = getIcost().getCost(getGrid(),xsolve);
      vec x=_gc->convert(xsolve);      
      vec y=_gc->convert(ysolve);      
      vec beta=_gc->convert(betasolve);      
      vec sd=_gc->convert(sdsolve);      
      vec zsol=_gc->convert(zsolve);      
      vec ones(Nm,1,fill::ones);

      cout<<"Cost: "<<genCost<<endl;
      
      
      vec pi = getA()*getCg()*beta;
      mat term = getA()*(getCg()*beta*ones.t() - getCm());    
      mat SIGy = term*_SIG*term.t();
      //      SIGy.print("sigy: ");
      vec z=_gc->risk(y,SIGy.diag(),_L,_p,_pc);
      x.t().print("x: ");
      y.t().print("y: ");
      beta.t().print("beta: ");
      pi.t().print("pi: ");
      sd.t().print("sd: ");
      z.t().print("z: ");
      zsol.t().print("z sol: ");
      cout<<"r: "<<accu(z)<<endl;


      if( accu(z) < _eps+tol){
	systemfail=false;
      }
      else  systemfail = postCC(y,z,beta,SIGy.diag(),&cplex);

      if(systemfail) cplex.solve();
      _betaSolve = beta;
      _sdSolve = SIGy.diag();
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

bool isj::postCC(vec y, vec z, vec beta, vec SIGy,IloCplex * cplex, int iteration){
  stringstream ss;
  //define tolerance for line risk > 0
  double tol = pow(10,-6);
  int Nl = getGrid()->numBranches();

  double account=0;
  
  
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
	cout<<"Line "<<i<<endl;
	if(_addCut(i)==0){
	  cout<<"Initialize Cutting Variables for \n \tLine "<<i<<endl;
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
	  
	}
	//Add cuts for each line with positive risk
	double y_i = abs(y(i));
	//	cout<<"z_"<<i<<": "<<z(i)<<", y_"<<i<<": "<<y_i<<endl;
	double U=getGrid()->getBranch(i).getRateA();
	double dmu=rv.deriveMu(_L,_p,_pc,y_i/U,sqrt(SIGy(i))/U);
	double dsigma=rv.deriveSigma(_L,_p,_pc,y_i/U,sqrt(SIGy(i))/U);
	IloRange cut(getEnv(),-IloInfinity,0);
	cut.setExpr( dmu/U*(_yplus[i] - y_i) + dsigma/U*(_sd[i] - sqrt(SIGy(i))) + z(i) - _z[i]);
	cout<<cut<<endl;
	getModel()->add(cut);
	_addCut(i)=_addCut(i)+1;
	
	//Add cuts to describe standard deviation of branch flow
	IloRange cut_sd(getEnv(),-IloInfinity,0);

	int Ng = getGrid()->numGens();
	double pi_i = dot(getA().row(i),getCg()*beta);
	double sd_i =  sqrt( pi_i*pi_i*_sig_delta - 2 * pi_i*_sig(i) + _sigger(i,i));	
	cut_sd.setExpr( sd_i - _sd[i] );
	double term = (pi_i*_sig_delta - _sig(i))/sd_i;
	cout<<"\n";
	for( int j=0; j<Ng;j++){
	  //double pf_j = term;
	  double pf_j = _A(i,_indexG(j))*term;
	  cut_sd.setExpr( cut_sd.getExpr() + pf_j*(_beta[j]-beta(j)) );
	}
	cout<<cut_sd<<endl;
	getModel()->add(cut_sd);
	
	cout<<"\n\n";
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
  cout<<"--sj-- ===Setup=== slack joint chance constraint"<<endl;
  stringstream ss;

  ranvar rv;
  double Ueps = rv.ginv(_eps,_L,_p,_pc);
  cout<<"------: "<<Ueps<<endl;

  IloEnv env = getEnv();
  int Nl = getGrid()->numBranches();
  int Nb = getGrid()->numBuses();
  int Ng = getGrid()->numGens();
  int Nm = _indexM.n_elem;

  mat Cm(Nb,Nm,fill::zeros);

  for(int i=0;i<Nm;i++){
    Cm(_indexM(i),i)=1;
  }
  _Cm = Cm;

  _indexG = _gc->getIndexG();

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

  _yup = IloRangeArray(env, Nl,0,IloInfinity);
  _ydown = IloRangeArray(env,Nl,0,IloInfinity);
  _genup = IloRangeArray(env, Ng,0,IloInfinity);
  _gendown = IloRangeArray(env,Ng,0,IloInfinity);

  _riskConstraint = IloRange(env,0,_eps,"riskconstraint");
  //  _riskConstraint.setExpr( IloSum(_z) );
  _betaSum = IloRange(env,1,1,"betasum");
  _betaSum.setExpr( IloSum(_beta) );

  for(int i=0;i<Nl;i++){
    double U = getGrid()->getBranch(i).getRateA();
    getF()[i].setBounds(-U*Ueps,U*Ueps);

    ss.str("");
    ss<<"sd"<<i<<"[0,inf]";
    _sd[i].setName( ss.str().c_str() );
    ss.str("");
    ss<<"z"<<i<<"[0,eps]";
    _z[i].setName( ss.str().c_str() );
  }
  cout<<"GenInfo (epsG="<<_epsG<<")"<<endl;

  for(int j=0;j<Ng;j++){
    ss.str("");
    ss<<"beta"<<j<<"[0,1]";
    _beta[j].setName( ss.str().c_str() );
    double pmax = getGrid()->getGen(j).getPmax();
    double pmin = getGrid()->getGen(j).getPmin();
    double eta;
    if(_epsG == 1) eta=0;
    else eta = rv.PHIInverse(1-_epsG);
    //    rv.demoInverse();
    cout<<"gen "<<j<<": "<<pmax<<endl;
    _genup[j].setExpr(getG()[j] - _beta[j]*sqrt(_sig_delta)*eta );
    _genup[j].setBounds(0,pmax);
    //    _gendown[j].setExpr( getG()[j] - _beta[j]*sqrt(_sig_delta)*eta );
    //    _gendown[j].setBounds(pmin,IloInfinity);

    cout<<"eta*sqrt(TV): "<<eta*sqrt(_sig_delta)<<endl;
    cout<<_genup[j]<<endl;
    getModel()->add(_genup[j]);
  }
  
  addCost(_beta,_sig_delta);

  getModel()->add(_sd);

  getModel()->add(_z);
  getModel()->add(_beta);
  getModel()->add(_genup);
  //  getModel()->add(_gendown);
  getModel()->add(_riskConstraint);
  getModel()->add(_betaSum);


}
