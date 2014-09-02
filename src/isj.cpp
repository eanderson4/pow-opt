#include "isj.h"

rgrid *  isj::solveModel( isolve * is){

  rgrid * rg = new rgrid();
  
  time_t tstart;
  tstart = clock();
    
  int Nl = getGrid()->numBranches();
  int Ng = getGrid()->numGens();
  //  int Nm = _Cm.n_cols;
  
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  int n=0;  

  if (cplex.solve()){
    bool systemfail=true;
    while(systemfail){
      n++; if (n>100) throw itlimit;
      cout<<"Checking Risk Constraint"<<endl;
      IloNumArray xsolve(getEnv(),Ng);
      IloNumArray betasolve(getEnv(),Ng);
      IloNumArray ysolve(getEnv(),Nl);

      cplex.getValues(xsolve,getG());
      cplex.getValues(ysolve,getF());
      cplex.getValues(betasolve,_beta);
      vec x=_gc->convert(xsolve);      
      vec y=_gc->convert(ysolve);      
      vec beta=_gc->convert(betasolve);      
      vec ones(Ng,1,fill::ones);

      
      vec pi = getA()*getCg()*beta*ones.t();
      mat SIGy = getA()*(getCg()*beta*ones.t() - getCm())*_SIG*(getCg()*beta*ones.t() - getCm())*getA().t();

      

      vec z=_gc->risk(y,SIGy.diag(),_L,_p,_pc);

      systemfail = postCC(y,z,&cplex);
      z.t().print("z: ");
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
    
    IloNumArray zout(getEnv(),Nl);
    IloNumArray fout(getEnv(),Nl);
    cplex.getValues(zout,_z);
    cplex.getValues(fout,_yplus);
    cout<<"\n\nRisk constraint satisfied\n\n"<<endl;
    cout<<"Iterations: "<<n<<endl;
    cout<<"Time: "<<total<<endl;
    cout<<"f: "<<fout<<endl;
    cout<<"z: "<<zout<<endl;    

  }

  else{
    cerr<<"Not solved"<<endl;
    cerr<<cplex.getStatus()<<endl;
  }
  cplex.end();
  
  return rg;

}      

bool isj::postCC(vec f, vec z,IloCplex * cplex, int iteration){
  /*  //define tolerance for line risk > 0
  double tol = pow(10,-5);
  int Nl = getGrid()->numBranches();
  double account=0;
  
  IloNumArray yplus(getEnv(),Nl);
  cplex->getValues(fp,getYplus());
  
  ranvar rv;
  double r = sum(z);
  cout<<"Risk: "<<r<<endl;
  if(r<=_eps+tol) return false;
  else{
    cout<<"CUTTING ----------"<<endl;
    for(int i=0; i<Nl; i++){
      if (z(i)>0){
	account += z(i);
	//Add cuts for each line with positive risk
	double y_i = abs(y(i));
	cout<<"z_"<<i<<": "<<z(i)<<", y_"<<i<<": "<<y_i<<endl;
	double U=getGrid()->getBranch(i).getRateA();
//-SIGm	double dz=rv.deriveMu(_L,_p,_pc,abs(f(i))/U,sqrt(_SIGy(i,i))/U);
	//	cout<<dz<<" "<<dz/U<<endl;
	//	cout<<"z_"<<i<<" >= "<<dz<<"(y_"<<i<<" - "<<y_i<<")/"<<U<<" + "<<z(i)<<endl;
	cout<<"z_"<<i<<" >= "<<dz/U<<" y_"<<i<<" + "<<z(i)-dz*y_i/U<<endl;
	IloRange cut(getEnv(),-IloInfinity,0);
	cut.setExpr( dz/U*(_fplus[i] - y_i) + z(i) - _z[i]);
	cout<<cut<<endl;
	getModel()->add(cut);
	_addCut(i)=_addCut(i)+1;
	
	if(i==28){
	  double a=z(i)-dz*y_i/U;
	  double b=dz/U;
	  double u=1/U;
	  double v=b;
	  //	  cerr<<iteration<<"\t"<<y_i/U<<"\t"<<z(i)<<"\t"<<z(i)-dz*y_i/U<<"\t"<<dz/U<<"\t"<<u<<"\t"<<v<<endl;
	}

	if(z(i)==1){
	  f.t().print("f: ");
	  cout<<"fplus: "<<fp<<endl;
	  cout<<"f: "<<getF()<<endl;
	  throw;
	}
	
      }
    }
    cout<<"Accounted: "<<account<<endl;
    return true;
  }     
  */
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
  //Calculate branch variance terms --------------------

  _sig_delta = accu(_SIG);
  _sig = vec(Nl);
  _sigger = mat(Nl,Nl);

  time_t tstart;
  tstart = clock();

  mat Ag(Nl,Ng,fill::zeros);
  for(int e=0;e<Nl;e++){
    for(int j=0;j<Ng;j++){
      Ag(e,j) = dot(getA().row(e),getCg().col(j)); 
    }
  }
  //  Ag.print("Ag: ");

  mat Ak(Nl,Nm,fill::zeros);
  for(int e=0;e<Nl;e++){
    for(int k=0;k<Nm;k++){
      Ak(e,k) = dot(getA().row(e),getCm().col(k)); 
    }
  }
  //  Ak.print("Ak: ");
  float total= float(clock() - tstart) ;  

  time_t tstart2;
  tstart2 = clock();


  mat Ag2 = getA()*getCg();
  (Ag - Ag2).print("error: ");

  mat Ak2 = getA()*getCm();
  (Ak - Ak2).print("error: ");

    float total2= float(clock() - tstart2);  

    cout<<"2 is better if >0: "<<total-total2<<endl;

  for(int e=0;e<Nl;e++){
    double total=0;
    for(int k=0;k<Nm;k++){
      for(int k2=0;k2<Nm;k2++){
	total=total+Ak(e,k)*_SIG(k,k2);
      }
    }
    _sig(e)=total;
  }
  for(int e1=0; e1<Nl;e1++){
    for(int e2=0;e2<Nl;e2++){
      double t2=0;
      for(int k1=0;k1<Nm;k1++){
	for(int k2=0;k2<Nm;k2++){
	  t2=t2+Ak(e1,k1)*Ak(e2,k2)*_SIG(k1,k2);
	}
      }
      _sigger(e1,e2)=t2;
    }
  }
  

  //Build CPLEX model ----------------------------------

  _z = IloNumVarArray(env,Nl,0,_eps);
  _yplus = IloNumVarArray(env,Nl,0,IloInfinity);
  _sd = IloNumVarArray(env,Nl,0,IloInfinity);
  _beta = IloNumVarArray(env,Nl,0,IloInfinity);
  _pi = IloNumVarArray(env,Nl,0,IloInfinity);

  _yup = IloRangeArray(env, Nl,0,IloInfinity);
  _ydown = IloRangeArray(env,Nl,0,IloInfinity);
  _ydown = IloRangeArray(env,Nl,0,IloInfinity);
  _pibeta = IloRangeArray(env,Nl,0,0);
  _sdfe = IloRangeArray(env,Nl,0,IloInfinity);

  _riskConstraint = IloRange(env,0,_eps,"riskconstraint");
  getModel()->add(_riskConstraint);

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
