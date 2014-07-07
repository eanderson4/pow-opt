#include "ijcc.h"

rgrid *  ijcc::solveModel( isolve * is){

  rgrid * rg = new rgrid();
  
  time_t tstart;
  tstart = clock();
    
  int Nl = getGrid()->numBranches();
  
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  gridcalc gc(getGrid());
  int n=0;
    

  if (cplex.solve()){
    bool systemfail=true;
    while(systemfail){
      n++; if (n>100) throw itlimit;
      cout<<"Checking Risk Constraint"<<endl;
      IloNumArray fsolve(getEnv(),Nl);
      cplex.getValues(fsolve,_fplus);
      vec f=gc.convert(fsolve);      
      vec z=gc.risk(f,_SIGy.diag(),_L,_p,_pc);

      systemfail = postCC(f,z,&cplex);
      cplex.solve();
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
    cplex.getValues(fout,_fplus);
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

bool ijcc::postCC(vec f, vec z,IloCplex * cplex){
  //define tolerance for line risk > 0
  double tol = pow(10,-6);
  int Nl = getGrid()->numBranches();
  
  ranvar rv;
  double r = sum(z);
  cout<<"Risk: "<<r<<endl;
  if(r<=_eps) return false;
  else{
    cout<<"Add cuts"<<endl;
    for(int i=0; i<Nl; i++){
      if (z(i)>tol){
	//Add cuts for each line with positive risk
	cout<<i<<": "<<z(i)<<endl;
	double y_i = abs(f(i));
	double U=getGrid()->getBranch(i).getRateA();
	double dz=rv.deriveMu(_L,_p,_pc,abs(f(i))/U,sqrt(_SIGy(i,i))/U);
	cout<<"z_"<<i<<" >= "<<dz<<"(y_"<<i<<" - "<<y_i<<") + "<<z(i)<<endl;
	IloRange cut(getEnv(),-IloInfinity,0);
	cut.setExpr( dz*(_fplus[i] - y_i) + z(i) - _z[i]);
	getModel()->add(cut);
      }
    }
    return true;
  }     
}

void ijcc::lineLimitStatus(bool status){
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


void ijcc::setup(){
  cout<<"Setup joint chance constraint"<<endl;
  stringstream ss;

  IloEnv env = getEnv();
  int Nl = getGrid()->numBranches();
  _z = IloNumVarArray(env,Nl,0,IloInfinity);
  _fplus = IloNumVarArray(env,Nl,0,IloInfinity);

  _fup = IloRangeArray(env, Nl,0,IloInfinity);
  _fdown = IloRangeArray(env,Nl,0,IloInfinity);
  
  for(int i=0;i<Nl;i++){
    _fup[i].setExpr( _fplus[i] - getF()[i] );
    _fdown[i].setExpr( _fplus[i] + getF()[i] );

      ss.str("");
      ss<<"fplus"<<i<<"[0,inf]";
      _fplus[i].setName( ss.str().c_str() );
      ss.str("");
      ss<<"z"<<i<<"[0,inf]";
      _z[i].setName( ss.str().c_str() );
      ss.str("");
      ss<<"fup"<<i;
      _fup[i].setName( ss.str().c_str() );
      ss.str("");
      ss<<"fdown"<<i;
      _fdown[i].setName( ss.str().c_str() );
  }

  _riskConstraint = IloRange(env,0,_eps,"riskconstraint");
  _riskConstraint.setExpr( IloSum(_z) );
  
  getModel()->add(_z);
  getModel()->add(_fplus);
  getModel()->add(_fup);
  getModel()->add(_fdown);
  getModel()->add(_riskConstraint);
}
