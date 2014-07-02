#include "ijcc.h"

rgrid *  ijcc::solveModel( isolve * is){

  rgrid * rg = new rgrid();
  
  time_t tstart;
  tstart = clock();
    
  //relax branch flow constraints
  //substitute absolute (branch flow )
  //add line risk z

  //solve at cut until risk constraint satisfied
  int Nl = getGrid()->numBranches();
  double tol = pow(10,-6);
  
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  int n=0;

  if (cplex.solve()){
    float tot= float(clock() - tstart) / CLOCKS_PER_SEC;  
    cout<<"\n - Solve info -"<<endl;
    cout<<"MODEL solved in "<<tot<<endl;
    rg->getSolveInfo(&cplex,tot);
    getBaseResults(&cplex, rg);
    if(isLoadShed()) getIshed().getLoadShed(&cplex, rg);
    
    cout<<"STATUS: "<<rg->getStatus()<<endl;
    cout<<"OBJECTIVE: "<<cplex.getObjValue()<<"\n"<<endl;
    double genCost = getIcost().getCost(getGrid(),rg->getG());
    rg->setGenCost(genCost);
    
    bool risk=false;
    gridcalc gc(getGrid());
    ranvar rv;
    while(!risk){
      n++; if (n>100) break;
      cout<<"Checking Risk Constraint"<<endl;
      vec f=gc.convert(rg->getF());

      vec z=gc.risk(f,_SIGy.diag(),_L,_p,_pc);
      double r = sum(z);
      cout<<"Risk: "<<r<<endl;
      if(r<=_eps){
	float total= float(clock() - tstart) / CLOCKS_PER_SEC;  
	cout<<"\n\nRisk constraint satisfied\n\n"<<endl;
	cout<<"Iterations: "<<n<<endl;
	cout<<"Time: "<<total<<endl;
	f.t().print("f: ");
	z.t().print("z: ");

	risk=true;
      }
      else{
	cout<<"Add cuts"<<endl;
	vector<int> m;
	vector<double> y_i;
	vector<double> z_i;
	vector<double> dz_i;
	for(int i=0; i<Nl; i++){
	  if (z(i)>tol){
	    cout<<i<<": "<<z(i)<<endl;
	    m.push_back(i);
	    y_i.push_back(abs(f(i)));
	    z_i.push_back(z(i));
	    double U=getGrid()->getBranch(i).getRateA();
	    double dz=rv.deriveMu(_L,_p,_pc,abs(f(i))/U,sqrt(_SIGy(i,i))/U);
	    dz_i.push_back(dz);
	  }
	}
	for(int i=0;i<(int)m.size();i++){
	  cout<<"z_"<<m[i]<<": "<<dz_i[i]<<"(y_"<<m[i]<<" - "<<y_i[i]<<") + "<<z_i[i]<<endl;
	  IloRange cut(getEnv(),-IloInfinity,0);
	  cut.setExpr( dz_i[i]*(_fplus[m[i]] - y_i[i]) + z_i[i] - _z[m[i]]);
	  cout<<cut<<endl;
	  getModel()->add(cut);
	}
	//cuts added, resolve
	time_t tstart2;
	tstart2 = clock();
	cplex.solve();
	float tot2= float(clock() - tstart2) / CLOCKS_PER_SEC;  
	cout<<"\n - Solve info -"<<endl;
	cout<<"MODEL solved in "<<tot2<<endl;
	rg->getSolveInfo(&cplex,tot2);
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
	//	cout<<fout<<endl;
	//	cout<<zout<<endl;

      }
    }
      
  }
  else{
    cerr<<"Not solved"<<endl;
    cerr<<cplex.getStatus()<<endl;
  }
  cplex.end();
  
  return rg;
}

void ijcc::lineLimitStatus(bool status){
  int Nl = getGrid()->numBranches();
  
  for(int i=0;i<Nl;i++){
    if(!status)    getF()[i].setBounds(-IloInfinity,IloInfinity);
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
