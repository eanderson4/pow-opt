#include "isjn.h"

rgrid *  isjn::solveModel( isolve * is){

  rgrid * rg = new rgrid();
  /*
  time_t tstart;
  tstart = clock();
    
  int Nl = getGrid()->numBranches();
  gridcalc gc(getGrid());
  
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  int n=0;
  double large=0;

  //  cerr<<"n\tyn\tzn\tan\tbn\tun\tvn\n";
  cout<<"Check: "<<sum(_check)<<" / "<<getGrid()->numBranches()<<endl;

  if (cplex.solve()){
    bool systemfail=true;
    while(systemfail){
      n++; if (n>500) throw itlimit;  
      cout<<" ----- Iteration - "<<n<<" ----- "<<endl;

      IloNumArray fsolve(getEnv(),Nl);
      cplex.getValues(fsolve,getF());
      vec f0=gc.convert(fsolve);
      vec z0=gc.risk(f0,_var0,getL(),getP(),getPc());      

      IloNumArray gsolve(getEnv(),Nl);
      cplex.getValues(gsolve,getG());
      vec g0=gc.convert(gsolve);
      
      cout<<"Base system"<<endl;
      systemfail = postCC(f0,z0,&cplex,n);
      
      if(n==2) large=5000;
      for(int i=0;i<Nl;i++){
	if(_check(i) && i<(30 + large)){
	  cout<<"Contingency: "<<i<<endl;
	  vec fn = getN1(i,f0,g0);
	  vec zn=gc.risk(fn,_var.row(i).t(),getL(),getP(),getPc());
	  if(!fn.is_finite()) fn.print("fn: ");
	  cout<<"Post eval"<<endl;
	  systemfail += postN1(i,fn,g0,zn,&cplex,n);
	}
	else cout<<"dont check Contingency "<<i<<endl;
      }

      
      if(!systemfail && n>2) break;
      if(!cplex.solve()) break;

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
    
    cout<<"Iterations: "<<n<<endl;
    cout<<"Cuts: "<<getTotalCuts()<<endl;
    cout<<"Time: "<<total<<endl;
        
  }

  else{
    cerr<<"Not solved"<<endl;
    cerr<<cplex.getStatus()<<endl;
  }
  cplex.end();
  */  
  return rg;
}


void isjn::setup(){
  cout<<"Setup N-1 joint chance constraint"<<endl;
  stringstream ss;
  /*
  IloEnv env = getEnv();

  int Nl = getGrid()->numBranches();
  gridcalc gc(getGrid());
  _addCut = mat(Nl,Nl,fill::zeros);
  _check = vec(Nl,fill::ones);

  _C = gc.getC();
  _Cg = gc.getCm();
  cout<<"Get Branch Sensitivities"<<endl;
  _Hb=_Hw*trans(_C);
  cout<<"Build L"<<endl;
  _L=gc.getL(_Hw);


  cout<<"Here"<<endl;;
  mat Sig = getSig();
  _var0 = Sig.diag();
  _var=mat(Nl,Nl);
  
  for(int n=0;n<Nl;n++){  //Line i outage
    for(int j=0;j<Nl;j++){
      if(_Hb(n,n)<=1-.000001 || _Hb(n,n)>=1+.0000001){
	_var(n,j)=Sig(j,j) + 2*_L(j,n)*Sig(n,j) + pow(_L(j,n),2)*Sig(n,n);
      }
      else{
	_var(n,j)=0;
	//	_var(n,j)=Sig(j,j) + 2*_Hb(j,n)*Sig(n,j) + pow(_Hb(j,n),2)*Sig(n,n);
	_check(n)=0;
      }
    }
  }
  cout<<sum(_check)<<" / "<<getGrid()->numBranches()<<endl;
  cout<<"there"<<endl;

  _in = mat(Nl,Nl,fill::zeros);

  _riskConstraint = IloRangeArray(env,Nl,0,_epsN);
  for(int n=0;n<Nl;n++){
    _z.push_back(IloNumVarArray(env,1,0,IloInfinity));
    _yplus.push_back(IloNumVarArray(env,1,0,IloInfinity));
    _yup.push_back(IloRangeArray(env, 1,0,IloInfinity));
    _ydown.push_back(IloRangeArray(env,1,0,IloInfinity));

  }
  getModel()->add(_riskConstraint);  
  cout<<"DONT BUILDING"<<endl;
  */
}

    /*
        _riskConstraint[n].setExpr( IloSum(_z[n]) );
        ss.str("");
        ss<<"rc"<<n<<"[0,"<<getEps()<<"]";
        _riskConstraint[n].setName( ss.str().c_str() );
    for(int e=0;e<Nl;e++){
	ss.str("");
	ss<<"z"<<n<<","<<e<<"[0,"<<_epsN<<"]";
	_z[n][e].setName( ss.str().c_str() );
      }
	  
        getModel()->add(_z[n]);
        getModel()->add(_riskConstraint[n]);
	        getModel()->add(_fplus[n]);  /// save
		getModel()->add(_fup);  //save
		getModel()->add(_fdown); ///save
    */



bool isjn::postN1(int n, vec f,vec g, vec z, IloCplex * cplex, int iteration){
  /*  if(_check(n)!=1) return false;
  //define tolerance for line risk > 0
  stringstream ss;
  grid * gr = getGrid();
  double tol = pow(10,-5);
  ranvar rv;
  int Nl = getGrid()->numBranches();



  //  f.t().print("f: ");
  //  z.t().print("z: ");
  //Calculate system risk
  double r = sum(z);
  cout<<"Risk: "<<r<<endl;
  if(r<=_epsN+tol) return false;  // SYSTEM SUCCEEDED
  else{
    cout<<"CUTTING ----------"<<endl;
    for(int i=0; i<Nl; i++){
      if (z(i)>tol){
	if(sum(_in.row(n))==0){
	  _z[n] = IloNumVarArray(getEnv(),Nl,0,IloInfinity);	  
	  _yplus[n] = IloNumVarArray(getEnv(),Nl,0,IloInfinity);	  
	  _yup[n] = IloRangeArray(getEnv(),Nl,0,IloInfinity);	  
	  _ydown[n] = IloRangeArray(getEnv(),Nl,0,IloInfinity);	  
	}

	if(_Hb(n,n)<=1-.000001 || _Hb(n,n)>=1+.0000001){
	  _yup[n][i].setExpr( _yplus[n][i] - (getF()[i] + _L(i,n)*(getF()[n])) );
	  _ydown[n][i].setExpr( _yplus[n][i] + (getF()[i] + _L(i,n)*getF()[n]) );
	}


	double U=getGrid()->getBranch(i).getRateA();
	double Ueps = rv.ginv(_epsN,getL(),getP(),getPc());
	_yplus[n][i].setBounds(0, U*Ueps);
	ss.str("");
	ss<<"y"<<n<<","<<i<<"[0,"<<U*Ueps<<"]";
	_yplus[n][i].setName( ss.str().c_str() );

	getModel()->add(_yplus[n][i]);
	getModel()->add(_yup[n][i]);
	getModel()->add(_ydown[n][i]);
	

	//Add cuts for each line with positive risk
	double y_i = abs(f(i));
	//	cout<<"z_"<<i<<": "<<z(i)<<", y_"<<i<<": "<<y_i<<endl;
	double dz=rv.deriveMu(getL(),getP(),getPc(),abs(f(i))/U,sqrt(getSig()(i,i))/U);
	//	cout<<"z_"<<i<<" >= "<<dz/U<<" y_"<<i<<" + "<<(z(i)-dz*y_i/U)<<endl;
	IloRange cut(getEnv(),-IloInfinity,0);
	cut.setExpr( dz*(_yplus[n][i] - y_i)/U + z(i) - _z[n][i]);
	//	cout<<cut<<endl;
	getModel()->add(cut);
	_addCut(n,i)=_addCut(n,i)+1;

	if(_in(n,i)==0){
	  _riskConstraint[n].setExpr( _riskConstraint[n].getExpr() + _z[n][i] );
	  _in(n,i)=1;
	}

	if(i==28)
	  {
	    double a=z(i)-dz*y_i/U;
	    double b=dz/U;
	    double u=1/U;
	    double v=b;
	    //	  cerr<<iteration<<"\t"<<y_i/U<<"\t"<<z(i)<<"\t"<<z(i)-dz*y_i/U<<"\t"<<dz/U<<"\t"<<u<<"\t"<<v<<endl;
	  }

	if(z(i)>.5) {
	  f.t().print("f: ");
	  z.t().print("z: ");	  
	  //	  throw;
	}

      }
    }
    return true;
  }     
  */
  return true;
}


vec isjn::getN1(int n, vec y0, vec g){
  grid * gr = getGrid();
  vec yn;
  /*
  if(_Hb(n,n)<=1-.000001 || _Hb(n,n)>=1+.0000001){
    //    cout<<"Line "<<n<<endl;
    yn = y0[n]*_L.col(n)+y0;
  }
  else{
    cout<<"Line: "<<n<<" --- HB = 1"<<endl;  //obvious island ?
    
    branch br = gr->getBranch(n);
    int from = gr->getBusNum(br.getFrom());
    int to = gr->getBusNum(br.getTo());
    int winner =-1;
    cout<<from<<" -> "<<to<<endl;
    
    if(sum(abs(_C.col(from))) == 1) {
      cout<<"winner winner: "<<from<<endl;
      winner=from;
    }
    else if(sum(abs(_C.col(to))) == 1) {
      cout<<"winner winner: "<<to<<endl;
      winner=to;
    }
    else {
      _C.col(from).t().print("from: ");
      _C.col(to).t().print("to: ");
      _L.col(n).t().print("L: ");
      throw nowin;
    }
    
    bus b = getGrid()->getBus(winner);
    double load = b.getP()+b.getGs();

    if(sum(_Cg.row(winner))>=1){
      cout<<"Have Generation"<<endl;
      vec gnode = _Cg.row(winner)*g;
      yn = _Hw.col(winner)*(load-gnode) + y0;
    }
    else{
      cout<<"Loadshed "<<load<<" at "<<winner<<endl;
      yn = _Hw.col(winner)*load+y0;
    }
    
  }
  */  
  return yn;

}
