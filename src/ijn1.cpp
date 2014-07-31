#include "ijn1.h"

rgrid *  ijn1::solveModel( isolve * is){

  rgrid * rg = new rgrid();

  time_t tstart;
  tstart = clock();
    
  int Nl = getGrid()->numBranches();
  gridcalc gc(getGrid());
  
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  int n=0;

  cerr<<"n\tyn\tzn\tan\tbn\tun\tvn\n";

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

      for(int i=0;i<Nl;i++){
	cout<<"Contingency: "<<i<<endl;
	vec fn = getN1(i,f0,g0);
	vec zn=gc.risk(fn,_var.row(i).t(),getL(),getP(),getPc());
	cout<<"Post eval"<<endl;
	systemfail += postN1(i,fn,g0,zn,&cplex,n);
      }
      
      if(!systemfail) break;
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
  
  return rg;
}


void ijn1::setup(){
  cout<<"Setup N-1 joint chance constraint"<<endl;
  stringstream ss;

  IloEnv env = getEnv();

  int Nl = getGrid()->numBranches();
  gridcalc gc(getGrid());
  ranvar rv;
  _addCut = mat(Nl,Nl,fill::zeros);

  _C = gc.getC();
  _Cg = gc.getCm();
  _L=gc.getL(_Hw);
  _Hb=_Hw*_C.t();
  cout<<"Here";
  mat Sig = getSig();
  _var0 = Sig.diag();
  _var=mat(Nl,Nl);
  
  for(int n=0;n<Nl;n++){  //Line i outage
    for(int j=0;j<Nl;j++){
      if(_Hb(n,n)<=1-.000001 || _Hb(n,n)>=1+.0000001){
	_var(n,j)=Sig(j,j) + 2*_L(j,n)*Sig(n,j) + pow(_L(j,n),2)*Sig(n,n);
      }
      else{
	_var(n,j)=Sig(j,j) + 2*_Hb(j,n)*Sig(n,j) + pow(_Hb(j,n),2)*Sig(n,n);
      }
    }
  }
  cout<<"there";

  _riskConstraint = IloRangeArray(env,Nl,0,_epsN);
  for(int n=0;n<Nl;n++){
    _z.push_back(IloNumVarArray(env,Nl,0,IloInfinity));
    _yplus.push_back(IloNumVarArray(env,Nl,0,IloInfinity));
    _yup.push_back(IloRangeArray(env, Nl,0,IloInfinity));
    _ydown.push_back(IloRangeArray(env,Nl,0,IloInfinity));
   
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
	//        getModel()->add(_fplus[n]);  /// save
	//	getModel()->add(_fup);  //save
	//	getModel()->add(_fdown); ///save

  }
}


bool ijn1::postN1(int n, vec f,vec g, vec z, IloCplex * cplex, int iteration){
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
	if(_Hb(n,n)<=1-.000001 || _Hb(n,n)>=1+.0000001){
	  _yup[n][i].setExpr( _yplus[n][i] - (getF()[i] + _L(i,n)*(getF()[n])) );
	  _ydown[n][i].setExpr( _yplus[n][i] + (getF()[i] + _L(i,n)*getF()[n]) );
	}
	else{
	  if(n==12){
	     f.t().print("f: ");
	     z.t().print("z: ");
	  }
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
	  else	    throw nowin;
	  bus b = getGrid()->getBus(winner);
	  double load = b.getP()+b.getGs();
	  _yup[n][i].setExpr( _yplus[n][i] - (getF()[i] + _Hw(i,winner)*load) );
	  _ydown[n][i].setExpr( _yplus[n][i] + (getF()[i] + _Hw(i,winner)*load) );

	  if(sum(_Cg.row(winner))>=1){
	    for(int j=0; j<(int)_Cg.row(winner).n_elem;j++){
	      _yup[n][i].setExpr( _yup[n][i].getExpr() + _Cg(winner,j)*getG()[j]);
	      _ydown[n][i].setExpr( _yup[n][i].getExpr() - _Cg(winner,j)*getG()[j]);
	    }
	  }
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
	cout<<"z_"<<i<<": "<<z(i)<<", y_"<<i<<": "<<y_i<<endl;
	double dz=rv.deriveMu(getL(),getP(),getPc(),abs(f(i))/U,sqrt(getSig()(i,i))/U);
	cout<<"z_"<<i<<" >= "<<dz/U<<" y_"<<i<<" + "<<(z(i)-dz*y_i/U)<<endl;
	IloRange cut(getEnv(),-IloInfinity,0);
	cut.setExpr( dz*(_yplus[n][i] - y_i)/U + z(i) - _z[n][i]);
	cout<<cut<<endl;
	getModel()->add(cut);
	_addCut(n,i)=_addCut(n,i)+1;

	if(i==28)
{	  double a=z(i)-dz*y_i/U;
	  double b=dz/U;
	  double u=1/U;
	  double v=b;
	  cerr<<iteration<<"\t"<<y_i/U<<"\t"<<z(i)<<"\t"<<z(i)-dz*y_i/U<<"\t"<<dz/U<<"\t"<<u<<"\t"<<v<<endl;
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
}


vec ijn1::getN1(int n, vec y0, vec g){
  grid * gr = getGrid();
  vec yn;

  if(_Hb(n,n)<=1-.000001 || _Hb(n,n)>=1+.0000001){
    cout<<"Line "<<n<<endl;
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
  
  return yn;

}


double ijn1::getTotalCuts(){
  double total=0;
  int Nl = getGrid()->numBranches();
  total += getNumBaseCuts();
  for(int i=0;i<Nl;i++){
    total += getNumCuts(i);
  }
  return total;
}
