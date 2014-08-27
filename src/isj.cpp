#include "isj.h"

rgrid *  isj::solveModel( isolve * is){
  rgrid * bad;
  throw nogo;
  return bad;
}      

bool isj::postCC(vec f, vec z,IloCplex * cplex, int iteration){
  //define tolerance for line risk > 0
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

  _z = IloNumVarArray(env,Nl,0,IloInfinity);
  _yplus = IloNumVarArray(env,Nl,0,IloInfinity);

  _yup = IloRangeArray(env, Nl,0,IloInfinity);
  _ydown = IloRangeArray(env,Nl,0,IloInfinity);
  
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
  }

  _riskConstraint = IloRange(env,0,_eps,"riskconstraint");
  _riskConstraint.setExpr( IloSum(_z) );
  
  getModel()->add(_z);
  getModel()->add(_yplus);
  getModel()->add(_yup);
  getModel()->add(_ydown);
  getModel()->add(_riskConstraint);

}
