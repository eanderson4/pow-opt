#include "test.h"

void test::run(){
  double etol=pow(10,-5); 
  bool fail=false;
  
  //Display Grid Info
  _gr->buildMap();
  _gr->printNums(cout);
  cout<<*_gr<<endl;
 
  try {
    //    PTDF(etol);
    LODF(etol);
  }  
  catch (IloException& e){
    cerr<<"Concert exception caught: "<<e<<endl; fail=true;
  }
  catch (exception& e){
    cerr<<"Exception: "<<e.what()<<endl; fail=true;
  }
  catch (...){
    cerr<<"Unknown exception caught"<<endl; fail=true;
  }
  
  cout<<"\n\n\nTest Over"<<endl;
  if(!fail) cout<<"No Problems"<<endl;


}


void test::LODF(double etol){ 
  grid * gr = _gr;
  int Nb=gr->numBuses();
  int Nl=gr->numBranches();

  //Solve Base System
  gridcalc gc(gr);
  igrid ig(gr);
  ig.allowLoadShed();
  ig.addCost();
  isolve is;
  is.setSolver(IloCplex::Dual,IloCplex::Dual);
  rgrid * rbase = ig.solveModel(&is);
  IloNumArray g_nom=rbase->getG();
  vec fbase= gc.convert(rbase->getF());
  
  //Set Slack
  IloNumArray slack(IloEnv(),Nb);
  for(int i=0;i<Nb;i++){
    slack[i]=0;
  }
  slack[1]=1;
  ig.addSlack(g_nom,slack);
  mat Hw = gc.getHw(gc.convert(slack));
  Hw.print("Hw: ");
  mat L = gc.getL(Hw);
  L.print("L: ");
  
  
  
  rgrid * rg;
  vector<rgrid *> _rg;
  del_g dg(gr);
  cout<<"N-1 contingencies [ ";
  for(int i=0;i<Nl;i++){
    dg.setStatus(i,false);
    cout<<dg<<" ";
    ig.modGrid( dg );
    rg = ig.solveModel();
    //    rg->displayOperatingPos(gr);
    rg->outputInfo(cout);
    _rg.push_back(rg);
    ig.unmodGrid( dg );
    dg.setStatus(i,true);
  }
  cout<<"]"<<endl;
  
  mat Hb=Hw*gc.getC().t();
  bool fail=false;
  for(int i=0;i<Nl;i++){
    //    cout<<Hb(i,i);
    if(Hb(i,i)<=1-.000001 || Hb(i,i)>=1+.0000001){
      cout<<"Line "<<i<<endl;
      vec fcalc = fbase[i]*L.col(i)+fbase;
      fcalc.print("fcalc: ");
      fail=false;
      for(int j=0;j<Nl;j++){
	double U = gr->getBranch(j).getRateA();
	if(fcalc[j]>=U)fail=true;  //need redispatch
	if(fcalc[j]<=-U)fail=true; //need redispatch
	}
      if(!fail){
	_rg[i]->outputInfo(cout);
	vec f = gc.convert(_rg[i]->getF());
	
	(f-fcalc).t().print("error: ");
	cout<<"sum: "<<sum(f-fcalc)<<endl;
	double error=sum(abs(f-fcalc));
	if(error>=etol) throw errtol;
      }
      else{
	//	  throw fse;
	cout<<"Need redispatch"<<endl;
      }
    }
  }
  
  //Build ilomodel out of armadillo matricies
  
  
  //Solve N-1 problem
  
}


