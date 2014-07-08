#include <stdlib.h>

#include "sqlinter.h"
#include "test.h"
#include "ijn1.h"

using namespace std;


int main(int argc, char* argv[]){

  if(argc<=1){
    cout<<"cmd: pow case/30.db\n"
	<<"\trun main for case30"<<endl;
    return 1;
  }

  sqlInter db;
  grid * gr = new grid;
  string db_name;
  
  db_name = argv[1];
  db.openDb(db_name);
  db.load(*gr);
  gr->buildMap();
  gr->printNums(cout);

  cout<<*gr<<endl;
  
  grid gr2 = *gr;

    test t(&gr2);
    t.run();
    cout<<"\n\n\n"<<endl;

    //    return 0;

  double multi=.8;
  int Nl = gr->numBranches();
  for(int i=0;i<Nl;i++){
    branch bi = gr->getBranch(i);
    double U = bi.getRateA();
    gr->setCapacity(i,multi*U);
  }

  int num=5;
  int index[5] = { 1, 3, 6, 20, 28 };
  double stdv[5] = {3,.5,5,3.75,.3};

  int Nb = gr->numBuses();
  mat Cm(Nb,num,fill::zeros);
  vec mu(num,fill::zeros);
  mat SIG(num,num,fill::zeros);
  for(int i=0;i<num;i++){
    Cm(index[i],i)=1;
    SIG(i,i)=pow(stdv[i],2);
  }
  SIG(0,2)=-14;SIG(2,0)=-14;
  SIG(0,3)=-9;SIG(3,0)=-9;
  SIG(2,3)=15;SIG(3,2)=15;
  
  vec slack(Nb,fill::zeros);
  slack[1]=1;
  gridcalc gc(gr);
  mat Hw = gc.getHw(slack);

  mat SIGy(Nl,Nl,fill::zeros);
  SIGy = Hw*Cm*SIG*(Hw*Cm).t();

  //Set Probability info
  ranvar rv;
  double L=.9;
  double p=.05;
  double pc=.5;


  igrid ig(gr);
  ig.addCost();
  isolve is;
  is.setSolver(IloCplex::Dual,IloCplex::Dual);
  rgrid * rbase = ig.solveModel(&is);
  rbase->displayOperatingPos(gr);

  try{
    double eps=.035;
    double epsN=.175;
    ijn1 n1(gr, SIGy,Hw,L,p,pc,eps,epsN);
    //ijcc n1(gr,SIGy,L,p,pc,eps);
    n1.addCost();
    rgrid * rn1_1 = n1.solveModel(&is);

    vec f0=gc.convert(rbase->getF());
    vec f1=gc.convert(rn1_1->getF());
    vec g0=gc.convert(rbase->getG());
    vec g1=gc.convert(rn1_1->getG());
    f0.t().print("f0: ");
    f1.t().print("f1: ");
    vec z0=gc.risk(f0,SIGy.diag(),L,p,pc);
    double r0 = sum(z0);
    vec z1=gc.risk(f1,SIGy.diag(),L,p,pc);
    double r1 = sum(z1);
    z0.t().print("z0: ");
    z1.t().print("z1: ");
    cout.precision(4);
    cout<<fixed<<r0<<" - "<<r1<<endl;
    double o0=rbase->getObjective();
    double o1=rn1_1->getObjective();
    cout<<o0<<" - "<<o1<<endl;
    //    gc.getL(Hw).col(35).print("L: ");
    //    cout<<f0(35)<<endl;
  }
  catch(IloException& e){
    cerr<<"Concert exception: "<<e<<endl; 
  }
  catch(exception& e){
    cerr<<"Exception: "<<e.what()<<endl; 
  }
  catch(...){
    cerr<<"Unknown Error "<<endl;
  }

  return 0;


}

   
