#include <stdlib.h>

#include "sqlinter.h"
#include "test.h"
#include "ijn1.h"

using namespace std;

int main(int argc, char* argv[]){

  if(argc<=2){
    cout<<"cmd: pow case/30.db <M>\n"
	<<"\trun main for case30"
	<<"\t<M> capacity multiplier"<<endl;
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
  

  double multi=atof(argv[2]);
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

  gridcalc gc(gr);  
  vec slack=gc.getSlack();

  mat Hw = gc.getHw(slack);

  mat SIGy(Nl,Nl,fill::zeros);
  SIGy = Hw*Cm*SIG*(Hw*Cm).t();

  //Set Probability info
  ranvar rv;
  double L=.9;
  double p=.01;
  double pc=.55;


  igrid ig(gr);
  ig.addCost();
  isolve is;
  is.setSolver(IloCplex::Dual,IloCplex::Dual);
  rgrid * rbase = ig.solveModel(&is);
  //  rbase->displayOperatingPos(gr);

  try{
    double eps=.05;
    double epsN=.08;
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
    cout<<"gen0: "<<rbase->getG()<<endl;
    cout<<"gen1: "<<rn1_1->getG()<<endl;
    //    gc.getL(Hw).col(35).print("L: ");
    //    cout<<f0(35)<<endl;


    running_stat<double> stats_r0;
    running_stat<double> stats_r1;

    vec check = n1.getCheck();
    for(int i=0;i<Nl;i++){
      if(check(i)==1){
	vec f0n = n1.getN1(i,f0,g0);
	vec z0n=gc.risk(f0n,SIGy.diag(),L,p,pc);
	double r0n = sum(z0n);
	stats_r0(r0n);
	vec f1n = n1.getN1(i,f1,g1);
	vec z1n=gc.risk(f1n,SIGy.diag(),L,p,pc);
	double r1n = sum(z1n);
	stats_r1(r1n);
      }
    }
    cout.precision(5);
    cout<<fixed<<endl;
    cout<<"C0: "<<o0<<endl;
    cout<<"r0 - "<<r0<<endl;
    cout << "count = " << stats_r0.count() << endl;
    cout << "mean = " << stats_r0.mean() << endl;
    cout << "stdv  = " << stats_r0.stddev()  << endl;
    cout << "min  = " << stats_r0.min()  << endl;
    cout << "max  = " << stats_r0.max()  << endl;
    cout<<endl;
    cout<<"C1: "<<o1<<endl;
    cout<<"r1 - "<<r1<<endl;
    cout << "count = " << stats_r1.count() << endl;
    cout << "mean = " << stats_r1.mean() << endl;
    cout << "stdv  = " << stats_r1.stddev()  << endl;
    cout << "min  = " << stats_r1.min()  << endl;
    cout << "max  = " << stats_r1.max()  << endl;
    

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

   
