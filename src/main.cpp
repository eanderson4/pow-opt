#include <stdlib.h>

#include "sqlinter.h"
#include "igrid.h"
#include "grid.h"
#include "gridcalc.h"
#include "rv.h"

using namespace std;

void test(grid * gr){
  gridcalc gc(gr);

  del_g dg(gr);
  dg.addDemand(4,5);
  dg.addDemand(6,7);
  
  cout<<"dg: "<<gc.getDelG(dg).t()<<endl;


  int Nb=gr->numBuses();
  vec delg(Nb,fill::zeros);
  vec slackdist(Nb,fill::zeros);
  delg(4)=5; delg(6)=7;
  slackdist(1)=1;

  vec del_f = gc.getDelF(delg,slackdist);


  igrid ig(gr);
  ig.addCost();
  rgrid * rg;
  rgrid * rg2;

  rg = ig.solveModel();
  IloNumArray g_nom=rg->getG();
  cout<<"Before"<<endl;
  cout<<"F: "<<rg->getF()<<endl;
  cout<<"G: "<<g_nom<<endl;
  

  ig.modGrid(dg);
    
  IloNumArray slack(IloEnv(),Nb);
  for(int i=0;i<Nb;i++) slack[i]=0;
  slack[1]=1;
  
  ig.addSlack(g_nom,slack);

  rg2 = ig.solveModel();
  cout<<"After"<<endl;
  cout<<"F: "<<rg2->getF()<<endl;
  cout<<"G: "<<rg2->getG()<<endl;

  cout<<"Total Demand: "<<gr->getTotalDemand()<<endl;
  cout<<"Total Gen: "<<IloSum(rg2->getG())<<endl;
  
  cout<<"delF (sim): "<<endl;
  for(int i=0;i<gr->numBranches();i++)
    cout<<" "<<(rg->getF()[i] - rg2->getF()[i])<<"   ";

  cout<<"\ndelF (sim) - delF (shift factor)"<<endl;
  double totalerror=0;
  for(int i=0;i<gr->numBranches();i++){
    cout<<" "<<(rg->getF()[i] - rg2->getF()[i] - del_f(i))<<"   ";
    totalerror=totalerror+(rg->getF()[i] - rg2->getF()[i] - del_f(i));
  }
  cout<<"\n";
  cout<<"Total Error: "<<totalerror<<endl;

  ranvar rv;
  rv.testRV();  

}

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
  gridcalc gc(gr);


  int Nb = gr->numBuses();
  int Nbr = gr->numBranches();

  int samples=750;
  vector<ranvar> rv_bus;
  int num=5;
  int index[5] = { 1, 3, 6, 20, 28 };
  double mean[5] = {0,0,0,0,0}; //p_0 = {21.7,7.6,22.8,17.5,2.4}
  double stdv[5] = {3,.5,5,3.75,.3};
  double totalmean=0;
  double totalstdv=0;

  //Sample for random variables
  cout<<"\nCreate samples for random del_g"<<endl;
  for(int i=0;i<num;i++){
    ranvar r(1052*i);
    r.createRV(samples,mean[i],stdv[i]);
    cout<<i<<":\t"<<mean[i]<<", "<<stdv[i]<<endl;
    cout<<"\t"<<r.sampleMean()<<", "<<r.sampleStDv()<<endl;
    rv_bus.push_back(r);
    totalmean=totalmean+mean[i];
    totalstdv=totalstdv+pow(stdv[i],2);
  }
  totalmean=totalmean/num;
  totalstdv=sqrt(totalstdv);
  

  //Solve Base System
  igrid ig(gr);
  ig.addCost();
  isolve is;
  is.setSolver(IloCplex::Dual,IloCplex::Dual);
  rgrid * rbase = ig.solveModel(&is);
  IloNumArray g_nom=rbase->getG();

  //Set Slack
  IloNumArray slack(IloEnv(),Nb);
  slack[1]=1;
  ig.addSlack(g_nom,slack);


  //Set Probability info
  ranvar rv;
  double L=.9;
  double p=.15;
  double pc=.5;
  
  rgrid * rg;
  //Run against samples and collect data
  double samplemean=0;
  double samplestdv=0;
  mat f(Nbr,samples);
  mat we(Nbr,samples);
  vec fbase = gc.convert(rbase->getF());
  vec fmean(Nbr,fill::zeros);
  vec wemean(Nbr,fill::zeros);
  vec fstdv(Nbr,fill::zeros);
  vec fanalstdv(Nbr,fill::zeros);
  vec fmeanerror;
  vec fstdverror(Nbr,fill::zeros);
  for(int n=0;n<samples;n++){
    del_g dg(gr);
    double total=0;
    for(int i=0;i<num;i++){
      dg.addDemand( index[i], rv_bus[i].getValue(n) );      
      total=total+rv_bus[i].getValue(n);
    }    
    samplemean=samplemean+total;
    ig.modGrid(dg);
    rg=ig.solveModel(&is);
    cout<<rg->getG()<<endl;
    ig.removeDemand(dg);
    f.col(n) = gc.convert(rg->getF());
    for(int j=0;j<Nbr;j++){
      we(j,n)=rv.g(f(j,n),L,p,pc);
    }
  }
  samplemean=samplemean/samples;
  //Calculate StDv
  for(int n=0; n<samples;n++){
    double total=0;
    for(int i=0;i<num;i++){
      total=total+rv_bus[i].getValue(n);
    }    
    samplestdv=samplestdv+pow(total-samplemean,2);
  }
  samplestdv=sqrt(samplestdv/(samples-1));

  //Calculate Branch Flow Statistics
  for(int i=0;i<Nbr;i++){
    fmean(i)=sum(f.row(i))/samples;
    wemean(i)=sum(we.row(i))/samples;
  }
  for(int i=0;i<Nbr;i++){
    fstdv(i)=sum(pow(f.row(i) - fmean(i),2));
    fstdv(i)=sqrt(fstdv(i)/(samples-1));
  }

  //Analytic Branch Flow Statistics
  fmeanerror= fmean - gc.convert(rbase->getF());
  mat Hw = gc.getHw(gc.convert(slack));
  
  for(int i=0;i<num;i++){
    fanalstdv=fanalstdv + pow(stdv[i]*Hw.col(index[i]),2);
  }
  fanalstdv=sqrt(fanalstdv);
  fanalstdv.print("sigma: ");

  fstdverror=fstdv-fanalstdv;


  cout<<"\n\nReport\n"<<endl;
  
  cout<<"del (total - sample) = error"<<endl;
  cout<<"mean: "<<totalmean<<" - "<<samplemean<<" = "<<totalmean-samplemean<<endl;
  cout<<"stdv: "<<totalstdv<<" - "<<samplestdv<<" = "<<totalstdv-samplestdv<<endl;

  cout<<"\nBranch Flow Statistics"<<endl;
  cout<<"Simulated"<<endl;
  fmean.t().print("Fmean: ");
  fstdv.t().print("Fstdv: ");

  cout<<"Analytic"<<endl;
  cout<<"Fmean: "<<rbase->getF()<<endl;
  cout<<"Fstdv: "<<fanalstdv.t()<<endl;

  cout<<"Error (Simulated - Analytic)"<<endl;
  cout<<"f mean error: sum()="<<sum(fmeanerror)<<endl;
  cout<<fmeanerror.t()<<endl;
  cout<<"f stdv error: sum()="<<sum(fstdverror)<<endl;
  cout<<fstdverror.t()<<endl;

  cout<<"Risk Analysis"<<endl;
  cout<<"we (simulation)"<<endl;
  wemean.print();
  cout<<"we (analytic)"<<endl;
  

  return 0;
}

//Analytic stdv not work
/*
  for(int i=0;i<Nbr;i++){
    for(int k=0;k<num;k++){
      fanalstdv(i)=pow(gc.getHw(gc.convert(slack))(i,index[k]),2)*pow(stdv[k],2);
      fanalstdv(i)=sqrt(fanalstdv(i));
    }
  }
*/ 
    //    cout<<"del_f sim:\n"<<f.col(n) - fbase<<endl;
    //    cout<<"del_f ana:\n"<<gc.getDelF(gc.getDelG(dg),gc.convert(slack))<<endl;
    //    cout<<"del_f error (sim - ana)"<<f.col(n) - fbase + gc.getDelF(gc.getDelG(dg),gc.convert(slack))<<endl;
    //    cout<<"n: "<<n<<"\n"<<dg<<endl;    
