#include <stdlib.h>

#include "sqlinter.h"
#include "igrid.h"
#include "grid.h"
#include "gridcalc.h"
#include "rv.h"

using namespace std;

int main(int argc, char* argv[]){

  if(argc<=1){
    cout<<"cmd: pow case/30.db\n"
	<<"\trun main for case30"<<endl;
    return 1;
  }

  //Load Grid
  sqlInter db;
  grid * gr = new grid;
  string db_name;
  db_name = argv[1];
  db.openDb(db_name);
  db.load(*gr);
 
  //Display Grid Info
  gr->buildMap();
  gr->printNums(cout);
  cout<<*gr<<endl;

  //Declare Parameters
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
    ranvar r(time(NULL)*i);
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
  gridcalc gc(gr);
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
  
  //Run against samples and collect data
  double samplemean=0;
  double samplestdv=0;
  double rsim=0;
  double ranal=0;
  double rerror;
  mat f(Nbr,samples);
  mat we(Nbr,samples);
  vec fbase = gc.convert(rbase->getF());
  vec fmean(Nbr,fill::zeros);
  vec wemean(Nbr,fill::zeros);
  vec weanal(Nbr,fill::zeros);
  vec weerror(Nbr,fill::zeros);
  vec fstdv(Nbr,fill::zeros);
  vec fanalstdv(Nbr,fill::zeros);
  vec fmeanerror;
  vec fstdverror(Nbr,fill::zeros);
  rgrid * rg;
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
      double U=.35*gr->getBranch(j).getRateA();
      we(j,n)=rv.g(abs(f(j,n))/U,L,p,pc);
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

  //Calculate analytic probabilities
  for(int j=0;j<Nbr;j++){
    double U=.35*gr->getBranch(j).getRateA();
    weanal(j)=rv.anaProb(L,p,pc,abs(fmean(j))/U,fstdv(j)/pow(U,1));
  }
  ranal=sum(weanal);
  rsim=sum(wemean);
  weerror=wemean-weanal;
  rerror=rsim-ranal;


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
  cout<<"we (simulation) : sum()="<<rsim<<endl;
  wemean.t().print();
  cout<<"we (analytic) : sum()="<<ranal<<endl;
  weanal.t().print();
  cout<<"we error(sim - ana) : sum()="<<rerror<<endl;
  weerror.t().print();

  return 0;
}
