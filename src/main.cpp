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
  int samples=75;
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
  for(int i=0;i<Nb;i++){
    slack[i]=0;
  }
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
  vec fanalmean(Nbr,fill::zeros);
  vec fmeanerror;
  vec fstdverror(Nbr,fill::zeros);
  vec fnormsim;
  vec fnormanal;
  vec fnormstdvsim;
  vec fnormstdvanal;
  rgrid * rg;
  mat Hw = gc.getHw(gc.convert(slack));
  for(int n=0;n<samples;n++){
    del_g dg(gr);
    double total=0;
    //Adjust Demand
    for(int i=0;i<num;i++){
      dg.addDemand( index[i], rv_bus[i].getValue(n) );      
      total=total+rv_bus[i].getValue(n);
    }    
    samplemean=samplemean+total;
    //Modify Implementation and Solve
    ig.modGrid(dg);
    rg=ig.solveModel(&is);
    cout<<rg->getG()<<endl;
    //Undo Modifications and store solve information
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
  fanalmean=gc.convert(rbase->getF());
  
  for(int i=0;i<num;i++){
    fanalstdv=fanalstdv + pow(stdv[i]*Hw.col(index[i]),2);
  }
  fanalstdv=sqrt(fanalstdv);
  fanalstdv.print("sigma: ");


  //Calculate analytic probabilities
  for(int j=0;j<Nbr;j++){
    double U=.35*gr->getBranch(j).getRateA();
    weanal(j)=rv.anaProb(L,p,pc,abs(fanalmean(j))/U,fanalstdv(j)/pow(U,1));
  }
  ranal=sum(weanal);
  rsim=sum(wemean);
  weerror=wemean-weanal;
  rerror=rsim-ranal;

  fnormsim=fmean;
  fnormanal=fanalmean;
  fnormstdvsim=fstdv;
  fnormstdvanal=fanalstdv;
  for(int i=0;i<Nbr;i++){
    double U=.35*gr->getBranch(i).getRateA();
    fnormsim(i)=fnormsim(i)/U;
    fnormanal(i)=fnormanal(i)/U;
    fnormstdvsim(i)=fnormstdvsim(i)/U;
    fnormstdvanal(i)=fnormstdvanal(i)/U;
  }
  fmeanerror=fnormsim - fnormanal;
  fstdverror=fnormstdvsim - fnormstdvanal;

  //Report Results
  cout<<"\n\nReport\n"<<endl;
  
  cout<<"del (total - sample) = error"<<endl;
  cout<<"mean: "<<totalmean<<" - "<<samplemean<<" = "<<totalmean-samplemean<<endl;
  cout<<"stdv: "<<totalstdv<<" - "<<samplestdv<<" = "<<totalstdv-samplestdv<<endl;

  cout<<"\nBranch Flow Statistics"<<endl;
  cout<<"Simulated"<<endl;
  fnormsim.t().print("Fmean: ");
  fnormstdvsim.t().print("Fstdv: ");

  cout<<"Analytic"<<endl;
  cout<<"Fmean: "<<fnormanal.t()<<endl;
  cout<<"Fstdv: "<<fnormstdvanal.t()<<endl;

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

  cout<<"\n\n\n"
      <<"Hessian Calculations for w"<<endl;

  p=.01;
  double eps=.0001;
  double epssig=.00001;
  double mustar=atof(argv[2]);
  double sigstar=atof(argv[3]);
  double w=rv.anaProb(L,p,pc,mustar,sigstar);
  double wmuup=rv.anaProb(L,p,pc,mustar+eps,sigstar);
  double wmudown=rv.anaProb(L,p,pc,mustar-eps,sigstar);
  double wsigmaup=rv.anaProb(L,p,pc,mustar,sigstar+epssig);
  double wsigmadown=rv.anaProb(L,p,pc,mustar,sigstar-epssig);
  double w1=rv.anaProb(L,p,pc,mustar+eps,sigstar+epssig);
  double w2=rv.anaProb(L,p,pc,mustar+eps,sigstar-epssig);
  double w3=rv.anaProb(L,p,pc,mustar-eps,sigstar+epssig);
  double w4=rv.anaProb(L,p,pc,mustar-eps,sigstar-epssig);
  double mufinitederiv=(wmuup-wmudown)/2/eps;
  double sigmafinitederiv=(wsigmaup-wsigmadown)/2/(epssig);
  double muanalderiv=rv.deriveMu(L,p,pc,mustar,sigstar);
  double sigmaanalderiv=rv.deriveSigma(L,p,pc,mustar,sigstar);
  double mufinitederiv2=(wmuup-2*w +wmudown)/eps/eps;
  double muanalderiv2=rv.d2Mu(L,p,pc,mustar,sigstar);
  double sigmafinitederiv2=(wsigmaup-2*w +wsigmadown)/epssig/epssig;
  double sigmaanalderiv2=rv.d2Sigma(L,p,pc,mustar,sigstar);
  double analderivesigmamu=rv.dSigmaMu(L,p,pc,mustar,sigstar);
  double finitederivesigmamu=(w1-w2-w3+w4)/4/eps/epssig;
  double d2mu=muanalderiv2;
  double d2sig=sigmaanalderiv2;
  double dsigmu=analderivesigmamu;
  double determinant=d2mu*d2sig-dsigmu*dsigmu;

  cout<<"\nCheck Functions"<<endl;
  cout<<"Mu: "<<mustar<<", Sigma: "<<sigstar<<endl;
  cout<<"\n partial w / partial mu"<<endl;
  cout<<"Finite Difference: "<<mufinitederiv<<endl;
  cout<<"Analytic: "<<muanalderiv<<endl;
  cout<<"\n partial w / partial sigma"<<endl;
  cout<<"Finite Difference: "<<sigmafinitederiv<<endl;
  cout<<"Analytic: "<<sigmaanalderiv<<endl;
  cout<<"\n partial2 w / partial mu2"<<endl;
  cout<<"Finite Difference: "<<mufinitederiv2<<endl;
  cout<<"Analytic: "<<muanalderiv2<<endl;
  cout<<"\n partial2 w / partial sigma2"<<endl;
  cout<<"Finite Difference: "<<sigmafinitederiv2<<endl;
  cout<<"Analytic: "<<sigmaanalderiv2<<endl;
  cout<<"\n partial2 w / partial sigma partial mu"<<endl;
  cout<<"Finite Difference: "<<finitederivesigmamu<<endl;
  cout<<"Analytic: "<<analderivesigmamu<<endl;

  cout<<"\n\n\tHessian"<<endl;
  cout<<setprecision(4)<<fixed;
  cout<<"\t"<<d2mu<<"\t"<<dsigmu<<endl;
  cout<<"\t"<<dsigmu<<"\t"<<d2sig<<endl;
  
  cout<<"\nDeterminant"<<endl;
  cout<<determinant<<endl;

  double a,b,c;
  a=1;
  b=-(d2mu+d2sig);
  c=d2mu*d2sig - dsigmu*dsigmu;
  double lambda1,lambda2;
  lambda1 = (-b + sqrt( b*b - 4*a*c ))/2/a;
  lambda2 = (-b - sqrt( b*b - 4*a*c ))/2/a;
  
  cout<<"Eigenvalues"<<endl;
  cout<<lambda1<<endl;
  cout<<lambda2<<endl;

  cout<<"\n";
  mat H(2,2);
  H(0,0)=d2mu;
  H(1,0)=dsigmu;
  H(0,1)=dsigmu;
  H(1,1)=d2sig;
  H.print("Hessian");

  vec eigval;
  mat eigvec;
  eig_sym(eigval,eigvec,H);
  eigval.print("eig val");
  eigvec.print("eig vec");
  
  return 0;
}
