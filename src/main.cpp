#include <stdlib.h>

#include "sqlinter.h"
#include "test.h"

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

  //    test t(gr);
  //    t.run();
  //    cout<<"\n\n\n"<<endl;

  int Nb = gr->numBuses();
  int Nl = gr->numBranches();
  int samples=175000;
  int num=5;
  int index[5] = { 1, 3, 6, 20, 28 };
  double stdv[5] = {3,.5,5,3.75,.3};
  
  arma_rng::set_seed_random();

  vec mu(num,fill::zeros);
  mat SIG(num,num,fill::zeros);
  SIG(0,0)=pow(3,2);
  SIG(1,1)=pow(.5,2);
  SIG(2,2)=pow(5,2);
  SIG(3,3)=pow(3.75,2);
  SIG(4,4)=pow(.3,2);

  mat Cm(Nb,num,fill::zeros);
  for(int i=0;i<num;i++){
    Cm(index[i],i)=1;
  }
  
       SIG(0,2)=-14;SIG(2,0)=-14;
          SIG(0,3)=-9;SIG(3,0)=-9;
          SIG(2,3)=15;SIG(3,2)=15;
  
      mat R=chol(SIG);

  //Set Probability info
  ranvar rv;
  double L=.9;
  double p=.15;
  double pc=.5;

  //Solve Base System
  gridcalc gc(gr);
  igrid ig(gr);
  ig.addCost();
  isolve is;
  is.setSolver(IloCplex::Dual,IloCplex::Dual);
  rgrid * rbase = ig.solveModel(&is);
  rbase->displayBranchPos(gr);
  vec fbase = gc.convert(rbase->getF());
  
  //Set Slack
  vec slack(Nb,fill::zeros);
  slack[1]=1;
  mat Hw = gc.getHw(slack);

  //Analytic Branch Flow Statistics
  mat SIGy(Nl,Nl,fill::zeros);
  SIGy = Hw*Cm*SIG*(Hw*Cm).t();
    vec fstdv(Nl,fill::zeros);
    for(int i=0;i<num;i++){
      fstdv=fstdv + pow(stdv[i]*Hw.col(index[i]),2);
    }
    //    fstdv=sqrt(fstdv);
    fstdv.print("var: ");
    SIGy.diag().print("diag s: ");
    (fstdv-SIGy.diag()).print("error: ");

  vec zprime(Nl,fill::zeros);
  vec zcheck(Nl,fill::zeros);
  //Calculate analytic probabilities
  for(int j=0;j<Nl;j++){
    double U=.75*gr->getBranch(j).getRateA();
    //    cout<<sqrt(SIGy(j,j))/U<<endl;
    zprime(j)=rv.anaProb(L,p,pc,abs(fbase(j))/U,sqrt(SIGy(j,j))/U);
    zcheck(j)=rv.anaProb(L,p,pc,abs(fbase(j))/U,sqrt(fstdv(j))/U);
    cout<<zprime(j)<<" - "<<zcheck(j)<<endl;
  }
  //  zprime.print("zprime: ");
  double rprime = sum(zprime);
  double rcheck = sum(zcheck);
  cout<<rprime<<" - "<<rcheck<<endl;
  


  vec r(samples,fill::zeros);

  
  running_stat_vec<vec> rsv_in(true);
  running_stat_vec<vec> rsv_delx(true);
  running_stat_vec<vec> rsv_y;
  cout<<"Run samples"<<endl;
  for(int n=0;n<samples;n++){
    vec randin(num,fill::randn);
    rsv_in(randin);
    vec x = mu + R.t()*randin;
    rsv_delx(x);
    vec delp(Nb,fill::zeros);
    delp = Cm*x;
    vec delf = Hw*delp;
    vec f = fbase + delf;
    rsv_y(f);

    vec z(Nl,fill::zeros);
    for(int j=0;j<Nl;j++){
      double U=.75*gr->getBranch(j).getRateA();
      z(j)=rv.g(abs(f(j))/U,L,p,pc);
    }
    r(n)=sum(z);
  }
  
  double ravg=sum(r)/samples;
  double rstdv=stddev(r);

  double error = ravg - rprime;

  std::cout << std::fixed;
  std::cout << std::setprecision(5);

  cout<<"RESULTS"<<endl;
  cout<<"sim - analytic"<<endl;
  cout<<ravg<<" - "<<rprime<<endl;
  cout<<error<<endl;
  cout<<"\n";
  cout<<"mean\t+- stdv\t\t(stderr)"<<endl;
  cout<<ravg<<"\t+-"<<rstdv<<"\t( "<<rstdv/sqrt(samples)<<" )"<<endl;

  return 0;
}




  /*
  cout<<"\n\nCheck stats"<<endl;
  
    rsv_delx.mean().print("delx mean: ");
  rsv_delx.cov().print("delx cov: ");
  rsv_y.var().print("y var: ");


    
  double delta=.01;
  double Up=10000;
  double a=.15;

  double mean=.9;
  double sigma=.3;

  double total=0;
  double expect=0;

  for(int i=1;i*delta<Up;i++){
    double z=i*delta;

    double alpha_plus = (sqrt(z/a) - mean)/sigma;
    double dalpha = 1/2/sigma/sqrt(a*z);

    double expr=rv.phi(alpha_plus)*dalpha;
    
    total = total + expr*delta;
    expect=expect+ z*expr*delta;
    
    }

  double check = rv.PHI(-mean/sigma);
  
  cout<<total<<" - "<<check<<endl;
  cout<<"mu: "<<expect<<endl;

  //N-1 calculations
  int s=100000;
  vec rn(s,fill::randn);
  for(int i=0;i<s;i++){
    rn(i)=rn(i)*rn(i);
  }
  double test = sum(rn)/s;
  cout<<test<<endl;
  */

   
