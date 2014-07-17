#include "test.h"

void test::run(){
  double etol=pow(10,-5); 
  bool fail=false;
  
  //Display Grid Info
  _gr->buildMap();
  _gr->printNums(cout);
  cout<<*_gr<<endl;
 
  try {
        RV(etol);
	PTDF(etol);
	//        RSLACK(etol);
	//        LODF(etol);
	//   RANDOM( .1 );
	//        HESSIAN( etol*100 );
	//    COVAR( etol*100 );
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
  else cout<<"FAILED FAILED FAILED FAILED FAILED FAILED FAILED FAILED FAILED FAILED FAILED FAILED FAILED"<<endl;
}

void test::COVAR(double etol){
  grid * gr = _gr;
  int Nb = gr->numBuses();
  int Nl = gr->numBranches();
  int samples=77500;
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
  mat Ln = gc.getL(Hw);
  //Analytic Branch Flow Statistics
  mat SIGy(Nl,Nl,fill::zeros);
  SIGy = Hw*Cm*SIG*(Hw*Cm).t();
    vec fstdv(Nl,fill::zeros);
    for(int i=0;i<num;i++){
      fstdv=fstdv + pow(stdv[i]*Hw.col(index[i]),2);
    }
    //    fstdv=sqrt(fstdv);  //DOESNT INCLUDE COVARIANCE CALCULATION
    //    fstdv.print("var: ");
    //SIGy.diag().print("diag s: ");
    //(fstdv-SIGy.diag()).print("error: ");
    vec mean5 = fbase + Ln.col(5)*fbase(5);
    vec sig5 = SIGy.diag();
    for( int i=0;i<Nl;i++){
      sig5(i)=sig5(i) + 2*Ln(i,5)*SIGy(5,i) + pow(Ln(i,5),2)*SIGy(5,5);
    }
    
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
  running_stat_vec<vec> rsv_y5(true);
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
    vec f5 = f[5]*Ln.col(5)+f;
    rsv_y5(f5);

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

  rsv_delx.mean().t().print("dx mean: ");
  rsv_delx.cov().print("dx cov: ");
  SIG.print("SIG: ");
  cout<<"error: "<<sum(rsv_delx.cov().diag()-SIG.diag())<<endl;

  std::cout << std::fixed;
  std::cout << std::setprecision(5);

  cout<<"RESULTS"<<endl;
  cout<<"sim - analytic"<<endl;
  cout<<ravg<<" - "<<rprime<<endl;
  cout<<error<<endl;
  cout<<"\n";
  cout<<"mean\t+- stdv\t\t(stderr)"<<endl;
  cout<<ravg<<"\t+-"<<rstdv<<"\t( "<<rstdv/sqrt(samples)<<" )"<<endl;

  if(error>=etol) throw errtol;  


  vec errormean = rsv_y5.mean() - mean5;
  vec errorvar = rsv_y5.cov().diag() - sig5;
  rsv_y5.mean().t().print("y5 mean: ");
  mean5.t().print("y5 exp: ");
  cout<<"error (mean - expected): "<<sum(errormean)<<endl;
  errormean.t().print("mean error: ");
  rsv_y5.cov().diag().t().print("y5 cov: ");
  sig5.t().print("sig5: ");
  cout<<"error (var - expected): "<<sum(errorvar)<<endl;
  (rsv_y5.cov().diag() -sig5).t().print("error: ");
  
}


void test::HESSIAN(double etol){
  ranvar rv;
  cout<<"\n\n\n"
      <<"Hessian of w"<<endl;

  double L=.90;
  double p=.15;
  double pc=.4;

  double eps=.0001;
  double epssig=.000001;
  double mustar=.93;
  double sigstar=.01;
  double wmuup=rv.anaProb(L,p,pc,mustar+eps,sigstar);
  double wmu=rv.anaProb(L,p,pc,mustar,sigstar);
  double wmudown=rv.anaProb(L,p,pc,mustar-eps,sigstar);
  double wsigmaup=rv.anaProb(L,p,pc,mustar,sigstar+epssig);
  double wsigma=rv.anaProb(L,p,pc,mustar,sigstar);
  double wsigmadown=rv.anaProb(L,p,pc,mustar,sigstar-epssig);
  double wupup=rv.anaProb(L,p,pc,mustar+eps,sigstar+epssig);
  double wdownup=rv.anaProb(L,p,pc,mustar-eps,sigstar+epssig);
  double wupdown=rv.anaProb(L,p,pc,mustar+eps,sigstar-epssig);
  double wdowndown=rv.anaProb(L,p,pc,mustar-eps,sigstar-epssig);
  
  double mufinitederiv=(wmuup-wmudown)/2/eps;
  double sigmafinitederiv=(wsigmaup-wsigmadown)/2/(epssig);
  double mufinited2=(wmuup-2*wmu+wmudown)/eps/eps;
  double sigmafinited2=(wsigmaup-2*wsigma+wsigmadown)/epssig/epssig;
  double musigmad=(wupup-wupdown-wdownup+wdowndown)/4/eps/epssig;
  double sigmamuderiv=(wupup-wmuup-wsigmaup+2*wsigma-wmudown-wsigmadown+wdowndown)/2/eps/epssig;
  

	       
  double muanalderiv=rv.deriveMu(L,p,pc,mustar,sigstar);
  double sigmaanalderiv=rv.deriveSigma(L,p,pc,mustar,sigstar);
  double md2=rv.d2Mu(L,p,pc,mustar,sigstar);
  double sd2=rv.d2Sigma(L,p,pc,mustar,sigstar);
  double smd=rv.dSigmaMu(L,p,pc,mustar,sigstar);

  double emud=abs(mufinitederiv-muanalderiv);
  double esigmad=abs(sigmafinitederiv-sigmaanalderiv);
  double emd2=abs(mufinited2-md2);
  double esd2=abs(sigmafinited2-sd2);
  double esmd=abs(musigmad-smd);

  cout<<"\nCheck Functions"<<endl;
  cout<<"Mu: "<<mustar<<", Sigma: "<<sigstar<<endl;
  cout<<"\n partial w / partial mu"<<endl;
  cout<<"Finite Difference: "<<mufinitederiv<<endl;
  cout<<"Analytic: "<<muanalderiv<<endl;
  cout<<"\n partial w / partial sigma"<<endl;
  cout<<"Finite Difference: "<<sigmafinitederiv<<endl;
  cout<<"Analytic: "<<sigmaanalderiv<<endl;

   cout<<"\n partial2 w / partial mu2"<<endl;
  cout<<"Finite Difference: "<<mufinited2<<endl;
  cout<<"Analytic: "<<md2<<endl;
  cout<<"\n partial2 w / partial sigma2"<<endl;
  cout<<"Finite Difference: "<<sigmafinited2<<endl;
  cout<<"Analytic: "<<sd2<<endl;
  cout<<"\n partial2 w / partial mu sigma"<<endl;
  cout<<"Finite Difference: "<<musigmad<<endl;
  cout<<"Analytic: "<<sigmamuderiv<<endl;
  cout<<"Errors:\n";
  cout<<emud<<endl;
  cout<<esigmad<<endl;
  cout<<emd2<<endl;
  cout<<esd2<<endl;
  cout<<esmd<<endl;

  cout<<"\n\nHessian"<<endl;
  cout<<"\t"<<md2<<"\t"<<smd<<endl;
  cout<<"\t"<<smd<<"\t"<<sd2<<endl;

  if(emud>=etol) throw errtol;
  if(esigmad>=etol) throw errtol;
  if(emd2>=etol) throw errtol;
  if(esd2>=etol) throw errtol;
  if(esmd>=etol) throw errtol;
  
}

void test::RANDOM(double etol){
  grid * gr=_gr;
  //Declare Parameters
  int Nb = gr->numBuses();
  int Nbr = gr->numBranches();
  int samples=5750;
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

  if(sum(abs(fmeanerror))>=etol) throw errtol;
  if(sum(abs(fstdverror))>=etol) throw errtol;
  if(abs(rerror)>=etol/10) throw errtol;
  
}




void test::RV(double etol){
  ranvar rv;
  cout<<"Max error in gaussian CDF function"<<endl;
  if(abs(rv.testPhi())>=etol) throw etol;

  cout<<"Run simulation and compare to analytic"<<endl;
  cout<<"\nSimulation Prob\n"<<endl;
  int samples=73500;
  double mean=.9;
  double variance=.1;
  double stdv=sqrt(variance);
  
  cout<<"Samples: "<<samples<<endl;
  cout<<"Mean: "<<mean<<endl;
  cout<<"Var: "<<variance<<endl;
  

  rv.createRV(samples,mean,stdv);

  double L=.9;
  double p=.15;
  double pc=.5;
  double totalProb=rv.simProb(L,p,pc);
  cout<<"Sim Prob: "<<totalProb<<endl;

  //Analytic

  double t_p = rv.anaProb(L,p,pc);
  double totalerror=t_p-totalProb;
  cout<<"Analytic Prob: "<<t_p<<endl;
  cout<<"\nTotal Error: "<<totalerror<<endl;
  if(abs(totalerror)>=.01) throw errtol;

}

void test::RSLACK(double etol){
  gridcalc gc(_gr);

  int Nb=_gr->numBuses();
  int Nl=_gr->numBranches();
  int Ng = _gr->numGens();


  mat del(Nb,1,fill::zeros);
  del(4,0)=5;
  del(6,0)=7;
  

  srand(time(NULL));
  mat slackdist(Nb,1,fill::zeros);


  double total=0;
  for(int i=0;i<Ng;i++){
    gen g = _gr->getGen(i);
    int buscon=_gr->getBusNum(g.getBus());
    double r = ((double) rand() / (RAND_MAX));
    
    if(i==1) r=1;
    else r=0; 

    slackdist(buscon,0)=r;
    
    total=total+r;
  }


  for(int i=0;i<Ng;i++){
    gen g = _gr->getGen(i);
    int buscon=_gr->getBusNum(g.getBus());
    double normalized = slackdist(buscon,0)/total;
    slackdist(buscon,0)=normalized;
  }    
  
  cout<<"Distribute slack from base bus to slack distribution"<<endl;
  slackdist.t().print("Slack: ");
  mat H=gc.getH();
  mat Hw(H);
  mat vt(Nl,1);
  vt=H*slackdist;
  vt.t().print("Vt: ");
  cout<<"Vt: Cols "<<vt.n_cols<<", Rows "<<vt.n_rows<<endl;
  for(int k=0;k<Nb;k++){
    for(int i=0;i<Nl;i++){
      if(abs(vt(i,0))>=0.0000005){
	Hw(i,k)=H(i,k)-vt(i,0);
      }
    }
  }

  vec del_f = Hw*del;

  del_f.t().print("delta f (shift factor): ");

  cout<<sum(del_f)<<endl;
  
  vec a,b,c;
  a=12*sum(vt)+sum(del_f);
  b=(5*sum(H.col(4))+7*sum(H.col(6)));
  c=(5*sum(H.col(4) - vt) + 7*sum(H.col(6)-vt))+12*sum(vt);

  cout<<a<<endl;
  cout<<b<<endl;
  cout<<c<<endl;

  if(abs(a[0]-b[0])>=etol) throw errtol;
  if(abs(b[0]-c[0])>=etol) throw errtol;
  if(abs(c[0]-a[0])>=etol) throw errtol;

}

void test::PTDF(double etol){
  grid * gr=_gr;
  gridcalc gc(gr);

  int Nb=gr->numBuses();
  vec delg(Nb,fill::zeros);
  vec slackdist(Nb,fill::zeros);
  delg(4)=5; delg(6)=7;
  slackdist(84)=1;

  vec del_f = gc.getDelF(delg,slackdist);
  
  del_g dg(gr);
  dg.addDemand(4,5);
  dg.addDemand(6,7);

  igrid ig(gr);
  ig.addCost();
  rgrid * rg;
  rgrid * rg2;

  double gen_b;
  double gen_a;



  rg = ig.solveModel();
  IloNumArray g_nom=rg->getG();
  gen_b=IloSum(g_nom);
  cout<<"Before"<<endl;
  cout<<"F: "<<rg->getF()<<endl;
  cout<<"G: "<<g_nom<<endl;
  
  ig.modGrid(dg);
    
  IloNumArray slack(IloEnv(),Nb);
  for(int i=0;i<Nb;i++) slack[i]=0;
  slack[5]=1;


  ig.addSlack(g_nom,slack);

  rg2 = ig.solveModel();
  rg2->outputInfo(cerr);
  cout<<"After"<<endl;
  cout<<"F: "<<rg2->getF()<<endl;
  cout<<"G: "<<rg2->getG()<<endl;

  cout<<"Total Demand: "<<gr->getTotalDemand()<<endl;
  cout<<"Total Gen: "<<IloSum(rg2->getG())<<endl;
  gen_a=IloSum(rg2->getG());
  
  cout<<"delF (sim): "<<endl;
  for(int i=0;i<gr->numBranches();i++)
    cout<<" "<<(rg->getF()[i] - rg2->getF()[i])<<"   ";

  cout<<"\ndelF (sim) - delF (shift factor)"<<endl;

  cout<<"\ngen (before) - gen (after) = error"<<endl;
  cout<<gen_b<<" - "<<gen_a<<" = "<<gen_b-gen_a<<endl;

  double totalerror=0;
  for(int i=0;i<gr->numBranches();i++){
    cout<<" "<<(rg->getF()[i] - rg2->getF()[i] - del_f(i))<<"   ";
    totalerror=totalerror+(rg->getF()[i] - rg2->getF()[i] - del_f(i));
  }
  cout<<"\n";
  cout<<"Total Error: "<<totalerror<<endl;
  if(abs(totalerror)>=etol) throw errtol;

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
      fail=false;
      for(int j=0;j<Nl;j++){
	double U = gr->getBranch(j).getRateA();
	if(fcalc[j]>=U)fail=true;  //need redispatch
	if(fcalc[j]<=-U)fail=true; //need redispatch
	}
      if(!fail){
	_rg[i]->outputInfo(cout);
	vec f = gc.convert(_rg[i]->getF());
	fcalc.t().print("fcalc: ");
	f.t().print("f: ");
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
    else{
      cout<<"Line: "<<i<<" --- HB = 1"<<endl;
      
      branch br = gr->getBranch(i);
      int from = gr->getBusNum(br.getFrom());
      int to = gr->getBusNum(br.getTo());
      int winner =-1;

      cout<<from<<" -> "<<to<<endl;

      mat C = gc.getC();
      C.col(from).t().print("from");
      C.col(to).t().print("to");

      if(sum(abs(C.col(from))) == 1) {
	cout<<"winner winner: "<<from<<endl;
	winner=from;
      }
      if(sum(abs(C.col(to))) == 1) {
	cout<<"winner winner: "<<to<<endl;
	winner=to;
      }
      bus b = _gr->getBus(winner);
      double load = b.getP()+b.getGs();
      mat Cm = gc.getCm();
      Cm.row(winner).print("Cm: ");
      int Ng=_gr->numGens();
      if(sum(Cm.row(winner))>=1){
	cout<<"Have Generation"<<endl;
	vec g = gc.convert(rbase->getG());
	vec gnode = Cm.row(winner)*g;
	vec fcalc = Hw.col(winner)*(load-gnode) + fbase;
	for(int j=0;j<Ng;j++){
	  if(Cm(winner,j)>0) {
	    ig.getSlack()->getGNom()[j]=0;
	    ig.slackMismatch();
	  }
	}
	rg = ig.solveModel();
	vec f = gc.convert(rg->getF());
	f.t().print("f actual: ");
	fcalc.t().print("f calc: ");
	(f-fcalc).t().print("error: ");
	if (sum(f-fcalc) > etol ) throw errtol;
	vec ge = gc.convert(rg->getG());
	ge.t().print("gen:");
	for(int j=0;j<Ng;j++){
	  if(Cm(winner,j)>0){
	    ig.getSlack()->getGNom()[j]=g[j];
	    ig.slackMismatch();
	  }
	}
      }
      else{
	cout<<"Loadshed "<<load<<" at "<<winner<<endl;
	vec fcalc = Hw.col(winner)*load+fbase;
	
	dg.setStatus(i,false);
	dg.addDemand(winner,-load);
	ig.modGrid( dg );
	rg = ig.solveModel();
	rg->outputInfo(cout);
	ig.unmodGrid( dg );
	dg.addDemand(winner,0);
	dg.setStatus(i,true);      
	
	vec f = gc.convert(rg->getF());
	f.t().print("f actual: ");
	fcalc.t().print("f calc: ");
	(f-fcalc).t().print("error: ");
	if (sum(f-fcalc) > etol ) throw errtol;
      }
      

    }
  }
  //Solve N-1 problem
  
  
}


