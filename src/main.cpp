#include <stdlib.h>

#include <ctime>
#include "sqlinter.h"

#include "icc.h"
#include "isj.h"
#include "isjn.h"
#include "ijn1.h"
#include "in1.h"

using namespace std;

int main(int argc, char* argv[]){

  if(argc<=12){
    cout<<"cmd: pow case/30.db <m0> <m1> <e0> <e1> <L> <p> <B> <T>\n"
	<<"\trun main for case30\n"
	<<"\t<m0> base capacity multiplier\n"
	<<"\t<m1> contingency capacity multiplier\n"
	<<"\t<e0> base risk\n"
	<<"\t<e1> contingency risk\n"
	<<"\t<pl> line fail prob\n"
	<<"\t<pg> generator fail prob\n"
	<<"\t<L> no risk intercept\n"
	<<"\t<p> probability of failure at nominal\n"
	<<"\t<B> Variance budget\n"
	<<"\t<T> num trials"<<endl;
    return 1;
  }

  double m0=atof(argv[2]);
  double m1=atof(argv[3]);
  double mg=atof(argv[4]);
  double eps=atof(argv[5]);
  double epsN=atof(argv[6]);
  double epsL=atof(argv[7]);
  double epsG=atof(argv[8]);
  double L=atof(argv[9]);
  double p=atof(argv[10]);
  double B=atof(argv[11]);
  double pc=.85;
  int T = atoi(argv[12]);

  ///  int sn=atoi(argv[7]); //standard deviation test

  sqlInter db;
  grid * gr = new grid;

  string db_name;
  
  db_name = argv[1];
  db.openDb(db_name);
  db.load(*gr);
  gr->buildMap();
  gr->printNums(cout);

  cout<<*gr<<endl;
  


  int Nl = gr->numBranches();
  for(int i=0;i<Nl;i++){
    branch bi = gr->getBranch(i);
    double U = bi.getRateA();
    gr->setCapacity(i,m0*U);
  }

  int Ng = gr->numGens();
  for(int j=0;j<Ng;j++){
    gen gj = gr->getGen(j);
    double U = gj.getPmax();
    gr->setGenPmax(j,mg*U);
  } 

  int Nb = gr->numBuses();
  vector<int> nodes;
  vector<int> demand;
  for(int i=0;i<Nb;i++){
    bus bn = gr->getBus(i);
    double p = bn.getP();
    
    if(p>1){ 
      nodes.push_back(i);
      demand.push_back(p);
    }  
  }

  int Nm=nodes.size();
  vec indexM(Nm,fill::zeros);


  mat Cm(Nb,Nm,fill::zeros);
  vec mu(Nm,fill::zeros);
  mat SIG(Nm,Nm,fill::zeros);
  for(int i=0;i<Nm;i++){
    Cm(nodes[i],i)=1;
    indexM(i)=nodes[i];
    SIG(i,i)=pow(.05*demand[i],2)*B;
  }
  SIG.diag().t().print("StDv Volatile Injects: ");
  double TV = accu(SIG);

  //Generate second demand scenario
  arma_rng::set_seed_random();
  mat as = chol(SIG);

  cout<<*gr<<endl;

  arma_rng::set_seed_random(); 
  gridcalc gc(gr);  
    

  mat Cg=gc.getCm();
  vec alpha(Ng,fill::zeros);

  for(int i=0;i<Ng;i++){
    gen gi=gr->getGen(i);
    double status = gi.getStatus();
    if (status>0)   alpha(i)=1;
  }
  int numSlack = accu(alpha);
  alpha = alpha/numSlack;
  alpha.t().print("alpha: ");

  vec slack=Cg*alpha;

  double randcost=0;
  for(int i=0;i<Ng;i++){
    double c2=gr->getGen(i).getC2();
    double genrandcost=c2*pow(alpha(i),2)*TV;
    randcost=randcost+genrandcost;
  }

  mat Hw = gc.getHw(slack);
  mat A = gc.getH();
  mat Lo = gc.getL(A);  


  mat SIGy(Nl,Nl,fill::zeros);
  SIGy = Hw*Cm*SIG*(Hw*Cm).t();
  vec sd = SIGy.diag();
  


  vec eN(Nl, fill::ones);
  eN = eN*epsN;
  

  //Set Probability info
  ranvar rv;


  igrid ig(gr);
  ig.addCost();
  isolve is;
  is.setSolver(IloCplex::Dual,IloCplex::Dual);
  //  rgrid * rbase = ig.solveModel(&is);

  double TG;
  
  
  
  running_stat<double> costs_opf;
  running_stat<double> costs_cc;
  running_stat<double> costs_sjc;

  running_stat<double> risk_opf;
  running_stat<double> risk_cc;
  running_stat<double> risk_sjc;

  running_stat<double> prob_mean_opf;
  running_stat<double> prob_mean_cc;
  running_stat<double> prob_mean_sjc;

  running_stat<double> prob_max_opf;
  running_stat<double> prob_max_cc;
  running_stat<double> prob_max_sjc;

  running_stat<double> time_opf;
  running_stat<double> time_cc;
  running_stat<double> time_sjc;
  running_stat<double> formtime_cc;
  running_stat<double> formtime_sjc;
  

  vec delta_store(T);


  ijn1 n1(gr,  SIGy,Hw,L,p,pc,eps,epsN);
  n1.addCost();
  vec check = n1.getCheck();
      
    for(int trial=0;trial<T;trial++){
      
      vec randd(Nm,fill::randu);
      vec d2;
      if(trial==0) d2=vec(Nm,fill::zeros);
      else d2 = as*randd;
      
      for(int m=0;m<Nm;m++){
	gr->addPd(indexM(m),d2(m));
      }      
      double delta=accu(d2);
      delta_store(trial) = accu(delta);
            
      



      try{
	clock_t start = clock();
	igrid nom(gr);
	nom.addCost();
	rgrid * rnom = nom.solveModel(&is);
	double time = (clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
	

	
	double o0=rnom->getObjective();
	vec f0=gc.convert(rnom->getF());
	vec g0=gc.convert(rnom->getG());
	vec z0=gc.risk(f0,SIGy.diag(),L,p,pc);
	double r0 = sum(z0);
	vec p0=gc.lineprob(f0,SIGy.diag());
	IloCplex::CplexStatus s0=rnom->getStatus();
	TG = accu(g0);
	
	running_stat<double> stats_r0;
	for(int i=0;i<Nl;i++){
	  if(check(i)==1){
	    vec f0n = n1.getN1(i,f0,g0);
	    vec z0n=gc.risk(f0n,SIGy.diag(),L,p,pc);
	    double r0n = sum(z0n);
	    stats_r0(r0n);
	  }
	}
	
	cout<<"OPF"<<"\t"<<s0<<endl;
	cout<<"C0: "<<o0<<endl;
	cout<<"r0 - "<<r0<<endl;
	cout << "count = " << stats_r0.count() << endl;
	cout << "mean = " << stats_r0.mean() << endl;
	cout << "stdv  = " << stats_r0.stddev()  << endl;
	cout << "min  = " << stats_r0.min()  << endl;
	cout << "max  = " << stats_r0.max()  << endl;
	cout<<endl;
	
	costs_opf(o0 + randcost);
	risk_opf(r0);
	time_opf(time);
	prob_mean_opf(mean(p0));
	prob_max_opf(max(p0));

	
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
      
      try{
	clock_t start = clock();	
	icc cc(gr, &gc, SIG, indexM, epsL, epsG);
	double formtime = (clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
	rgrid * rcc = cc.solveModel(&is);
	double time = (clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
	
	double o=rcc->getObjective();
	vec f=gc.convert(rcc->getF());
	vec g=gc.convert(rcc->getG());
	vec beta=cc.getBeta();
	vec sd=cc.getSD();
	vec z=gc.risk(f,sd,L,p,pc);
	double r = sum(z);
	vec p1=gc.lineprob(f,sd);
	IloCplex::CplexStatus s=rcc->getStatus();
	
	running_stat<double> stats_risk;
	for(int i=0;i<Nl;i++){
	  if(check(i)==1){
	    vec fn = n1.getN1(i,f,g);
	    vec zn=gc.risk(fn,sd,L,p,pc);
	    double rn = sum(zn);
	    stats_risk(rn);
	  }
	}
	
	cout<<"CC"<<"\t"<<s<<endl;
	cout<<"C: "<<o<<endl;
	cout<<"r - "<<r<<endl;
	cout << "count = " << stats_risk.count() << endl;
	cout << "mean = " << stats_risk.mean() << endl;
	cout << "stdv  = " << stats_risk.stddev()  << endl;
	cout << "min  = " << stats_risk.min()  << endl;
	cout << "max  = " << stats_risk.max()  << endl;
	cout<<endl;
	
	costs_cc(o);
	risk_cc(r);
	prob_mean_cc(mean(p1));
	prob_max_cc(max(p1));
	time_cc(time);	
	formtime_cc(formtime);	
	//	costs(4) = o4;
	
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
      

      try{
	clock_t start = clock();	
	isj sj(gr, &gc, SIG, indexM, L, p, pc, eps);
	double formtime = (clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
	rgrid * rsj = sj.solveModel(&is);
	double time = (clock() - start) / (double)(CLOCKS_PER_SEC / 1000);
	
	double o4=rsj->getObjective();
	vec f4=gc.convert(rsj->getF());
	vec g4=gc.convert(rsj->getG());
	vec beta4=sj.getBeta();
	vec sd4=sj.getSD();
	vec z4=gc.risk(f4,sd4,L,p,pc);
	double r4 = sum(z4);
	vec p4=gc.lineprob(f4,sd4);
	IloCplex::CplexStatus s4=rsj->getStatus();
	
	running_stat<double> stats_r4;
	for(int i=0;i<Nl;i++){
	  if(check(i)==1){
	    vec f4n = n1.getN1(i,f4,g4);
	    vec z4n=gc.risk(f4n,sd4,L,p,pc);
	    double r4n = sum(z4n);
	    stats_r4(r4n);
	  }
	}
	
	cout<<"SJ"<<"\t"<<s4<<endl;
	cout<<"C4: "<<o4<<endl;
	cout<<"r4 - "<<r4<<endl;
	cout << "count = " << stats_r4.count() << endl;
	cout << "mean = " << stats_r4.mean() << endl;
	cout << "stdv  = " << stats_r4.stddev()  << endl;
	cout << "min  = " << stats_r4.min()  << endl;
	cout << "max  = " << stats_r4.max()  << endl;
	cout<<endl;
	
	costs_sjc(o4);
	risk_sjc(r4);
	prob_mean_sjc(mean(p4));
	prob_max_sjc(max(p4));
	time_sjc(time);	
	formtime_sjc(formtime);	
	//	costs(4) = o4;
	
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


      for(int m=0;m<Nm;m++){
	gr->addPd(indexM(m),-d2(m));
      }      
      
    }
    
    

    
    
    cout<<"2nd stage Comparison"<<endl;
    //    cout<<*gr<<endl;
    SIG.diag().t().print("Cov.diag m: ");
    alpha.t().print("slack: ");
    

    cout<<"Parameters: "<<endl;
    cout<<"m0\tm1\tmg\teps\tepsN\tpL\tpG\tL\tp\tB\tT"<<endl;
    cout<<m0<<"\t"<<m1<<"\t"<<mg<<"\t"<<eps<<"\t"<<epsN<<"\t"<<epsL<<"\t"<<epsG<<"\t"<<L<<"\t"<<p<<"\t"<<B<<"\t"<<T<<endl;

    cout<<"Total Gen: "<<TG<<endl;
    cout<<"Total Variance: "<<TV<<" ( "<<sqrt(TV)<<" )"<<endl;
    cout<<"Total Random Cost: "<<randcost<<endl;
    cout<<"G inv: "<<rv.ginv(eps,L,p,pc)<<endl;

    cout<<"\n\n";
    cout<<" --- \t --- RISK --- \t ---\n"<<endl;
    cout<<"OPF"<<endl;
    cout << "risk: "<<endl;
    cout << "count = " << risk_opf.count() << endl;
    cout << "mean = " << risk_opf.mean() << endl;
    cout << "stdv  = " << risk_opf.stddev()  << endl;
    cout << "min  = " << risk_opf.min()  << endl;
    cout << "max  = " << risk_opf.max()  << endl;
    cout << "prob: "<<endl;
    cout << "mean  = " << prob_mean_opf.mean()  << endl;
    cout << "max  = " << prob_max_opf.mean()  << endl;
    cout<<endl;

    cout<<"CC"<<endl;
    cout << "risk: "<<endl;
    cout << "count = " << risk_cc.count() << endl;
    cout << "mean = " << risk_cc.mean() << endl;
    cout << "stdv  = " << risk_cc.stddev()  << endl;
    cout << "min  = " << risk_cc.min()  << endl;
    cout << "max  = " << risk_cc.max()  << endl;
    cout << "prob: "<<endl;
    cout << "mean  = " << prob_mean_cc.mean()  << endl;
    cout << "max  = " << prob_max_cc.mean()  << endl;
    cout<<endl;


    cout<<"SJ"<<endl;
    cout << "risk: "<<endl;
    cout << "count = " << risk_sjc.count() << endl;
    cout << "mean = " << risk_sjc.mean() << endl;
    cout << "stdv  = " << risk_sjc.stddev()  << endl;
    cout << "min  = " << risk_sjc.min()  << endl;
    cout << "max  = " << risk_sjc.max()  << endl;
    cout << "prob: "<<endl;
    cout << "mean  = " << prob_mean_sjc.mean()  << endl;
    cout << "max  = " << prob_max_sjc.mean()  << endl;
    cout<<endl;


    cout<<"\n\n";
    cout<<" --- \t --- TIME --- \t ---\n"<<endl;
    cout<<"OPF"<<endl;
    cout << "time: "<<endl;
    cout << "count = " << time_opf.count() << endl;
    cout << "mean = " << time_opf.mean() << endl;
    cout << "stdv  = " << time_opf.stddev()  << endl;
    cout << "min  = " << time_opf.min()  << endl;
    cout << "max  = " << time_opf.max()  << endl;
    cout<<endl;

    cout<<"CC"<<endl;
    cout << "time: "<<endl;
    cout << "count = " << time_cc.count() << endl;
    cout << "mean = " << time_cc.mean() << endl;
    cout << "stdv  = " << time_cc.stddev()  << endl;
    cout << "min  = " << time_cc.min()  << endl;
    cout << "max  = " << time_cc.max()  << endl;
    cout<<endl;

    cout<<"SJ"<<endl;
    cout << "time: "<<endl;
    cout << "count = " << time_sjc.count() << endl;
    cout << "mean = " << time_sjc.mean() << endl;
    cout << "stdv  = " << time_sjc.stddev()  << endl;
    cout << "min  = " << time_sjc.min()  << endl;
    cout << "max  = " << time_sjc.max()  << endl;
    cout<<endl;


    cout<<"SJ"<<endl;
    cout << "formtime: "<<endl;
    cout << "count = " << formtime_sjc.count() << endl;
    cout << "mean = " << formtime_sjc.mean() << endl;
    cout << "stdv  = " << formtime_sjc.stddev()  << endl;
    cout << "min  = " << formtime_sjc.min()  << endl;
    cout << "max  = " << formtime_sjc.max()  << endl;
    cout<<endl;



    cout<<"\n\n";
    cout<<" --- \t --- COST --- \t ---\n"<<endl;
    cout<<"OPF"<<endl;
    cout << "costs: "<<endl;
    cout << "count = " << costs_opf.count() << endl;
    cout << "mean = " << costs_opf.mean() << endl;
    cout << "stdv  = " << costs_opf.stddev()  << endl;
    cout << "min  = " << costs_opf.min()  << endl;
    cout << "max  = " << costs_opf.max()  << endl;
    cout << "band = " << costs_opf.max() - costs_opf.min() <<endl;
    cout<<endl;

    cout<<"CC"<<endl;
    cout << "costs: "<<endl;
    cout << "count = " << costs_cc.count() << endl;
    cout << "mean = " << costs_cc.mean() << endl;
    cout << "stdv  = " << costs_cc.stddev()  << endl;
    cout << "min  = " << costs_cc.min()  << endl;
    cout << "max  = " << costs_cc.max()  << endl;
    cout << "band = " << costs_cc.max() - costs_cc.min() <<endl;
    cout<<endl;


    cout<<"SJ"<<endl;
    cout << "costs: "<<endl;
    cout << "count = " << costs_sjc.count() << endl;
    cout << "mean = " << costs_sjc.mean() << endl;
    cout << "stdv  = " << costs_sjc.stddev()  << endl;
    cout << "min  = " << costs_sjc.min()  << endl;
    cout << "max  = " << costs_sjc.max()  << endl;
    cout << "band = " << costs_sjc.max() - costs_sjc.min() <<endl;
    cout<<endl;







  return 0;



}

   
