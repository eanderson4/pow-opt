#include <stdlib.h>

#include "sqlinter.h"
#include "test.h"
#include "isj.h"
#include "isjn.h"
#include "ijn1.h"
#include "in1.h"

using namespace std;

int main(int argc, char* argv[]){

  if(argc<=11){
    cout<<"cmd: pow case/30.db <m0> <m1> <e0> <e1> <L> <p> <B> <T>\n"
	<<"\trun main for case30\n"
	<<"\t<m0> base capacity multiplier\n"
	<<"\t<m1> contingency capacity multiplier\n"
	<<"\t<e0> base risk\n"
	<<"\t<e1> contingency risk\n"
	<<"\t<eg> generator risk\n"
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
  double epsG=atof(argv[7]);
  double L=atof(argv[8]);
  double p=atof(argv[9]);
  double B=atof(argv[10]);
  double pc=.85;
  int T = atoi(argv[11]);

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

  int num=5;
  int Nm=num;
  int index[5] = { 1, 3, 6, 20, 28 };
  vec indexM(Nm,fill::zeros);
  double stdv[5] = {3,.5,5,3.75,.3};

  int Nb = gr->numBuses();
  mat Cm(Nb,num,fill::zeros);
  vec mu(num,fill::zeros);
  mat SIG(num,num,fill::zeros);
  for(int i=0;i<num;i++){
    Cm(index[i],i)=1;
    indexM(i)=index[i];
    SIG(i,i)=pow(stdv[i],2)*B;
  }

  SIG(0,2)=-14*B;SIG(2,0)=-14*B;
  SIG(0,3)=-9*B;SIG(3,0)=-9*B;
  SIG(2,3)=15*B;SIG(3,2)=15*B;
  
  double TV = accu(SIG);

  //Generate second demand scenario
  arma_rng::set_seed_random();
  mat as = chol(SIG);


        cout<<*gr<<endl;

  arma_rng::set_seed_random(); 
  gridcalc gc(gr);  
    
  //  vec slack=gc.getSlack()
  mat Cg=gc.getCm();
  vec alpha(Ng,fill::randu);
  vec yes(Ng,fill::randu);
  for(int i=0;i<Ng;i++){
    if(yes(i) <= .5) alpha(i)=0;
  }
  double total=accu(alpha);//  alpha(sn)=1;
  if(total==0) {
    for(int i=0;i<Ng;i++){
      alpha(i)=1;
      total=total+1;
    }
  }
  alpha=alpha/total;

  /*  alpha(0)=atof(argv[2]);
  alpha(1)=atof(argv[3]);
  alpha(2)=atof(argv[4]);
  alpha(3)=atof(argv[5]);
  alpha(4)=atof(argv[6]);
  alpha(5)=1-atof(argv[2]);*/

  for(int i=0;i<Ng;i++){
    alpha(i)=1/Ng;
  }
  alpha = vec(Ng,fill::ones)* (1/(double)Ng);
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
  

  try{

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

  //Set Probability info
  ranvar rv;


  igrid ig(gr);
  ig.addCost();
  isolve is;
  is.setSolver(IloCplex::Dual,IloCplex::Dual);
  rgrid * rbase = ig.solveModel(&is);
  //  rbase->displayOperatingPos(gr);

  /*  mat xnom(Ng,2);
  mat xsjn(Ng,2);
  mat betanom(Ng,2);
  mat betasjn(Ng,2);
  mat cnom(1,2);
  mat csjn(1,2);
  mat risknom(1,2);
  mat risksjn(1,2);
  */
  double TG;
 

  
    running_stat<double> costs_r0;
    running_stat<double> costs_r1;
    running_stat<double> costs_r2;
    running_stat<double> costs_r3;
    running_stat<double> costs_r4;
    running_stat<double> costs_r5;  

    running_stat<double> risk_r0;
    running_stat<double> risk_r1;
    running_stat<double> risk_r2;
    running_stat<double> risk_r3;
    running_stat<double> risk_r4;
    running_stat<double> risk_r5;  


    //    running_stat_vec<vec> ct(true);
    mat cost_store(T,6);
    vec delta_store(T);

    vec gen_opf;

    vec beta_main;
    vec gen_main;
    running_stat_vec<vec> bt_stat;
    running_stat_vec<vec> gen_dev_opf;
    running_stat_vec<vec> gen_dev;
    running_stat_vec<vec> gen_2nd;
    running_stat_vec<vec> gen_exp;

    for(int trial=0;trial<T;trial++){
      
  vec randd(Nm,fill::randu);
  vec d2;
  if(trial==0) d2=vec(Nm,fill::zeros);
  else d2 = as*randd;

  for(int m=0;m<Nm;m++){
    gr->addPd(indexM(m),d2(m));
  }      
  double delta=accu(d2);
  delta_store(trial) = accu(d2);

  vec costs(6,fill::zeros);

    ijn1 n1(gr, SIGy,Hw,L,p,pc,eps,epsN);
    n1.addCost();

  try{
    igrid nom(gr);
    nom.addCost();
    rgrid * rnom = nom.solveModel(&is);

    double o0=rnom->getObjective();
    vec f0=gc.convert(rnom->getF());
    vec g0=gc.convert(rnom->getG());
    vec z0=gc.risk(f0,SIGy.diag(),L,p,pc);
    double r0 = sum(z0);
    IloCplex::CplexStatus s0=rnom->getStatus();
    TG = accu(g0);

    running_stat<double> stats_r0;
    vec check = n1.getCheck();
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

    costs_r0(o0 + randcost);
    risk_r0(r0);
    costs(0) = o0;

    if(trial==0) gen_opf = g0;
    else{
      vec gen_expect = gen_opf + alpha*delta;
      gen_dev_opf( gen_expect - g0 );

    }

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

    in1 bn1(gr, SIGy, Hw, m1); 
    bn1.addCost();
    rgrid * rbn1 = bn1.solveModel(&is);

    double o3=rbn1->getObjective();
    vec f3=gc.convert(rbn1->getF());
    vec g3=gc.convert(rbn1->getG());
    vec z3=gc.risk(f3,SIGy.diag(),L,p,pc);
    double r3 = sum(z3);
    IloCplex::CplexStatus s3=rbn1->getStatus();

    running_stat<double> stats_r3;
    vec check = n1.getCheck();
    for(int i=0;i<Nl;i++){
      if(check(i)==1){
	vec f3n = n1.getN1(i,f3,g3);
	vec z3n=gc.risk(f3n,SIGy.diag(),L,p,pc);
	double r3n = sum(z3n);
	stats_r3(r3n);
      }
    }

    cout<<"OPF N-1"<<"\t"<<s3<<endl;
    cout<<"C3: "<<o3<<endl;
    cout<<"r3 - "<<r3<<endl;
    cout << "count = " << stats_r3.count() << endl;
    cout << "mean = " << stats_r3.mean() << endl;
    cout << "stdv  = " << stats_r3.stddev()  << endl;
    cout << "min  = " << stats_r3.min()  << endl;
    cout << "max  = " << stats_r3.max()  << endl;
    cout<<endl;

    costs_r3(o3 + randcost);
    risk_r3(r3);

    costs(3) = o3;
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

    //    ijn1 n1(gr, SIGy,Hw,L,p,pc,eps,epsN);
    //    n1.addCost();
    rgrid * rn1_1 = n1.solveModel(&is);

    double o1=rn1_1->getObjective();
    vec f1=gc.convert(rn1_1->getF());
    vec g1=gc.convert(rn1_1->getG());
    vec z1=gc.risk(f1,SIGy.diag(),L,p,pc);
    double r1 = sum(z1);
    IloCplex::CplexStatus s1=rn1_1->getStatus();

    running_stat<double> stats_r1;
    vec check = n1.getCheck();
    for(int i=0;i<Nl;i++){
      if(check(i)==1){
	vec f1n = n1.getN1(i,f1,g1);
	vec z1n=gc.risk(f1n,SIGy.diag(),L,p,pc);
	double r1n = sum(z1n);
	stats_r1(r1n);
      }
    }

    cout<<"JCC N-1"<<"\t"<<s1<<endl;
    cout<<"C1: "<<o1<<endl;
    cout<<"r1 - "<<r1<<endl;
    cout << "count = " << stats_r1.count() << endl;
    cout << "mean = " << stats_r1.mean() << endl;
    cout << "stdv  = " << stats_r1.stddev()  << endl;
    cout << "min  = " << stats_r1.min()  << endl;
    cout << "max  = " << stats_r1.max()  << endl;
    cout<<endl;

    costs_r1(o1 + randcost);
    risk_r1(r1);

    costs(1) = o1;


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

    ijn1 jcc(gr, SIGy,Hw,L,p,pc,eps,1);
    jcc.addCost();
    rgrid * rjcc = jcc.solveModel(&is);

    double o2=rjcc->getObjective();
    vec f2=gc.convert(rjcc->getF());
    vec g2=gc.convert(rjcc->getG());
    vec z2=gc.risk(f2,SIGy.diag(),L,p,pc);
    double r2 = sum(z2);
    IloCplex::CplexStatus s2=rjcc->getStatus();

    running_stat<double> stats_r2;
    vec check = n1.getCheck();
    for(int i=0;i<Nl;i++){
      if(check(i)==1){
	vec f2n = n1.getN1(i,f2,g2);
	vec z2n=gc.risk(f2n,SIGy.diag(),L,p,pc);
	double r2n = sum(z2n);
	stats_r2(r2n);
      }
    }

    cout<<"JCC"<<"\t"<<s2<<endl;
    cout<<"C2: "<<o2<<endl;
    cout<<"r2 - "<<r2<<endl;
    cout << "count = " << stats_r2.count() << endl;
    cout << "mean = " << stats_r2.mean() << endl;
    cout << "stdv  = " << stats_r2.stddev()  << endl;
    cout << "min  = " << stats_r2.min()  << endl;
    cout << "max  = " << stats_r2.max()  << endl;
    cout<<endl;    

    costs_r2(o2 + randcost);
    risk_r2(r2);

    costs(2) = o2;

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

    isj sj(gr, &gc, SIG, indexM, L, p, pc, eps);
    rgrid * rsj = sj.solveModel(&is);

    double o4=rsj->getObjective();
    vec f4=gc.convert(rsj->getF());
    vec g4=gc.convert(rsj->getG());
    vec beta4=sj.getBeta();
    vec sd4=sj.getSD();
    vec z4=gc.risk(f4,sd4,L,p,pc);
    double r4 = sum(z4);
    IloCplex::CplexStatus s4=rsj->getStatus();

    running_stat<double> stats_r4;
    vec check = n1.getCheck();
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

    costs_r4(o4);
    risk_r4(r4);

        costs(4) = o4;

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

    isjn sjn(gr, &gc, SIG, indexM, L, p, pc, eps,eN,epsG);
    rgrid * rsjn = sjn.solveModel(&is);

    double o5=rsjn->getObjective();
    vec f5=gc.convert(rsjn->getF());
    vec g5=gc.convert(rsjn->getG());
    vec beta5=sjn.getBeta();
    vec sd5=sjn.getSD();
    vec z5=gc.risk(f5,sd5,L,p,pc);
    double r5 = sum(z5);
    IloCplex::CplexStatus s5=rsjn->getStatus();

    running_stat<double> stats_r5;
    vec check = n1.getCheck();
    for(int i=0;i<Nl;i++){
      if(check(i)==1){
	vec f5n = n1.getN1(i,f5,g5);
	vec z5n=gc.risk(f5n,sd5,L,p,pc);
	double r5n = sum(z5n);
	stats_r5(r5n);
      }
    }

    cout<<"SJ N-1"<<"\t"<<s5<<endl;
    cout<<"C5: "<<o5<<endl;
    cout<<"r5 - "<<r5<<endl;
    cout << "count = " << stats_r5.count() << endl;
    cout << "mean = " << stats_r5.mean() << endl;
    cout << "stdv  = " << stats_r5.stddev()  << endl;
    cout << "min  = " << stats_r5.min()  << endl;
    cout << "max  = " << stats_r5.max()  << endl;
    cout<<endl;

    costs_r5(o5);
    risk_r5(r5);
    
    costs(5) = o5;

    bt_stat(beta5);

    if(trial==0){
      beta_main=beta5;
      gen_main=g5;
    }
    else {
      vec gen_expect = gen_main + beta_main*delta;

      gen_dev( gen_expect-gen_main );
      gen_2nd( g5 );
      gen_exp( gen_expect);
      
    }

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

  
    /*
    cerr<<"\t"<<randcost<<"\t"<<TV<<"\t\t";
    cerr<<o0<<"\t"<<r0<<"\t"<<stats_r0.mean()<<"\t"<<stats_r0.max()<<"\t\t";
    cerr<<o2<<"\t"<<r2<<"\t"<<stats_r2.mean()<<"\t"<<stats_r2.max()<<"\t\t";
    cerr<<o3<<"\t"<<r3<<"\t"<<stats_r3.mean()<<"\t"<<stats_r3.max()<<"\t\t";
    cerr<<o1<<"\t"<<r1<<"\t"<<stats_r1.mean()<<"\t"<<stats_r1.max()<<endl;*/
  //  ct(costs);
  cost_store.row(trial) = costs.t();

  for(int m=0;m<Nm;m++){
    gr->addPd(indexM(m),-d2(m));
  }      

  }





  cout<<"2nd stage Comparison"<<endl;


    cout<<"Total Gen: "<<TG<<endl;
    cout<<"Total Variance: "<<TV<<" ( "<<sqrt(TV)<<" )"<<endl;
    cout<<"Total Random Cost: "<<randcost<<endl;
    cout<<"G inv: "<<rv.ginv(eps,L,p,pc)<<endl;
    SIG.print("Cov m: ");
    alpha.t().print("slack: ");
    cout<<"\n\n";
    cout<<"OPF"<<endl;
    cout << "costs: "<<endl;
    cout << "count = " << costs_r0.count() << endl;
    cout << "mean = " << costs_r0.mean() << endl;
    cout << "stdv  = " << costs_r0.stddev()  << endl;
    cout << "min  = " << costs_r0.min()  << endl;
    cout << "max  = " << costs_r0.max()  << endl;
    cout << "band = " << costs_r0.max() - costs_r0.min() <<endl;
    cout<<endl;
    cout<<"JCC"<<endl;
    cout << "costs: "<<endl;
    cout << "count = " << costs_r2.count() << endl;
    cout << "mean = " << costs_r2.mean() << endl;
    cout << "stdv  = " << costs_r2.stddev()  << endl;
    cout << "min  = " << costs_r2.min()  << endl;
    cout << "max  = " << costs_r2.max()  << endl;
    cout << "band = " << costs_r2.max() - costs_r2.min() <<endl;
    cout<<endl;    
    cout<<"SJ"<<endl;
    cout << "costs: "<<endl;
    cout << "count = " << costs_r4.count() << endl;
    cout << "mean = " << costs_r4.mean() << endl;
    cout << "stdv  = " << costs_r4.stddev()  << endl;
    cout << "min  = " << costs_r4.min()  << endl;
    cout << "max  = " << costs_r4.max()  << endl;
    cout << "band = " << costs_r4.max() - costs_r4.min() <<endl;
    cout<<endl;
    cout<<"OPF N-1"<<endl;
    cout << "costs: "<<endl;
    cout << "count = " << costs_r3.count() << endl;
    cout << "mean = " << costs_r3.mean() << endl;
    cout << "stdv  = " << costs_r3.stddev()  << endl;
    cout << "min  = " << costs_r3.min()  << endl;
    cout << "max  = " << costs_r3.max()  << endl;
    cout << "band = " << costs_r3.max() - costs_r3.min() <<endl;
    cout<<endl;
    cout<<"JCC N-1"<<endl;
    cout << "costs: "<<endl;
    cout << "count = " << costs_r1.count() << endl;
    cout << "mean = " << costs_r1.mean() << endl;
    cout << "stdv  = " << costs_r1.stddev()  << endl;
    cout << "min  = " << costs_r1.min()  << endl;
    cout << "max  = " << costs_r1.max()  << endl;
    cout << "band = " << costs_r1.max() - costs_r1.min() <<endl;
    cout<<endl;
    cout<<"SJ N-1"<<endl;
    cout << "costs: "<<endl;
    cout << "count = " << costs_r5.count() << endl;
    cout << "mean = " << costs_r5.mean() << endl;
    cout << "stdv  = " << costs_r5.stddev()  << endl;
    cout << "min  = " << costs_r5.min()  << endl;
    cout << "max  = " << costs_r5.max()  << endl;
    cout << "band = " << costs_r5.max() - costs_r5.min() <<endl;
    cout<<endl;


    cout<<"\n\n";
    cout<<"OPF"<<endl;
    cout << "risk: "<<endl;
    cout << "count = " << risk_r0.count() << endl;
    cout << "mean = " << risk_r0.mean() << endl;
    cout << "stdv  = " << risk_r0.stddev()  << endl;
    cout << "min  = " << risk_r0.min()  << endl;
    cout << "max  = " << risk_r0.max()  << endl;
    cout<<endl;
    cout<<"JCC"<<endl;
    cout << "risk: "<<endl;
    cout << "count = " << risk_r2.count() << endl;
    cout << "mean = " << risk_r2.mean() << endl;
    cout << "stdv  = " << risk_r2.stddev()  << endl;
    cout << "min  = " << risk_r2.min()  << endl;
    cout << "max  = " << risk_r2.max()  << endl;
    cout<<endl;    
    cout<<"SJ"<<endl;
    cout << "risk: "<<endl;
    cout << "count = " << risk_r4.count() << endl;
    cout << "mean = " << risk_r4.mean() << endl;
    cout << "stdv  = " << risk_r4.stddev()  << endl;
    cout << "min  = " << risk_r4.min()  << endl;
    cout << "max  = " << risk_r4.max()  << endl;
    cout<<endl;
    cout<<"OPF N-1"<<endl;
    cout << "risk: "<<endl;
    cout << "count = " << risk_r3.count() << endl;
    cout << "mean = " << risk_r3.mean() << endl;
    cout << "stdv  = " << risk_r3.stddev()  << endl;
    cout << "min  = " << risk_r3.min()  << endl;
    cout << "max  = " << risk_r3.max()  << endl;
    cout<<endl;
    cout<<"JCC N-1"<<endl;
    cout << "risk: "<<endl;
    cout << "count = " << risk_r1.count() << endl;
    cout << "mean = " << risk_r1.mean() << endl;
    cout << "stdv  = " << risk_r1.stddev()  << endl;
    cout << "min  = " << risk_r1.min()  << endl;
    cout << "max  = " << risk_r1.max()  << endl;
    cout<<endl;
    cout<<"SJ N-1"<<endl;
    cout << "risk: "<<endl;
    cout << "count = " << risk_r5.count() << endl;
    cout << "mean = " << risk_r5.mean() << endl;
    cout << "stdv  = " << risk_r5.stddev()  << endl;
    cout << "min  = " << risk_r5.min()  << endl;
    cout << "max  = " << risk_r5.max()  << endl;
    cout<<endl;


    bt_stat.mean().t().print("bt mean: ");
    bt_stat.stddev().t().print("bt stdev: ");


    gen_dev_opf.mean().t().print("gen dev opf mean: ");
    gen_dev_opf.stddev().t().print("gen dev opf stdv: ");


    gen_dev.mean().t().print("gen dev mean: ");
    gen_dev.stddev().t().print("gen dev stdv: ");


    gen_2nd.mean().t().print("gen 2nd mean: ");
    gen_2nd.stddev().t().print("gen 2nd stdv: ");


    gen_exp.mean().t().print("gen exp mean: ");
    gen_exp.stddev().t().print("gen exp stdv: ");


    /*    ct.mean().t().print("mean: ");
    ct.cov().print("cov: ");
    mat CV = ct.cov();
    mat COR = CV;
    for(int i=0; i<6;i++){
      for(int j=0; j<6; j++){
	double sdi = sqrt(CV(i,i));
	double sdj = sqrt(CV(j,j));
	COR(i,j)=CV(i,j)/sdi/sdj;
      }
    }

    COR.print("cor: ");
    */
    mat sort_costs = sort(cost_store);

    sort_costs.print("costs: ");

    for(int i=0; i<T; i++){
      //      cerr<<i<<"\t"<<sort_costs(i,0)<<"\t"<<sort_costs(i,5)<<endl;
      cerr<<i<<"\t"<<delta_store(i)<<"\t"<<cost_store(i,0)<<"\t"<<cost_store(i,5)<<endl;
    }

    /*    
    if(cnom(0)>1 && cnom(1)>1 && csjn(0)>1 && csjn(1)>1){
      cerr<<cnom(0)<<"\t"<<cnom(1)<<"\t"<<csjn(0)<<"\t"<<csjn(1)<<"\t"<<difnom<<"\t"<<difsjn<<"\t"<<difsjn-difnom<<"\t";
      cerr<<risknom(0)<<"\t"<<risknom(1)<<"\t"<<risksjn(0)<<"\t"<<risksjn(1)<<endl;
    }

  */

  return 0;



}

   
