#include <stdlib.h>

#include "sqlinter.h"
#include "test.h"
#include "isj.h"
#include "ijn1.h"
#include "in1.h"

using namespace std;

int main(int argc, char* argv[]){

  if(argc<=9){
    cout<<"cmd: pow case/30.db <m0> <m1> <e0> <e1> <L> <p> <B>\n"
	<<"\trun main for case30\n"
	<<"\t<m0> base capacity multiplier\n"
	<<"\t<m1> contingency capacity multiplier\n"
	<<"\t<e0> base risk\n"
	<<"\t<e1> contingency risk\n"
	<<"\t<L> no risk intercept\n"
	<<"\t<p> probability of failure at nominal\n"
	<<"\t<B> Variance budget"<<endl;
    return 1;
  }

  double m0=atof(argv[2]);
  double m1=atof(argv[3]);
  double eps=atof(argv[4]);
  double epsN=atof(argv[5]);
  double L=atof(argv[6]);
  double p=atof(argv[7]);
  double B=atof(argv[8]);
  int sn=atoi(argv[9]);
  double pc=.85;


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

  int num=5;
  int Nm=num;
  int index[5] = { 1, 3, 6, 20, 28 };
  double stdv[5] = {3,.5,5,3.75,.3};

  int Nb = gr->numBuses();
  mat Cm(Nb,num,fill::zeros);
  vec mu(num,fill::zeros);
  mat SIG(num,num,fill::zeros);
  for(int i=0;i<num;i++){
    Cm(index[i],i)=1;
    SIG(i,i)=pow(stdv[i],2)*B;
  }
  SIG(0,2)=-14*B;SIG(2,0)=-14*B;
  SIG(0,3)=-9*B;SIG(3,0)=-9*B;
  SIG(2,3)=15*B;SIG(3,2)=15*B;
  
  double TV = accu(SIG);

  arma_rng::set_seed_random(); 
  gridcalc gc(gr);  
  //  vec slack=gc.getSlack()
  int Ng = gr->numGens();;
  mat Cg=gc.getCm();
  vec alpha(Ng,fill::randu);
  vec yes(Ng,fill::randu);
  for(int i=0;i<Ng;i++){
    if(yes(i) <= .5) alpha(i)=0;
  }
  double total=accu(alpha);//  alpha(sn)=1;
  if(total==0) {alpha(sn)=1;total=1;}
  alpha=alpha/total;
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

  mat Ag(Nl,Ng,fill::zeros);
  for(int e=0;e<Nl;e++){
    for(int j=0;j<Ng;j++){
      Ag(e,j) = dot(A.row(e),Cg.col(j)); 
    }
  }
  //  Ag.print("Ag: ");

  mat Ak(Nl,Nm,fill::zeros);
  for(int e=0;e<Nl;e++){
    for(int k=0;k<Nm;k++){
      Ak(e,k) = dot(A.row(e),Cm.col(k)); 
    }
  }
  //  Ak.print("Ak: ");

  vec sig(Nl,fill::zeros);
  mat sigger(Nl,Nl,fill::zeros);
  for(int e=0;e<Nl;e++){
    double total=0;
    for(int k=0;k<Nm;k++){
      for(int k2=0;k2<Nm;k2++){
	total=total+Ak(e,k)*SIG(k,k2);
      }
    }
    sig(e)=total;
  }
  for(int e1=0; e1<Nl;e1++){
    for(int e2=0;e2<Nl;e2++){
      double t2=0;
      for(int k1=0;k1<Nm;k1++){
	for(int k2=0;k2<Nm;k2++){
	  t2=t2+Ak(e1,k1)*Ak(e2,k2)*SIG(k1,k2);
	}
      }
      sigger(e1,e2)=t2;
    }
  }
  //  Ak.print("Ak: ");
  //  sig.print("sig: ");
  //  sigger.diag().print("sigger: ");
  //  SIG.print("SIG: ");

  mat SIGy(Nl,Nl,fill::zeros);
  SIGy = Hw*Cm*SIG*(Hw*Cm).t();


  mat cy(Nl,Nl,fill::zeros);

  for(int i=0;i<Nl;i++){
    for(int j=0;j<Nl;j++){
      double pi_i = dot(Ag.row(i),alpha);
      double pi_j = dot(Ag.row(j),alpha);
      double total = pi_i*pi_j*TV - pi_i*sig(j)-pi_j*sig(i)+sigger(i,j);
      cy(i,j)=total;
      
    }
  }
  //  cy.print("cy: ");
  //  SIGy.print("SIGy: ");
  mat err=cy - SIGy;
  //  err.diag().print("err: ");
  cout<<"diag error: "<<accu(err.diag())<<endl;
  cout<<"total error: "<<accu(err)<<endl;


  ijn1 n1(gr, SIGy,Hw,L,p,pc,eps,epsN);


  mat var = n1.getVar();

  mat outvar(Nl,Nl,fill::zeros);
  for(int e=0;e<Nl;e++){
    for(int n=0;n<Nl;n++){
      double len = Lo(e,n);
      double pi_e = dot(Ag.row(e),alpha);
      double pi_n = dot(Ag.row(n),alpha);
      double phi_en = pi_e + Lo(e,n)*pi_n;

      //double vare = pi_e*pi_e*TV - 2*pi_e*sig(e) + sigger(e,e);
      //double varn = pi_n*pi_n*TV - 2*pi_n*sig(n) + sigger(n,n);
      //double covar = pi_e*pi_n*TV - pi_e*sig(n) -pi_n*sig(e) + sigger(e,n);
      //double total = vare + len*len*varn + 2*len*covar;
      double total = phi_en*phi_en*TV -2*phi_en*(sig(e)+len*sig(n))+ sigger(e,e) + 2*len*sigger(e,n)+len*len*sigger(n,n);

      if(!isfinite(total)) total=0;
      outvar(n,e)=total;
    }
  }
  mat errvar=var-outvar;
  //  var.row(0).print("var: ");
  //  outvar.row(0).print("var: ");
  //  errvar.row(0).print("errvar: ");
    cout<<"total error: "<<accu(errvar)<<endl;


    vec eN(Nl, fill::ones);
    eN = eN*epsN;

    isj sj(gr, &gc, SIG, Cm, L, p, pc, eps);
    

  return 0;
  //Set Probability info
  ranvar rv;


  igrid ig(gr);
  ig.addCost();
  isolve is;
  is.setSolver(IloCplex::Dual,IloCplex::Dual);
  rgrid * rbase = ig.solveModel(&is);
  //  rbase->displayOperatingPos(gr);

  try{
    in1 bn1(gr, SIGy, Hw, m1); 
    ijn1 n1(gr, SIGy,Hw,L,p,pc,eps,epsN);
    ijn1 jcc(gr, SIGy,Hw,L,p,pc,eps,1);

    n1.addCost();
    jcc.addCost();
    bn1.addCost();

    rgrid * rbn1 = bn1.solveModel(&is);
    rgrid * rjcc = jcc.solveModel(&is);
    rgrid * rn1_1 = n1.solveModel(&is);

    double o0=rbase->getObjective();
    vec f0=gc.convert(rbase->getF());
    vec g0=gc.convert(rbase->getG());
    vec z0=gc.risk(f0,SIGy.diag(),L,p,pc);
    double r0 = sum(z0);
    IloCplex::CplexStatus s0=rbase->getStatus();

    double o1=rn1_1->getObjective();
    vec f1=gc.convert(rn1_1->getF());
    vec g1=gc.convert(rn1_1->getG());
    vec z1=gc.risk(f1,SIGy.diag(),L,p,pc);
    double r1 = sum(z1);
    IloCplex::CplexStatus s1=rn1_1->getStatus();

    double o2=rjcc->getObjective();
    vec f2=gc.convert(rjcc->getF());
    vec g2=gc.convert(rjcc->getG());
    vec z2=gc.risk(f2,SIGy.diag(),L,p,pc);
    double r2 = sum(z2);
    IloCplex::CplexStatus s2=rjcc->getStatus();

    double o3=rbn1->getObjective();
    vec f3=gc.convert(rbn1->getF());
    vec g3=gc.convert(rbn1->getG());
    vec z3=gc.risk(f3,SIGy.diag(),L,p,pc);
    double r3 = sum(z3);
    IloCplex::CplexStatus s3=rbn1->getStatus();
    

    f0.t().print("f0: ");
    f1.t().print("f1: ");
    z0.t().print("z0: ");
    z1.t().print("z1: ");

    cout.precision(4);
    cout<<fixed<<r0<<" - "<<r1<<endl;


    cout<<o0<<" - "<<o1<<endl;
    cout<<"gen0: "<<rbase->getG()<<endl;
    cout<<"gen1: "<<rn1_1->getG()<<endl;
    
    double TG = accu(g0);

    running_stat<double> stats_r0;
    running_stat<double> stats_r1;
    running_stat<double> stats_r2;
    running_stat<double> stats_r3;

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
	vec f2n = n1.getN1(i,f2,g2);
	vec z2n=gc.risk(f2n,SIGy.diag(),L,p,pc);
	double r2n = sum(z2n);
	stats_r2(r2n);
	vec f3n = n1.getN1(i,f3,g3);
	vec z3n=gc.risk(f3n,SIGy.diag(),L,p,pc);
	double r3n = sum(z3n);
	stats_r3(r3n);
      }
    }
    cout.precision(5);
    cout<<fixed<<endl;
    cout<<"Total Gen: "<<TG<<endl;
    cout<<"Total Variance: "<<TV<<" ( "<<sqrt(TV)<<" )"<<endl;
    cout<<"Total Random Cost: "<<randcost<<endl;
    slack.t().print("slack: ");
    cout<<"\n\n";
    cout<<"OPF"<<"\t"<<s0<<endl;
    cout<<"C0: "<<o0<<endl;
    cout<<"r0 - "<<r0<<endl;
    cout << "count = " << stats_r0.count() << endl;
    cout << "mean = " << stats_r0.mean() << endl;
    cout << "stdv  = " << stats_r0.stddev()  << endl;
    cout << "min  = " << stats_r0.min()  << endl;
    cout << "max  = " << stats_r0.max()  << endl;
    cout<<endl;
    cout<<"JCC"<<"\t"<<s2<<endl;
    cout<<"C2: "<<o2<<endl;
    cout<<"r2 - "<<r2<<endl;
    cout << "count = " << stats_r2.count() << endl;
    cout << "mean = " << stats_r2.mean() << endl;
    cout << "stdv  = " << stats_r2.stddev()  << endl;
    cout << "min  = " << stats_r2.min()  << endl;
    cout << "max  = " << stats_r2.max()  << endl;
    cout<<endl;    
    cout<<"OPF N-1"<<"\t"<<s3<<endl;
    cout<<"C3: "<<o3<<endl;
    cout<<"r3 - "<<r3<<endl;
    cout << "count = " << stats_r3.count() << endl;
    cout << "mean = " << stats_r3.mean() << endl;
    cout << "stdv  = " << stats_r3.stddev()  << endl;
    cout << "min  = " << stats_r3.min()  << endl;
    cout << "max  = " << stats_r3.max()  << endl;
    cout<<endl;
    cout<<"JCC N-1"<<"\t"<<s1<<endl;
    cout<<"C1: "<<o1<<endl;
    cout<<"r1 - "<<r1<<endl;
    cout << "count = " << stats_r1.count() << endl;
    cout << "mean = " << stats_r1.mean() << endl;
    cout << "stdv  = " << stats_r1.stddev()  << endl;
    cout << "min  = " << stats_r1.min()  << endl;
    cout << "max  = " << stats_r1.max()  << endl;
    cout<<endl;

    for(int i=0;i<Ng;i++)   cerr<<alpha(i)<<"\t";
    cerr<<"\t"<<randcost<<"\t"<<TV<<"\t\t";
    cerr<<o0<<"\t"<<r0<<"\t"<<stats_r0.mean()<<"\t"<<stats_r0.max()<<"\t\t";
    cerr<<o2<<"\t"<<r2<<"\t"<<stats_r2.mean()<<"\t"<<stats_r2.max()<<"\t\t";
    cerr<<o3<<"\t"<<r3<<"\t"<<stats_r3.mean()<<"\t"<<stats_r3.max()<<"\t\t";
    cerr<<o1<<"\t"<<r1<<"\t"<<stats_r1.mean()<<"\t"<<stats_r1.max()<<endl;
    
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

   
