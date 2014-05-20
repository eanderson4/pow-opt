#include <stdlib.h>

#include "sqlinter.h"
#include "grid.h"
#include "igrid.h"
#include "rgrid.h"
#include "rv.h"

using namespace std;



int main(int argc, char* argv[]){
  if(argc<=3){
    return 1;
  }

  int samples=atoi( argv[1] );
  double mean=atof( argv[2] );
  double variance=atof( argv[3] );
  double stdv=sqrt(variance);
  
  cout<<"Samples: "<<samples<<endl;
  cout<<"Mean: "<<mean<<endl;
  cout<<"Var: "<<variance<<endl;
  
  cerr<<argv[4]<<endl;

  ranvar rv(atoi(argv[4]));
  rv.createRV(samples,mean,stdv);

  //  cerr<<"flow,fail"<<endl;
  //  cerr<<rv<<endl;

  double totalProb=0;
  double probInt=1/double(samples);

  double L=.9;
  double p=.15;
  double pc=.5;
  double Uc=L + pc*(1-L)/p;
  double a=-p*L/(1-L);
  double b=p/(1-L);

  for(int i=0; i<rv.getNum(); i++){
    double r=rv.getValue(i);
    double g;

    if(r<=L)   g=0;
    else if (r<Uc){
      g= a + b*r;
    }
    else g=1;
    //    cerr<<r<<","<<g<<endl;

    totalProb=totalProb+probInt*g;

  }
  cout<<"Total Prob: "<<totalProb<<endl;
  //  cerr<<"Total Prob: "<<totalProb<<endl;
  
  


  //Analytic
  double mu = mean;
  double sigma = stdv;
  double low = L;
  double high = Uc;
  
  double alpha_low=(low - mu)/sigma;
  double alpha_high=(high - mu)/sigma;
  double p_iL = rv.PHI(alpha_low);
  double p_LU = rv.PHI(alpha_high)-rv.PHI(alpha_low);
  double ef_LU;
  if(p_LU==0) ef_LU = 0; 
  else ef_LU = mu+sigma*((rv.phi(alpha_low)-rv.phi(alpha_high))/(rv.PHI(alpha_high)-rv.PHI(alpha_low)));
  cout<<rv.phi(alpha_low)<<" "<<rv.phi(alpha_high)<<endl;
  cout<<rv.PHI(alpha_low)<<" "<<rv.PHI(alpha_high)<<endl;
  cout<<rv.PHI(alpha_high)-rv.PHI(alpha_low)<<endl;
  double p_Ui = 1-rv.PHI(alpha_high);
  
  double t_p = a*p_LU + b*ef_LU*p_LU + p_Ui;

  cout<<"\nAnalytic Prob"<<endl;
  cout<<"\nmu: "<<mu<<endl;
  cout<<"sigma: "<<sigma<<endl;
  cout<<"L: "<<L<<endl;
  cout<<"Uc: "<<Uc<<endl;
  cout<<"alpha_l: "<<alpha_low<<endl;
  cout<<"alpha_h: "<<alpha_high<<endl;
  cout<<"\nProb f in [-inf,L]: "<<p_iL<<endl;
  cout<<"Prob f in [L,Uc]: "<<p_LU<<endl;
  cout<<"Prob f in [Uc,inf]: "<<p_Ui<<endl;
  cout<<"\nE[f|f in [L,Uc]]: "<<ef_LU<<endl;
  cout<<"\na*Prob f in [L,Uc]: "<<a*p_LU<<endl;
  cout<<"b*E[f|f in [L,Uc]*Prob f in [L,Uc]: "<<b*ef_LU*p_LU<<endl;
  cout<<"1*Prob f in [Uc,inf]: "<<p_Ui<<endl;

  cout<<"\nTotal: "<<t_p<<endl;

  cout<<"\nTotal Error: "<<t_p - totalProb<<endl;
  
 

  //   rv.testPhi();





  /*
  sqlInter db;
  grid * gr = new grid;

  string db_name;

  db_name = argv[1];
  db.openDb(db_name);
  db.load(*gr);
  gr->buildMap();	
  gr->printNums(cout);      	
  
  igrid * ig = new igrid( gr);
  ig->addCost();
  
  rgrid * rg = ig->solveModel();
  rg->displayOperatingPos( gr );
  
  
  del_g del;
  del.baseTopo(gr);
      
  ig->modGrid( del );     
  
  rg = ig->solveModel( );
  rg->displayOperatingPos( gr );
  

  cout<<ig->getNodalBalance()[139]<<endl;
  cout<<ig->getNodalBalance()[140]<<endl;
  cout<<ig->getNodalBalance()[141]<<endl;
  cout<<rg->getF()[241]<<endl;
  cout<<ig->getNodalBalance()[275]<<endl;
  cout<<ig->getBranchFlow()[19]<<endl;
  cout<<ig->getNodalBalance()[277]<<endl;
  cout<<ig->getNodalBalance()[278]<<endl;
  
  int nB=gr->numBuses();
  int nR=gr->numBranches();
  int nG=gr->numGens();

  
  
  //Outputing resolts to JSON format
      cerr<<"{\n"
	  <<"\t\"dataset\": { \n"
	  <<"\t\t\"nodes\": [ \n";
      for(int i =0; i<nB; i++){
	bus b = gr->getBus(i);
	double ni=0; //net inject
	cerr<<"\t\t\t { \"name\": \""<<b.getNum()<<"\",\n"
	    <<"\t\t\t   \"index\": "<<i<<",\n"
	    <<"\t\t\t   \"demand\": "<<b.getP()<<",\n";
	ni=ni+b.getP();

	for(int j =0; j<nG; j++){
	  gen g =gr->getGen(j); 
	  int bus = g.getBus();
	  if(bus-1==i){
	    double p = rg->getG()[j];
	    if( p > 0.001 ){
	      cerr <<"\t\t\t   \"p\": "<<p<<",\n";
	      ni=ni-p;

	    }
	  }
	}
	cerr <<"\t\t\t   \"ni\": "<<ni<<" ";
	if(i!=nB-1) cerr<<"},"<<endl;
	else  cerr<<"} \n"
		 <<"\t\t],"<<endl;

      }
      cerr<<"\t\t \"edges\": [ \n";
      for(int i =0; i<nR; i++){
	branch b = gr->getBranch(i);
	int from= gr->getBusNum(b.getFrom());
	int to= gr->getBusNum(b.getTo());
	double X=gr->getBranch(i).getX();
	double flow=double(rg->getF()[i]);
	cerr<<"\t\t\t { \"source\": "<<from<<", \"target\": "<<to<<","<<endl;
	cerr<<"\t\t\t   \"index\": "<<i<<",\n";
	cerr<<"\t\t\t   \"X\": "<<X<<",\n";
	cerr<<"\t\t\t   \"flow\": "<<flow;
	if (i!=nR-1) cerr<<" },"<<endl;
	else cerr<<" }"<<endl;
      }
           cerr<<"\t\t ],\n"
	  <<"\t\t \"gens\": [ \n";
      for(int i =0; i<nG; i++){
	gen g =gr->getGen(i); 
	int bus = g.getBus();
	double p = g.getP();
	if( p > 0.001 ){
	  cerr<<"\t\t\t { \"bus\": "<<bus<<",\n"
	      <<"\t\t\t   \"p\": "<<p;
	  if (i!=nG-1) cerr<<" },"<<endl;
	  else cerr<<" }"<<endl;
	}
	}
      cerr<<"\t\t ]\n"
	  <<"\t} \n"
	  <<"}"<<endl;
      



    }
    catch (exception& e) {
      cerr << e.what() << "\n";
    }
  }
  
  delete gr;
  */
  return 0;
}
   
