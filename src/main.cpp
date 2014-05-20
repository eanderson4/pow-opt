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
  
 
  // Test purposes
 //   rv.testPhi();


  return 0;
}
   
