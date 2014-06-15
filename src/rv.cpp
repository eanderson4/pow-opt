#include "rv.h"

double ranvar::g(double fnormalized, double L, double p, double pc){
  double f=fnormalized;
  double prob=0;
  double Uc=L + pc*(1-L)/p;
  double a=-p*L/(1-L);
  double b=p/(1-L);

  if(f<=L) prob=0;
  else if (f<Uc) prob=a+b*f;
  else prob=1;

  return prob;
}

double ranvar::simProb(double L, double p, double pc){
  double totalProb=0;
  double probInt=1/double(_N);

  double Uc=L + pc*(1-L)/p;
  double a=-p*L/(1-L);
  double b=p/(1-L);

  for(int i=0; i<getNum(); i++){
    double r=getValue(i);
    double g;

    if(r<=L)   g=0;
    else if (r<Uc){
      g= a + b*r;
    }
    else g=1;

    totalProb=totalProb+probInt*g;

  }
  return totalProb;
}

double ranvar::deriveMu(double L, double p, double pc, double mean, double stdv){
  double Uc=L + pc*(1-L)/p;
  double a=-p*L/(1-L);
  double b=p/(1-L);

  double mu = mean;
  double sigma = stdv;

  double alpha_low=(L - mu)/sigma;
  double alpha_high=(Uc - mu)/sigma;

  double beta=PHI(alpha_high) - PHI(alpha_low);
  double delta=phi(alpha_low)-phi(alpha_high);
  double tau=alpha_low*phi(alpha_low)-alpha_high*phi(alpha_high);
  


  double one=a*delta/sigma;
  double two = b*mu*delta/sigma + b*beta + b*tau;
  double three=phi(alpha_high)/sigma;

  double derive=one+two+three;

  return derive;

}

double ranvar::deriveSigma(double L, double p, double pc, double mean, double stdv){
  double Uc=L + pc*(1-L)/p;
  double a=-p*L/(1-L);
  double b=p/(1-L);

  double mu = mean;
  double sigma = stdv;

  double alpha_low=(L - mu)/sigma;
  double alpha_high=(Uc - mu)/sigma;

  double beta=PHI(alpha_high) - PHI(alpha_low);
  double delta=phi(alpha_low)-phi(alpha_high);
  double xi=mu + sigma*delta/beta;
  double tau=alpha_low*phi(alpha_low)-alpha_high*phi(alpha_high);
  double rho=alpha_low*alpha_low*phi(alpha_low) - alpha_high*alpha_high*phi(alpha_high);
  
  double one=a*tau/sigma;
  double two = b*xi*tau/sigma + b*delta+b*rho-b*delta*tau/beta;
  double three=alpha_high*phi(alpha_high)/sigma;

  double derive=one+two+three;

  return derive;

}

double ranvar::anaProb(double L, double p, double pc, double mean, double stdv){
  double Uc=L + pc*(1-L)/p;
  double a=-p*L/(1-L);
  double b=p/(1-L);

  double mu = mean;
  double sigma = stdv;
  double low = L;
  double high = Uc;
  
  double alpha_low=(low - mu)/sigma;
  double alpha_high=(high - mu)/sigma;
  double p_iL = PHI(alpha_low);
  double p_LU = PHI(alpha_high)-PHI(alpha_low);
  double ef_LU;
  if(p_LU==0) ef_LU = 0; 
  else ef_LU = mu+sigma*(( phi(alpha_low)- phi(alpha_high))/( PHI(alpha_high)- PHI(alpha_low)));

  double p_Ui = 1- PHI(alpha_high);
  double t_p = a*p_LU + b*ef_LU*p_LU + p_Ui;


  if(1==0){
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
  }

  return t_p;
}

double ranvar::sampleMean(){
  double mean=0;
  for(int i=0;i<getNum();i++){
    mean=mean+getValue(i);    
  }
  mean=mean/getNum();
  
  return mean;
}

double ranvar::sampleStDv(){
  double mean=sampleMean();

  double s2=0;
  for(int i=0;i<getNum();i++){
    s2=s2+pow((getValue(i)-mean),2);
  }
  s2=s2/(getNum()-1);
  s2=sqrt(s2);
  
  return s2;
}

void ranvar::createRV(int N,double mean, double stdv){
  _N=N;_mean=mean;_stdv=stdv;
  //Draw 0-1
  double zero_one = (rand() % RAND_MAX);
  zero_one = zero_one/RAND_MAX;

  //Gaussian
  std::default_random_engine generator;
  generator.seed(_seed);
  std::normal_distribution<double> distribution(mean,stdv);

  for(int i =0; i< N; i++){
    random_variable.push_back(distribution(generator));
  }
  
}

double ranvar::phi(double x){ 
  if(x >= std::numeric_limits<double>::max()) return 0;
  else return 1/sqrt(2*M_PI)*exp(-.5*pow(x,2)); 
}

double ranvar::PHI(double x){
  if(x >= std::numeric_limits<double>::max()) return 1;
  if(x <= std::numeric_limits<double>::lowest()) return 0;
  else   return 0.5 * erfc(-x * M_SQRT1_2);
}

void ranvar::testRV(){
  //Test purposes
  cout<<"Max error in gaussian CDF function"<<endl;
  testPhi();

  cout<<"Run simulation and compare to analytic"<<endl;
  cout<<"\nSimulation Prob\n"<<endl;
  int samples=53700;
  double mean=.9;
  double variance=.1;
  double stdv=sqrt(variance);
  
  cout<<"Samples: "<<samples<<endl;
  cout<<"Mean: "<<mean<<endl;
  cout<<"Var: "<<variance<<endl;
  

  createRV(samples,mean,stdv);

  double L=.9;
  double p=.15;
  double pc=.5;
  double totalProb=simProb(L,p,pc);
  cout<<"Sim Prob: "<<totalProb<<endl;

  //Analytic

  double t_p = anaProb(L,p,pc);
  cout<<"Analytic Prob: "<<t_p<<endl;
  cout<<"\nTotal Error: "<<t_p - totalProb<<endl;
  


}

void ranvar::testPhi()
{
  //CREDIT TO http://www.johndcook.com/cpp_phi.html
  // Select a few input values
  double x[] =
    {
      -3,
      -1,
      0.0,
      0.5,
      2.1
    };
  
  // Output computed by Mathematica
  // y = Phi[x]
  double y[] =
    {
      0.00134989803163,
      0.158655253931,
      0.5,
      0.691462461274,
      0.982135579437
    };
  
  int numTests = sizeof(x)/sizeof(double);
  
  double maxError = 0.0;
  for (int i = 0; i < numTests; ++i)
    {
      double error = fabs(y[i] - PHI(x[i]));
      if (error > maxError)
	maxError = error;
    }
  
  std::cout << "Maximum error: " << maxError << "\n";
} 

ostream& operator<<(ostream& os, const ranvar& rv)
{

  int nE=rv.random_variable.size();

  for(int i =0; i<nE; i++){
    os<<rv.random_variable[i]<<endl;
  }
  return os;
}
