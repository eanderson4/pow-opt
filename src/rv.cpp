#include "rv.h"

void ranvar::createRV(int N,double mean, double stdv){

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
