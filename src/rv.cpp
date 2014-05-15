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

ostream& operator<<(ostream& os, const ranvar& rv)
{

  int nE=rv.random_variable.size();

  for(int i =0; i<nE; i++){
    os<<rv.random_variable[i]<<endl;
  }
  return os;
}
