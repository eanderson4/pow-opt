#ifndef RV_H
#define RV_H


#include <iostream>
#include <stdlib.h>
#include <random>
#include <vector>

using namespace std;



class ranvar
{

 public:
 ranvar(int seed=372):_seed(seed) { }
  ~ranvar() {}

  void createRV(int N,double mean, double stdv);
  int getNum(){ return random_variable.size(); }
  double getValue(int N){ return random_variable[N]; }

  friend ostream& operator<<(ostream& os, const ranvar& rv);

 private:

  vector<double> random_variable;
  int _seed;
  
};

#endif
