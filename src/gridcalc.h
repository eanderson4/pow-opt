#ifndef GRIDCALC_H
#define GRIDCALC_H

#include <iostream>
#include <stdlib.h>
#include "grid.h"
#include "del_g.h"

#define ARMA_DONT_USE_WRAPPER
#include <armadillo>


using namespace arma;
using namespace std;

class gridcalc
{

  public:
  gridcalc(grid * gr);
  ~gridcalc() {}

  vec getDelF(vec delg, vec slack);
  void testSlack();

 private:
  grid * _gr;

  mat _H;
  

};


#endif
