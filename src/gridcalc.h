#ifndef GRIDCALC_H
#define GRIDCALC_H

#include <iostream>
#include <stdlib.h>
#include "grid.h"

#define ARMA_DONT_USE_WRAPPER
#include <armadillo>


using namespace arma;
using namespace std;

class gridcalc
{

  public:
  gridcalc(grid * gr);
  ~gridcalc() {}

 private:
  grid * _gr;

};


#endif
