#ifndef GRIDCALC_H
#define GRIDCALC_H

#include <iostream>
#include <stdlib.h>
#include <ilcplex/ilocplex.h>
#include "grid.h"
#include "igrid.h"
#include "rgrid.h"
#include "rv.h"
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

  vec getDelG(del_g dg);
  vec getDelF(vec delg, vec slack);
  mat getC(){ return _C; };  //Admittance Matrix
  mat getH(){ return _H; };  //Linear Shift Factor
  mat getHw(vec slack);      //with arbitray slack distribution
  mat getL(mat Hw);          //Line outage distribution factor
  vec convert(IloNumArray na);

 private:
  grid * _gr;

  mat _H;
  mat _C;

};


#endif
