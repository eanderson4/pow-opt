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
  mat getH(){ return _H; };
  mat getHw(vec slack);
  vec convert(IloNumArray na);
  void testSlack();
  void test();

 private:
  grid * _gr;

  mat _H;
  

};


#endif
