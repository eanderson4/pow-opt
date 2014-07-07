#ifndef IN1_H
#define IN1_H

#include "gridcalc.h"
#include "igrid.h"

using namespace std;


class in1 : public igrid {

 public:
 in1(grid * gr,mat SIGy, mat Hw): igrid(gr), _SIGy(SIGy), _Hw(Hw) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  
 private:
  vector<IloNumVarArray> _yplus;
  vector<IloRangeArray> _fup;
  vector<IloRangeArray> _fdown;

  mat _SIGy;
  mat _Hw;  

};
#endif
