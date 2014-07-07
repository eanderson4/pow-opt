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
  bool postN1(int n,vec yn, IloCplex * cplex);
  vec getN1(int n, vec y0, vec g0,bool add=false);

 private:
  vector<IloNumVarArray> _z;
  vector<IloNumVarArray> _yplus;
  vector<IloRangeArray> _yup;
  vector<IloRangeArray> _ydown;

  IloRangeArray _riskConstraint;

  mat _SIGy;
  mat _Hw;  

};
#endif
