#ifndef IN1_H
#define IN1_H

#include "ijcc.h"
#include "igrid.h"

using namespace std;


class in1 : public ijcc {

 public:
 in1(grid * gr,mat SIGy, double L, double p, double pc, double eps, double epsN): ijcc(gr,SIGy,L,p,pc,eps), _epsN(epsN) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  
 private:
  vector<IloNumVarArray> _z;
  vector<IloNumVarArray> _yplus;
  IloRangeArray _riskConstraint;
  vector<IloRangeArray> _fup;
  vector<IloRangeArray> _fdown;

  double _epsN;

};
#endif
