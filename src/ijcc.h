#ifndef IJCC_H
#define IJCC_H

#include "igrid.h"
#include "gridcalc.h"

using namespace std;


class ijcc : public igrid {

 public:
 ijcc(grid * gr,mat SIGy, double L, double p, double pc, double eps): igrid(gr), _SIGy(SIGy), _L(L), _p(p), _pc(pc), _eps(eps) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  void lineLimitStatus(bool status);

  //Risk functions
  
  class itlimit: public exception 
  {
    virtual const char* what() const throw()
    {
      return "Stopped due to iteration limit";
    }
  }itlimit;

 private:
  IloNumVarArray _z;
  IloNumVarArray _fplus;
  IloRange _riskConstraint;
  IloRangeArray _fup;
  IloRangeArray _fdown;

  mat _SIGy;
  double _L;
  double _p;
  double _pc;
  double _eps;
  
};
#endif
