#ifndef ISJ_H
#define ISJ_H

#include "igrid.h"
#include "gridcalc.h"

using namespace std;


class isj : public igrid {

 public:
 isj(grid * gr, gridcalc * gc, mat SIGm, mat Cm,double L, double p, double pc, double eps): igrid(gr), _gc(gc), _SIGm(SIGm), _Cm(Cm) , _L(L), _p(p), _pc(pc), _eps(eps) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  void lineLimitStatus(bool status);
  bool postCC(vec f, vec z,IloCplex * cplex, int iteration=0);
  mat getSIGm(){ return _SIGm; }
  mat getA(){ return _gc->getH(); }
  mat getCb(){ return _gc->getC(); }
  mat getCm(){ return _Cm; }
  mat getCg(){ return _gc->getCm(); }
  double getEps(){ return _eps; }
  double getL(){ return _L; }
  double getP(){ return _p; }
  double getPc(){ return _pc; }
  IloNumVarArray getZ(){ return _z; }
  IloNumVarArray getYplus(){ return _yplus; }

  //Risk functions
  
  class iterlimit: public exception 
  {
    virtual const char* what() const throw()
    {
      return "Cant solve this version";
    }
  }itlimit;

 private:
  gridcalc * _gc;

  IloNumVar _risk;
  IloRange _riskConstraint;
  IloNumVarArray _z;
  IloNumVarArray _yplus;
  IloNumVarArray _sd;
  IloNumVarArray _beta;
  IloNumVarArray _pi;
  IloRangeArray _yup;
  IloRangeArray _ydown;
  IloRangeArray _pibeta;
  IloRangeArray _sdfe;


  mat _SIGm;
  mat _Cm;
  double _L;
  double _p;
  double _pc;
  double _eps;

  vec sig_e;
  vec sigger_ee;
  
  vec _in;
  
};
#endif
