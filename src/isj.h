#ifndef ISJ_H
#define ISJ_H

#include "igrid.h"
#include "gridcalc.h"

using namespace std;


class isj : public igrid {

 public:
 isj(grid * gr, gridcalc * gc, mat SIGm, vec indexM,double L, double p, double pc, double eps): igrid(gr), _gc(gc), _SIG(SIGm), _indexM(indexM) , _L(L), _p(p), _pc(pc), _eps(eps) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  void lineLimitStatus(bool status);
  bool postCC(vec f, vec z, vec beta, vec SIGy,IloCplex * cplex, int iteration=0);
  mat getSIGm(){ return _SIG; }
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
  
  vec getBeta(){ return _betaSolve; }
  vec getSD(){ return _sdSolve; }

  //Risk functions
  
  class iterlimit: public exception 
  {
    virtual const char* what() const throw()
    {
      return "Hit ITERATION LIMIT";
    }
  }itlimit;

 private:
  gridcalc * _gc;

  IloNumVar _risk;
  IloRange _riskConstraint;
  IloRange _betaSum;
  IloRange _budgetConstrain;
  IloNumVarArray _nu;
  IloNumVarArray _z;
  IloNumVarArray _yplus;
  IloNumVarArray _sd;
  IloNumVarArray _beta;
  IloNumVarArray _pi;
  IloRangeArray _yup;
  IloRangeArray _ydown;
  IloRangeArray _pibeta;
  IloRangeArray _sdfe;


  mat _SIG;
  vec _indexM;
  mat _Cm;
  double _L;
  double _p;
  double _pc;
  double _eps;

  mat _A;

  double _sig_delta;
  vec _sig;
  mat _sigger;
  
  vec _addCut;

  vec _betaSolve;
  vec _sdSolve;
  
};
#endif
