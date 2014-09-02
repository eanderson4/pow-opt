#ifndef ISJN_H
#define ISJN_H

#include "isj.h"
#include "igrid.h"

using namespace std;


class isjn : public isj {

 public:
 isjn(grid * gr,gridcalc * gc, mat SIGm,mat Cm, double L, double p, double pc, double eps, vec epsN): isj( gr, gc, SIGm, Cm, L, p, pc, eps ), _epsN(epsN) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  bool postN1(int n,vec f, vec g,vec z, IloCplex * cplex, int iteration=0);
  vec getN1(int n, vec y0, vec g);
  vec getCheck(){ return _check; }

  class nowinner: public exception 
  {
    virtual const char* what() const throw()
    {
      return "No islanded bus found!";
    }
  }nowin;


 private:
  IloNumVarArray _risk;
  IloRangeArray _riskConstraint;
  vector<IloNumVarArray> _z;
  vector<IloNumVarArray> _yplus;
  vector<IloNumVarArray> _sd;
  vector<IloNumVarArray> _psi;
  vector<IloRangeArray> _yup;
  vector<IloRangeArray> _ydown;
  vector<IloRangeArray> _psipi;
  vector<IloRangeArray> _sdfe;
  
  vec _epsN;

  mat _L;
  mat _Hb;

  vec _check;
  mat _in;


};
#endif