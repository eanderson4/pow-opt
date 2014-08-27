#ifndef ISJN_H
#define ISJN_H

#include "isj.h"
#include "igrid.h"

using namespace std;


class isjn : public isj {

 public:
 isjn(grid * gr,mat SIGm, mat A, double L, double p, double pc, double eps, double epsN): isj( gr, SIGm ), _A(A),_L(L),_p_epsN(epsN) { setup(); }
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
  
  vec _epsN;

  mat _L;
  mat _Hb;

  vec _check;
  mat _in;
  mat _addCut;

};
#endif
