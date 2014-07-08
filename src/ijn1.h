#ifndef IJN1_H
#define IJN1_H

#include "ijcc.h"
#include "igrid.h"

using namespace std;


class ijn1 : public ijcc {

 public:
 ijn1(grid * gr,mat SIGy, mat Hw, double L, double p, double pc, double eps, double epsN): ijcc(gr,SIGy,L,p,pc,eps), _Hw(Hw),_epsN(epsN) { setup(); }
  void setup();
  rgrid * solveModel( isolve * is=NULL);
  bool postN1(int n,vec f, vec g,vec z, IloCplex * cplex);
  vec getN1(int n, vec y0, vec g);
  
  class nowinner: public exception 
  {
    virtual const char* what() const throw()
    {
      return "No islanded bus found!";
    }
  }nowin;


 private:
  vector<IloNumVarArray> _z;
  vector<IloNumVarArray> _yplus;
  vector<IloRangeArray> _yup;
  vector<IloRangeArray> _ydown;
  IloRangeArray _riskConstraint;

  mat _Hw;
  mat _L;
  mat _C;
  mat _Cg;
  mat _Hb;
  mat _var;
  vec _var0;
  double _epsN;

};
#endif
