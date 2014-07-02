#ifndef TEST_H
#define TEST_H

#include "gridcalc.h"
#include "rv.h"

using namespace std;

class test {

 public:
 test( grid * gr ): _gr(gr) {}

  void run();
  void RSLACK(double etol);
  void LODF(double etol);
  void PTDF(double etol);
  void RV(double etol);
  void RANDOM(double etol);
  void HESSIAN(double etol);
  void COVAR(double etol);

  class errortol: public exception 
  {
    virtual const char* what() const throw()
    {
      return "Error out of tolerance";
    }
  }errtol;

  
 private:
  grid * _gr;

};

#endif
