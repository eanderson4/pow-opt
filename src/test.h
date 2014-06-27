#ifndef TEST_H
#define TEST_H

#include <exception>
#include "gridcalc.h"

using namespace std;

class errortol: public exception
{
  virtual const char* what() const throw()
  {
    return "Error out of tolerance";
  }
} errtol;

class test {

 public:
 test( grid * gr ): _gr(gr) {}

  void run();
  void LODF(double etol);
  
 private:
  grid * _gr;

};

#endif
