#include "gen.h"

ostream& operator<<(ostream& os, const gen& g)
{
  os<<g._bus<<": "<<g._status<<"\t"<<g._pg<<" in ["<<g._pmin<<","<<g._pmax<<"]\t";
  os<<g._c2<<"x^2+"<<g._c1<<"x+"<<g._c0<<endl;
  return os;
}
