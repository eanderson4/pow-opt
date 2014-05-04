#ifndef DELG_H
#define DELG_H

#include <vector>
#include "grid.h"
using namespace std;


class del_g
{
 public:
  del_g(): have_demand(false), have_topo(false), have_cap(false){ }
  ~del_g(){}
  
  void setOutage(int num);
  void baseTopo(grid * gr);

  bool haveTopo(){ return have_topo; }
  double getStatus(int num){ return del_topo[num]; }



 private:

  bool have_demand;
  bool have_topo;
  bool have_cap;

  vector< double > del_demand;
  vector< bool > del_topo;
  vector< double > del_cap;

};

#endif
