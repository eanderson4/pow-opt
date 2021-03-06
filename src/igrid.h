#ifndef IGRID_H
#define IGRID_H

#include "ibase.h"
#include "rgrid.h"
#include "del_g.h"

class igrid : public ibase {

 public:

 igrid(grid * gr ) : ibase ( ), have_cost(false),have_loadshed(false){
    _gr = gr;
    if(buildModel( gr  )) cout<<"IloModel Built\n";
    else cout<<"Could not build IloModel\n";
  }
  virtual ~igrid() {  }

  virtual rgrid * solveModel( isolve * is=NULL);
  int addCost();
  void allowLoadShed();
  void modGrid( del_g mod );

  double getCost(IloNumArray g){ return _ic.getCost(_gr,g);  }
  grid * getGrid(){ return _gr; }
  icost getIcost(){ return _ic; }
  ished getIshed(){ return _ils; }

  void relayInfo(rgrid * rg);  

 private:
  
  grid * _gr;

  icost _ic;
  ished _ils;

  bool have_cost;
  bool have_loadshed;

};


#endif
