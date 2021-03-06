#ifndef GRID_H
#define GRID_H

#include <iostream>
#include <vector>
#include <map>
#include "bus.h"
#include "branch.h"
#include "gen.h"

using namespace std;



class grid
{
  friend class bus;
  friend class branch;
  friend class gen;

 public:
 grid():mapBus( false ), loadShedPenalty(1000) {}
  ~grid() {}

  void setCapacity(int n, double U){    _branches[n].setRateA(U);  }
  void addBus(bus b){    _buses.push_back(b);  }
  void addBranch(branch b){    _branches.push_back(b);  }
  void addGen(gen g){    _gens.push_back(g);  }
  void addGenCost(int num, int model, double startup, double shutdown, int ncost, double c2, double c1, double c0);
  branch getBranch(int num){    return _branches[num];  }
  bus getBus(int num){    return _buses[num];  }
  gen getGen(int num){    return _gens[num];  }
  int numBuses(){    return _buses.size();  }
  int numBranches(){    return _branches.size();  }
  int numGens(){ return _gens.size();  }
  int getBusNum(int num);
  void buildMap();
  void printNums(ostream& stream);
  double getLoadShedPenalty(){ return loadShedPenalty; }

  friend ostream& operator<<(ostream& os, const grid& gr);

 private:
  vector<bus> _buses;
  vector<branch> _branches;
  vector<gen> _gens;

  static map<int, int> map_busNum;

  bool mapBus;

  double loadShedPenalty;

};

#endif
