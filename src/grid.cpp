#include "grid.h"

map<int, int> grid::map_busNum;


void grid::addGenCost(int num, int model, double startup, double shutdown, int ncost, double c2, double c1, double c0){
    _gens[num-1].addCost(model,startup,shutdown,ncost,c2,c1,c0);
}

int grid::getBusNum(int num){    
    if(!mapBus) buildMap();    
    return map_busNum[num];
}


void grid::buildMap(){
    for(int i=0;i<numBuses();i++){
      map_busNum[_buses[i].getNum()] = i;
    }
    mapBus=true;
}

void grid::printNums(ostream& stream){
    stream<<"Buses: "<<_buses.size()<<endl;
    stream<<"Branches: "<<_branches.size()<<endl;
    stream<<"Generators: "<<_gens.size()<<endl;
}
