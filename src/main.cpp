#include <stdlib.h>

#include "sqlinter.h"
#include "igrid.h"
#include "grid.h"
#include "gridcalc.h"

using namespace std;



int main(int argc, char* argv[]){
  if(argc<1){
    return 1;
  }

  sqlInter db;
  grid * gr = new grid;
  string db_name;

  db_name = argv[1];
  db.openDb(db_name);
  db.load(*gr);
  gr->buildMap();
  gr->printNums(cout);

  gridcalc gc(gr);
  vec del_f = gc.getDelF();

  del_g dg(gr);
  dg.addDemand(4,5);
  dg.addDemand(6,7);
  
  igrid ig(gr);
  ig.addCost();
  rgrid * rg;
  rgrid * rg2;

  rg = ig.solveModel();
  IloNumArray g_nom=rg->getG();
  cout<<"Before"<<endl;
  cout<<"F: "<<rg->getF()<<endl;
  cout<<"G: "<<g_nom<<endl;
  

  //gr->modGrid(dg);
  ig.modGrid(dg);
    
  int Nb=gr->numBuses();
  IloNumArray slack(IloEnv(),Nb);
  for(int i=0;i<Nb;i++) slack[i]=0;
  slack[1]=1;
  
  ig.addSlack(g_nom,slack);

  rg2 = ig.solveModel();
  cout<<"After"<<endl;
  cout<<"F: "<<rg2->getF()<<endl;
  cout<<"G: "<<rg2->getG()<<endl;

  cout<<"Total Demand: "<<gr->getTotalDemand()<<endl;
  cout<<"Total Gen: "<<IloSum(rg2->getG())<<endl;
  
  cout<<"delF: "<<endl;
  for(int i=0;i<gr->numBranches();i++)
    cout<<" "<<(rg->getF()[i] - rg2->getF()[i])<<"   ";

  cout<<"\n\n\n";
  for(int i=0;i<gr->numBranches();i++)
    cout<<" "<<(rg->getF()[i] - rg2->getF()[i] - del_f(i))<<"   ";
  cout<<"\n";

  return 0;
}
   
