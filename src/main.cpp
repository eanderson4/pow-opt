#include <stdlib.h>

#include "sqlinter.h"
#include "igrid.h"
#include "grid.h"
#include "gridcalc.h"
#include "rv.h"

using namespace std;



int main(int argc, char* argv[]){

  if(argc<=1){
    cout<<"cmd: pow case/30.db\n"
	<<"\trun main for case30"<<endl;
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

  cout<<*gr<<endl;

  gridcalc gc(gr);

  int Nb=gr->numBuses();
  vec delg(Nb,fill::zeros);
  vec slackdist(Nb,fill::zeros);
  delg(4)=5; delg(6)=7;
  slackdist(1)=1;

  vec del_f = gc.getDelF(delg,slackdist);
  
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
  

  ig.modGrid(dg);
    
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
  
  cout<<"delF (sim): "<<endl;
  for(int i=0;i<gr->numBranches();i++)
    cout<<" "<<(rg->getF()[i] - rg2->getF()[i])<<"   ";

  cout<<"\ndelF (sim) - delF (shift factor)"<<endl;
  double totalerror=0;
  for(int i=0;i<gr->numBranches();i++){
    cout<<" "<<(rg->getF()[i] - rg2->getF()[i] - del_f(i))<<"   ";
    totalerror=totalerror+(rg->getF()[i] - rg2->getF()[i] - del_f(i));
  }
  cout<<"\n";
  cout<<"Total Error: "<<totalerror<<endl;

  ranvar rv;
  rv.testRV();  

  return 0;
}
   
