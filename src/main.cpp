#include <stdlib.h>

#include "sqlinter.h"
#include "grid.h"
#include "igrid.h"
#include "rgrid.h"

using namespace std;

int main(int argc, char* argv[]){

  sqlInter db;
  grid * gr = new grid;

  string db_name;

  if(argc>1){
    try {
      db_name = argv[1];
      db.openDb(db_name);
      db.load(*gr);
      gr->buildMap();	
      gr->printNums(cout);      	
      
      igrid * ig = new igrid( gr);
      ig->addCost();

      rgrid * rg = ig->solveModel();
      rg->displayOperatingPos( gr );

      
      del_g del;
      del.baseTopo(gr);
      del.setOutage(15);
      
      ig->modGrid( del );     
     
      rg = ig->solveModel( );
      rg->displayOperatingPos( gr );
      
      
      int nB=gr->numBuses();
      int nR=gr->numBranches();
      int nG=gr->numGens();
      for(int i =0; i<nB; i++){
	cout<<gr->getBus(i)<<endl;
      }
      for(int i =0; i<nR; i++){
	cout<<gr->getBranch(i)<<endl;
      }
      for(int i =0; i<nG; i++){
	cout<<gr->getGen(i)<<endl;
      }

    }
    catch (exception& e) {
      cerr << e.what() << "\n";
    }
  }
  
  delete gr;

  return 0;
}
   
