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
      
      ig->modGrid( del );     
     
      rg = ig->solveModel( );
      rg->displayOperatingPos( gr );
     

    catch (exception& e) {
      cerr << e.what() << "\n";
    }
  }
  
  delete gr;

  return 0;
}
   
