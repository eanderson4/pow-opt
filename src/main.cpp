#include <stdlib.h>

#include "sqlinter.h"
#include "grid.h"
#include "igrid.h"
#include "rgrid.h"

using namespace std;

int main(int argc, char* argv[]){

  if(argc>1){
    try {
      sqlInter db;
      grid * gr = new grid;
      
      string db_name;

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

      delete gr;     
    }
    catch (exception& e) {
      cerr << e.what() << "\n";
    }
  }
  else
    cout<<"Command: pow case/30.db\n"
	<<"\trun power flow for case30"<<endl;
  
  


  return 0;
}
   
