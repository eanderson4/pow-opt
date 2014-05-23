#include <stdlib.h>

#include "sqlinter.h"
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




  return 0;
}
   
