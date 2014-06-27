#include <stdlib.h>

#include "sqlinter.h"
#include "igrid.h"
#include "grid.h"
#include "gridcalc.h"
#include "rv.h"
#include "test.h"

using namespace std;

int main(int argc, char* argv[]){
  
  if(argc<=1){
    cout<<"cmd: pow case/30.db\n"
	<<"\trun main for case30"<<endl;
    return 1;
  }

  //Load Grid
  sqlInter db;
  grid * gr = new grid;
  string db_name;
  db_name = argv[1];
  db.openDb(db_name);
  db.load(*gr);
  
  test t(gr);
  t.run();
    
  return 0;
}
