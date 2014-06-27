#include <stdlib.h>

#include "sqlinter.h"
#include "test.h"

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

  test t(gr);
  t.run();
  
  cout<<"\n\n\n"<<endl;

  //N-1 calculations

  return 0;
}
   
