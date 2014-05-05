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
     

      /*
      cout<<ig->getNodalBalance()[139]<<endl;
      cout<<ig->getNodalBalance()[140]<<endl;
      cout<<ig->getNodalBalance()[141]<<endl;
      cout<<rg->getF()[241]<<endl;
      cout<<ig->getNodalBalance()[275]<<endl;
      cout<<ig->getBranchFlow()[19]<<endl;
      cout<<ig->getNodalBalance()[277]<<endl;
      cout<<ig->getNodalBalance()[278]<<endl;
      */

      int nB=gr->numBuses();
      int nR=gr->numBranches();
      int nG=gr->numGens();



      //Outputing resolts to JSON format
      cerr<<"{\n"
	  <<"\t\"dataset\": { \n"
	  <<"\t\t\"nodes\": [ \n";
      for(int i =0; i<nB; i++){
	bus b = gr->getBus(i);
	double ni=0; //net inject
	cerr<<"\t\t\t { \"name\": \""<<b.getNum()<<"\",\n"
	    <<"\t\t\t   \"index\": "<<i<<",\n"
	    <<"\t\t\t   \"demand\": "<<b.getP()<<",\n";
	ni=ni+b.getP();

	for(int j =0; j<nG; j++){
	  gen g =gr->getGen(j); 
	  int bus = g.getBus();
	  if(bus-1==i){
	    double p = rg->getG()[j];
	    if( p > 0.001 ){
	      cerr <<"\t\t\t   \"p\": "<<p<<",\n";
	      ni=ni-p;

	    }
	  }
	}
	cerr <<"\t\t\t   \"ni\": "<<ni<<" ";
	if(i!=nB-1) cerr<<"},"<<endl;
	else  cerr<<"} \n"
		 <<"\t\t],"<<endl;

      }
      cerr<<"\t\t \"edges\": [ \n";
      for(int i =0; i<nR; i++){
	branch b = gr->getBranch(i);
	int from= gr->getBusNum(b.getFrom());
	int to= gr->getBusNum(b.getTo());
	double X=gr->getBranch(i).getX();
	double flow=double(rg->getF()[i]);
	cerr<<"\t\t\t { \"source\": "<<from<<", \"target\": "<<to<<","<<endl;
	cerr<<"\t\t\t   \"index\": "<<i<<",\n";
	cerr<<"\t\t\t   \"X\": "<<X<<",\n";
	cerr<<"\t\t\t   \"flow\": "<<flow;
	if (i!=nR-1) cerr<<" },"<<endl;
	else cerr<<" }"<<endl;
      }
           cerr<<"\t\t ],\n"
	  <<"\t\t \"gens\": [ \n";
      for(int i =0; i<nG; i++){
	gen g =gr->getGen(i); 
	int bus = g.getBus();
	double p = g.getP();
	if( p > 0.001 ){
	  cerr<<"\t\t\t { \"bus\": "<<bus<<",\n"
	      <<"\t\t\t   \"p\": "<<p;
	  if (i!=nG-1) cerr<<" },"<<endl;
	  else cerr<<" }"<<endl;
	}
	}
      cerr<<"\t\t ]\n"
	  <<"\t} \n"
	  <<"}"<<endl;
      



    }
    catch (exception& e) {
      cerr << e.what() << "\n";
    }
  }
  
  delete gr;

  return 0;
}
   
