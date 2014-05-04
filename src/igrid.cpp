#include "igrid.h"

rgrid *  igrid::solveModel( isolve * is){

  rgrid * rg = new rgrid();
  
  time_t tstart;
  tstart = clock();
    
  IloCplex cplex(*getModel());
  
  if(is!=NULL) is->setCplexParams(&cplex);
  
  if (cplex.solve()){
      float tot= float(clock() - tstart) / CLOCKS_PER_SEC;
      cerr<<"MODEL solved in "<<tot<<endl;
      rg->getSolveInfo(&cplex,tot);
      getBaseResults(&cplex, rg);
      if(have_loadshed) getIshed().getLoadShed(&cplex, rg);
      
      cerr<<"OBJECTIVE: "<<cplex.getObjValue()<<endl;
      double genCost = getIcost().getCost(_gr,rg->getG());
      rg->setGenCost(genCost);
      
  }
  else{
    cerr<<"Not solved"<<endl;
  }
  cplex.end();
  
  return rg;
}

void igrid::allowLoadShed(){
  have_loadshed = _ils.buildLoadShed(_gr,getModel(),getNodalBalance());
  if(have_loadshed) cout<<"Built load shed variables"<<endl;
}

int igrid::addCost(){ 
  if(have_loadshed){
    have_cost=_ic.buildCostWithLoadShed(_gr,getModel(),getG(),_ils.getLS());
  }
  else{
    have_cost=_ic.buildCost(_gr,getModel(),getG());
  }
  if(have_cost) cout<<"Built cost object"<<endl;
  return have_cost; 
}


void igrid::modGrid( del_g mod ){
  
  //  make modification
  int nB = _gr->numBranches();

  if(mod.haveTopo()){
    
    for(int i=0; i<nB; i++){
      if(mod.getStatus(i)==false){
	// line i is outaged
	cout<<"Line "<<i<<" - ";
	getF()[i].setBounds(0,0);
	getBranchFlow()[i].setBounds(-IloInfinity,IloInfinity );
	getPhaseAngle()[i].setBounds(-IloInfinity,IloInfinity );
	cout<<"f=0, phase angles relaxed"<<endl;
      }
    }
    
  }


}

