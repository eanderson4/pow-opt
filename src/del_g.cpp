#include "del_g.h"

void del_g::setOutage(int num){ 
  if(!have_topo){
    cerr<<"Set base topology before allowing changes to it"<<endl;
    return;
  }
  del_topo[num]=false; 
}
void del_g::baseTopo(grid * gr){
  int nB = gr->numBranches();
  
  for(int i=0; i<nB;i++){
    del_topo.push_back(true);
  }

  have_topo=true;
}
