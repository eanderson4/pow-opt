#include "gridcalc.h"

gridcalc::gridcalc(grid * gr) { 
  cout<<"gridcalc constructor"<<endl;
  cout<<"create useful objects"<<endl;

  _gr=gr;


  int Nb=gr->numBuses();
  int Nl=gr->numBranches();

  mat C(Nl,Nb);
  mat Bff(Nl,Nl);
  int slack;

  cout<<"Line: To bus - From bus  (indicies)"<<endl;
  for(int i =0; i<Nl; i++){
    branch br = gr->getBranch(i);
    int from = br.getFrom();
    int fn = gr->getBusNum(from);
    int to = br.getTo();
    int tn = gr->getBusNum(to);
    
    cout<<i<<": "<<from<<" - "<<to<<endl;;

    C(i,fn)=1;
    C(i,tn)=-1;

    double x = br.getX();
    double tap = br.getTap();
    double btmp;
    if (tap == 0) btmp = 1/x;
    else btmp = 1/(x*tap);
    Bff(i,i)=btmp;

    //IGNORING PHASE SHIFT FOR NOW
  }
  
  
  cout<<"C: "<<sp_mat(C)<<endl;

  mat Bf(Nl,Nb);
  mat Bbus(Nb,Nb);

  Bf=Bff*C;
  Bbus=trans(C)*Bf;
  
  cout<<"Bf: "<<sp_mat(Bf)<<endl;
  cout<<"Bbus: "<<sp_mat(Bbus)<<endl;

  for(int i=0;i<Nb;i++){
    bus bi = gr->getBus(i);
    int type=bi.getType();
    if(type==3) slack=i;
  }
  cout<<"Base Slack Bus: "<<slack<<endl;
 

  mat H(Nl,Nb-1);
  mat Hp(Nl,Nb-1);

  Hp = Bf.submat(0,1,Nl-1,Nb-1)*inv(Bbus.submat(1,1,Nb-1,Nb-1));
  H = Hp;
  H.insert_cols(0,1);
  cout<<"H calculated using inverse of Bbus matrix"<<endl;
  _H=H;

}

vec gridcalc::getDelG(del_g dg){
  int Nb=_gr->numBuses();
  vec delg(Nb,fill::zeros);

  for(int i=0;i<Nb;i++){
    delg(i)=dg.getDemand(i);
  }

  return delg;
}

vec gridcalc::convert(IloNumArray in){
  int N = in.getSize();
  vec out(N);
  for(int i=0;i<N;i++){
    out(i)=in[i];
  }
  return out;
}

vec gridcalc::getDelF(vec delg, vec slack) {
  int Nb=_gr->numBuses();
  int Nl=_gr->numBranches();
  mat Hw(_H);
  vec vt(Nl);
  vt=_H*slack;
  //  vt.t().print("Vt: ");
  //  cout<<"Vt: Cols "<<vt.n_cols<<", Rows "<<vt.n_rows<<endl;
  for(int k=0;k<Nb;k++){
    for(int i=0;i<Nl;i++){
      if(abs(vt(i,0))>=0.0000005){
	//	cout<<H(i,k)<<" - "<<vt(i,0)<<endl;
	Hw(i,k)=_H(i,k)-vt(i,0);
      }
    }
  }
  vec del_f = Hw*delg;

  return del_f;

}

mat gridcalc::getHw(vec slack){
 int Nb=_gr->numBuses();
  int Nl=_gr->numBranches();
  mat Hw(_H);
  vec vt(Nl);
  vt=_H*slack;
  //  vt.t().print("Vt: ");
  //  cout<<"Vt: Cols "<<vt.n_cols<<", Rows "<<vt.n_rows<<endl;
  for(int k=0;k<Nb;k++){
    for(int i=0;i<Nl;i++){
      if(abs(vt(i,0))>=0.0000005){
	//	cout<<H(i,k)<<" - "<<vt(i,0)<<endl;
	Hw(i,k)=_H(i,k)-vt(i,0);
      }
    }
  }
  return Hw;
}

void gridcalc::testSlack(){
  int Nb=_gr->numBuses();
  int Nl=_gr->numBranches();
  int Ng = _gr->numGens();


  mat del(Nb,1,fill::zeros);
  del(4,0)=5;
  del(6,0)=7;
  

  srand(time(NULL));
  mat slackdist(Nb,1,fill::zeros);


  double total=0;
  for(int i=0;i<Ng;i++){
    gen g = _gr->getGen(i);
    int buscon=_gr->getBusNum(g.getBus());
    double r = ((double) rand() / (RAND_MAX));
    
    if(i==1) r=1;
    else r=0; 

    slackdist(buscon,0)=r;
    
    total=total+r;
  }


  for(int i=0;i<Ng;i++){
    gen g = _gr->getGen(i);
    int buscon=_gr->getBusNum(g.getBus());
    double normalized = slackdist(buscon,0)/total;
    slackdist(buscon,0)=normalized;
  }    
  
  cout<<"Distribute slack from base bus to slack distribution"<<endl;
  slackdist.t().print("Slack: ");

  mat Hw(_H);
  mat vt(Nl,1);
  vt=_H*slackdist;
  vt.t().print("Vt: ");
  cout<<"Vt: Cols "<<vt.n_cols<<", Rows "<<vt.n_rows<<endl;
  for(int k=0;k<Nb;k++){
    for(int i=0;i<Nl;i++){
      if(abs(vt(i,0))>=0.0000005){
	//	cout<<H(i,k)<<" - "<<vt(i,0)<<endl;
	Hw(i,k)=_H(i,k)-vt(i,0);
      }
    }
  }

  //WRONG
  /*  for(int i=0;i<Nb;i++){
    if(abs(vt(i,0)*del(i,0))>=0.000007){
      cout<<i<<": "<<vt(i,0)*del(i,0)<<endl;
      cout<<vt(i,0)<<"=";
      for(int j=0;j<Ng;j++){
	gen g = _gr->getGen(j);
	int buscon=_gr->getBusNum(g.getBus());
	cout<<H(i,buscon)<<"*"<<slackdist(buscon,0)<<endl; USEFUL INFORMATION
      }
    }
    } */
  //END WRONG

  vec del_f = Hw*del;

  del_f.t().print("delta f (shift factor): ");

  cout<<sum(del_f)<<endl;

  cout<<12*sum(vt)<<endl;
  cout<<(5*sum(_H.col(4))+7*sum(_H.col(6)))<<endl;

  cout<<(5*sum(_H.col(4) - vt) + 7*sum(_H.col(6)-vt))<<endl;

  
}


void gridcalc::test(){

  del_g dg(_gr);
  dg.addDemand(4,5);
  dg.addDemand(6,7);
  
  cout<<"dg: "<< getDelG(dg).t()<<endl;

  int Nb=_gr->numBuses();
  vec delg(Nb,fill::zeros);
  vec slackdist(Nb,fill::zeros);
  delg(4)=5; delg(6)=7;
  slackdist(1)=1;

  vec del_f =  getDelF(delg,slackdist);


  igrid ig(_gr);
  ig.addCost();
  rgrid * rg;
  rgrid * rg2;

  rg = ig.solveModel();
  IloNumArray g_nom=rg->getG();
  cout<<"Before"<<endl;
  cout<<"F: "<<rg->getF()<<endl;
  cout<<"G: "<<g_nom<<endl;
  

    
  IloNumArray slack(IloEnv(),Nb);
  for(int i=0;i<Nb;i++) slack[i]=0;
  slack[1]=1;
  
  ig.addSlack(g_nom,slack);

  ig.modGrid(dg);

  rg2 = ig.solveModel();
  cout<<"After"<<endl;
  cout<<"F: "<<rg2->getF()<<endl;
  cout<<"G: "<<rg2->getG()<<endl;

  cout<<"Total Demand: "<<_gr->getTotalDemand()<<endl;
  cout<<"Total Gen: "<<IloSum(rg2->getG())<<endl;
  
  cout<<"delF (sim): "<<endl;
  for(int i=0;i<_gr->numBranches();i++)
    cout<<" "<<(rg->getF()[i] - rg2->getF()[i])<<"   ";

  cout<<"\ndelF (sim) - delF (shift factor)"<<endl;
  double totalerror=0;
  for(int i=0;i<_gr->numBranches();i++){
    cout<<" "<<(rg->getF()[i] - rg2->getF()[i] - del_f(i))<<"   ";
    totalerror=totalerror+(rg->getF()[i] - rg2->getF()[i] - del_f(i));
  }
  cout<<"\n";
  cout<<"Total Error: "<<totalerror<<endl;

  ranvar rv;
  rv.testRV();  

}
