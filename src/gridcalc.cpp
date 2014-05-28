#include "gridcalc.h"

gridcalc::gridcalc(grid * gr) { 
  cout<<"gridcalc constructor"<<endl;
  cout<<"create useful objects"<<endl;

  _gr=gr;


  int Nb=gr->numBuses();
  int Nl=gr->numBranches();

  cout<<*gr<<endl;

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

  mat del(Nb,1,fill::zeros);
  del(4,0)=5;
  del(6,0)=7;
  

  srand(time(NULL));
  mat slackdist(Nb,1,fill::zeros);

  int Ng = gr->numGens();
  double total=0;
  for(int i=0;i<Ng;i++){
    gen g = gr->getGen(i);
    int buscon=gr->getBusNum(g.getBus());
    double r = ((double) rand() / (RAND_MAX));
    
    if(i==1) r=1;
    else r=0; 

    slackdist(buscon,0)=r;
    
    total=total+r;
  }


  for(int i=0;i<Ng;i++){
    gen g = gr->getGen(i);
    int buscon=gr->getBusNum(g.getBus());
    double normalized = slackdist(buscon,0)/total;
    slackdist(buscon,0)=normalized;
  }    
  
  cout<<"Distribute slack from base bus to slack distribution"<<endl;
  slackdist.t().print("Slack: ");

  mat Hw(H);
  mat vt(Nl,1);
  vt=H*slackdist;
  vt.t().print("Vt: ");
  cout<<"Vt: Cols "<<vt.n_cols<<", Rows "<<vt.n_rows<<endl;
  for(int k=0;k<Nb;k++){
    for(int i=0;i<Nl;i++){
      if(abs(vt(i,0))>=0.0000005){
	//	cout<<H(i,k)<<" - "<<vt(i,0)<<endl;
	Hw(i,k)=H(i,k)-vt(i,0);
      }
    }
  }

  //WRONG
  /*  for(int i=0;i<Nb;i++){
    if(abs(vt(i,0)*del(i,0))>=0.000007){
      cout<<i<<": "<<vt(i,0)*del(i,0)<<endl;
      cout<<vt(i,0)<<"=";
      for(int j=0;j<Ng;j++){
	gen g = gr->getGen(j);
	int buscon=gr->getBusNum(g.getBus());
	cout<<H(i,buscon)<<"*"<<slackdist(buscon,0)<<endl;
      }
    }
    } */
  //END WRONG


  vec del_f = Hw*del;

  del_f.t().print("delta f: ");
  _del_f = del_f;

  cout<<sum(del_f)<<endl;

  cout<<12*sum(vt)<<endl;
  cout<<(5*sum(H.col(4))+7*sum(H.col(6)))<<endl;

  cout<<(5*sum(H.col(4) - vt) + 7*sum(H.col(6)-vt))<<endl;

  
}
