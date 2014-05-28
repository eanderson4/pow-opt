#include "gridcalc.h"

gridcalc::gridcalc(grid * gr) { 
  cout<<"gridcalc constructor"<<endl;
  cout<<"create useful objects"<<endl;

  _gr=gr;


  int Nb=gr->numBuses();
  int Nl=gr->numBranches();
  int Nk=gr->numGens();
  double v[2*Nl];
  double bff[2*Nl];
  int c[2*Nl];
  int r[2*Nl];

  cout<<*gr<<endl;
  cout<<Nl<<endl;

  mat C(Nl,Nb);
  mat Bff(Nl,Nl);
  //  C.fill(0);
  int slack;

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
  
  
  cout<<C<<endl;

  mat Bf(Nl,Nb);
  mat Bbus(Nb,Nb);

  Bf=Bff*C;
  Bbus=trans(C)*Bf;

  cout<<Bf<<endl;
  cout<<Bbus<<endl;

  for(int i=0;i<Nb;i++){
    bus bi = gr->getBus(i);
    int type=bi.getType();
    if(type==3) slack=i;
  }
  cout<<"slack: "<<slack<<endl;
 

  mat H(Nl,Nb-1);
  mat Hp(Nl,Nb-1);

  Hp = Bf.submat(0,1,Nl-1,Nb-1)*inv(Bbus.submat(1,1,Nb-1,Nb-1));
  H = Hp;
  H.insert_cols(0,1);


  mat del(Nb,1,fill::zeros);
  //  del.fill(0);
  del(4,0)=10;
  

  cout<<H<<endl;
  
  cout<<H*del<<endl;

  cout<<sum(H*del)<<endl;

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
  

  

  cout<<slackdist<<endl;

  

  mat Hw(H);
  mat vt(Nl,1);
  vt=H*slackdist;
  cout<<vt<<endl;
  for(int k=0;k<Nb;k++){
    for(int i=0;i<Nl;i++){
      if(abs(vt(i,0))>=0.0000005){
	//	cout<<H(i,k)<<" - "<<vt(i,0)<<endl;
	Hw(i,k)=H(i,k)-vt(i,0);
      }
    }
  }
  for(int i=0;i<Nb;i++){
    if(abs(vt(i,0)*del(i,0))>=0.000007){
      cout<<i<<": "<<vt(i,0)*del(i,0)<<endl;
      cout<<vt(i,0)<<"=";
      for(int j=0;j<Ng;j++){
	gen g = gr->getGen(j);
	int buscon=gr->getBusNum(g.getBus());
	cout<<H(i,buscon)<<"*"<<slackdist(buscon,0)<<endl;
      }
    }
  }
  cout<<Hw*del<<endl;
  cout<<sum(Hw*del)<<endl;


}
