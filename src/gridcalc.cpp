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

  H = Bf.submat(0,1,Nl-1,Nb-1)*inv(Bbus.submat(1,1,Nb-1,Nb-1));

  mat del(Nb-1,1);
  del.fill(0);
  del(3,0)=10;
  

  cout<<H<<endl;
  
  cout<<H*del<<endl;

  cout<<sum(H*del)<<endl;


}
