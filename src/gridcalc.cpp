#include "gridcalc.h"

gridcalc::gridcalc(grid * gr) { 
  cout<<"gridcalc constructor"<<endl;
  cout<<"create useful objects"<<endl;

  _gr=gr;


  int Nb=gr->numBuses();
  int Nl=gr->numBranches();

  mat C(Nl,Nb,fill::zeros);
  mat Bff(Nl,Nl,fill::zeros);

  //  cout<<"Line: To bus - From bus  (indicies)"<<endl;
  for(int i =0; i<Nl; i++){
    branch br = gr->getBranch(i);
    int from = br.getFrom();
    int fn = gr->getBusNum(from);
    int to = br.getTo();
    int tn = gr->getBusNum(to);
    
    //    cout<<i<<": "<<from<<" - "<<to<<endl;;

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
  _C=C;
  //cout<<"C: "<<sp_mat(C.col(3))<<endl;

  mat Bf(Nl,Nb);
  mat Bbus(Nb,Nb);

  Bf=Bff*C;
  Bbus=trans(C)*Bf;
  
  //  cout<<"Bf: "<<sp_mat(Bf.col(3))<<endl;
  //  cout<<"Bbus: "<<sp_mat(Bbus)<<endl;

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
  vec del_f = getHw(slack)*delg;

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

mat gridcalc::getL(mat Hw){
  mat H = Hw*trans(_C);
  mat L = H;
  
  int Nl=_gr->numBranches();
  for(int i=0;i<Nl;i++){
    for(int j=0;j<Nl;j++){
      if(i!=j)      L(i,j)=H(i,j)/(1-H(j,j));
      else L(i,j)=-1;
    }
  }

  return L;  
}

vec gridcalc::risk(vec f,vec varf,double L, double p, double pc){
  int N=f.n_elem;
  vec z(N,fill::zeros);
  ranvar rv;
  for(int i=0;i<N;i++){
    double U = _gr->getBranch(i).getRateA();
    z(i)=rv.anaProb(L,p,pc,abs(f(i))/U,sqrt(varf(i))/U);
  }

  return z;
}
