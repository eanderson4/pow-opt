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
    int fn = gr->getFromBus(i);
    int tn = gr->getToBus(i);

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

  cout<<"Calculating Bbus inverse"<<endl;
  Hp = Bf.submat(0,1,Nl-1,Nb-1)*inv(Bbus.submat(1,1,Nb-1,Nb-1));
  H = Hp;
  H.insert_cols(0,1);
  cout<<"H calculated using inverse of Bbus matrix"<<endl;
  _H=H;

  int Ng = gr->numGens();
  vec slackdist(Nb,1,fill::zeros);
  
  cout<<"slack find"<<endl;

  double total=0;
  for(int i=0;i<Ng;i++){
    int buscon=_gr->getGenBus(i);
    int r;
    if(i==1) r=1;
    else r=0; 
    slackdist(buscon,0)=r;
    total=total+r;
  }

  cout<<"normalize"<<endl;
  if(total!=1){
    for(int i=0;i<Ng;i++){
      int buscon=_gr->getGenBus(i);
      double normalized = slackdist(buscon)/total;
      slackdist(buscon,0)=normalized;
    } 
  }
  _slack=slackdist;
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
      if(H(j,j)<=1+.0000001 && H(j,j)>=1-.0000001) {
	L(i,j)=datum::nan;
      }
      else{
	if(i!=j)      L(i,j)=H(i,j)/(1-H(j,j));
	else L(i,j)=-1;
      }
    }
  }

  cout<<"L built"<<endl;
  return L;  
}

mat gridcalc::getCm(){
  int Nb=_gr->numBuses();
  int Ng=_gr->numGens();
  mat Cm(Nb,Ng,fill::zeros);
  for(int i=0;i<Ng;i++){
    gen g = _gr->getGen(i);
    int bus=_gr->getBusNum(g.getBus());
    Cm(bus,i)=1;
  }

  return Cm;
}

vec gridcalc::risk(vec f,vec varf,double L, double p, double pc){
  int N=f.n_elem;
  vec z(N,fill::zeros);
  ranvar rv;
  for(int i=0;i<N;i++){
    double U = _gr->getBranch(i).getRateA();
    if(varf(i)<=0.0000001)    z(i)=rv.anaProb(L,p,pc,abs(f(i))/U,0);
    else    z(i)=rv.anaProb(L,p,pc,abs(f(i))/U,sqrt(varf(i))/U);
    if(isnan(float(z(i)))){
      cout<<z(i)<<endl;
      cout<<L<<"\t"<<p<<"\t"<<pc<<"\t"<<abs(f(i))/U<<"\t"<<sqrt(varf(i))/U<<endl;
    }
  }

  return z;
}


vec gridcalc::getN1(int n, vec y0, vec g, mat Hw){
  grid * gr = _gr;
  vec yn;
  mat C = getC();
  mat L=getL(Hw);
  mat Hb=Hw*C.t();

  if(Hb(n,n)<=1-.000001 || Hb(n,n)>=1+.0000001){
    cout<<"Line "<<n<<endl;
    yn = y0[n]*L.col(n)+y0;
  }
  else{
    cout<<"Line: "<<n<<" --- HB = 1"<<endl;  //obvious island ?
    
    branch br = gr->getBranch(n);
    int from = gr->getBusNum(br.getFrom());
    int to = gr->getBusNum(br.getTo());
    int winner =-1;
    cout<<from<<" -> "<<to<<endl;
    
    if(sum(abs(C.col(from))) == 1) {
      cout<<"winner winner: "<<from<<endl;
      winner=from;
    }
    else if(sum(abs(C.col(to))) == 1) {
      cout<<"winner winner: "<<to<<endl;
      winner=to;
    }
    else {
      throw nowin;
    }
    
    bus b = _gr->getBus(winner);
    double load = b.getP()+b.getGs();
    mat Cm = getCm();
    if(sum(Cm.row(winner))>=1){
      cout<<"Have Generation"<<endl;
      vec gnode = Cm.row(winner)*g;
      yn = Hw.col(winner)*(load-gnode) + y0;
    }
    else{
      cout<<"Loadshed "<<load<<" at "<<winner<<endl;
      yn = Hw.col(winner)*load+y0;
    }
    
  }
  
  return yn;

}
