#ifndef _SAMPLER_H_
#define _SAMPLER_H_

#include <string>
#include <iostream>
#include <vector>
#include <iomanip>

using namespace itensor;

class Sampler {

private:
  
  int N_;
  int nsamples_;
  MPS psi_;
  std::mt19937 rgen_;

  std::vector<int> state_;
  std::vector<ITensor> partial_tensors_;
  std::vector<ITensor> trail_tensors_;
 
  std::vector<std::vector<std::string> > bases_;
  std::complex<double> I_;

public:
    
  //Default constructor
  Sampler(int nsites,int nsamples,MPS & psi):N_(nsites),nsamples_(nsamples),psi_(psi),I_(0,1) {
    state_.assign(N_,0);
  } 
  
  void GetPartialTensors(MPS &psi){
    partial_tensors_.clear();
    ITensor pt; //partial tensor tmp variable
    
    pt = psi_.A(1);
    pt *= dag(prime(psi.A(1),Link));
    partial_tensors_.push_back(pt);
    for (int j=2; j<N_; j++){
        pt *= psi.A(j);
        pt *= dag(prime(psi.A(j),Link));
        partial_tensors_.push_back(pt);
    }
  }

  void Sample(){
    std::vector<std::string> basis;
    int basis_id;
    //If sampling in the computational basis
    std::uniform_int_distribution<int> dist(0,bases_.size()-1);
    for(int j=0;j<nsamples_;j++){
      basis_id = dist(rgen_);
      basis = bases_[basis_id];
      MPS psi = RotateMPS(basis);
      OneSample(psi);
    }
  }

  void OneSample(MPS &psi){

    GetPartialTensors(psi); 
    double prob;
    double joint_prob= 1.0;
    std::uniform_real_distribution<double> distr(0,1);

    trail_tensors_.clear();
    ITensor rho;
    ITensor trail,trail_dag;
    ITensor tmp;
    
    // Sample spin at site N
    rho = partial_tensors_[N_-2];
    rho *= psi.A(N_);
    tmp = dag(prime(psi.A(N_),Site,Link));
    rho *= tmp;
    prob = rho.real(2,2);
    //std::cout<<prob<<std::endl;
    state_[N_-1] = distr(rgen_) < prob;
    if (state_[N_-1] == 1){
      joint_prob *= prob;
    }
    else {
      joint_prob *=(1.0-prob);
    }
    
    // Sample spin at site N-1
    rho = partial_tensors_[N_-3];
    rho = rho * psi.A(N_-1);
    tmp = dag(prime(psi.A(N_-1),Site,Link));
    rho = rho * tmp;
    
    trail  = psi_.A(N_);
    trail *=BasisState(state_[N_-1],psi_.A(N_));
    trail *= dag(prime(psi.A(N_),Site,Link));
    trail *= prime(BasisState(state_[N_-1],psi.A(N_)),Site);
    rho   *= trail;
    rho    = rho / joint_prob;

    prob = rho.real(2,2);
    state_[N_-2] = distr(rgen_) < prob;
    if (state_[N_-2] == 1){
      joint_prob *= prob;
    }
    else {
      joint_prob *=(1.0-prob);
    }
    //std::cout<<"state="<<state_[N_-1]<<" prob  = "<<prob<<std::endl; 
    // Sample the other spins
    for(int j= N_-2; j>0; j--){
      if(j>1){
        rho = partial_tensors_[j-2];
        rho = rho * psi.A(j);
      }
      else {
        rho= psi.A(j);
      }
      tmp = dag(prime(psi.A(j),Site,Link));
      rho = rho * tmp;
      trail *= psi.A(j+1);
      trail *= BasisState(state_[j],psi.A(j+1));
      trail *= dag(prime(psi.A(j+1),Site,Link));
      trail *= prime(BasisState(state_[j],psi.A(j+1)),Site);
      rho   *= trail;
      rho = rho / joint_prob;
      prob = rho.real(2,2);
      state_[j-1] = distr(rgen_) < prob;
      if (state_[j-1] == 1){
        joint_prob *= prob;
      }
      else {
        joint_prob *=(1.0-prob);
      }
      //std::cout<<"state="<<state_[j]<<" prob  = "<<prob<<std::endl; 
    }
  }
   
  ITensor BasisState(int & s_value, const ITensor & A){
    auto index = findtype(A,Site);
    ITensor s(index);
    if (s_value == 0){
      s.set(index(1),1.0);
      s.set(index(2),0.0);
    }
    else {
      s.set(index(1),0.0);
      s.set(index(2),1.0);
    }
    return s;
  }

  std::complex<double> CollapseWavefunction(MPS &psi,std::vector<int> & sigma){
    ITensor psi_sigma;
    ITensor s;
    s = BasisState(sigma[0],psi.A(1));
    psi_sigma = psi.A(1);
    psi_sigma *= s;

    for (int j=2;j<=N_;j++){
      s = BasisState(sigma[j-1],psi.A(j));
      psi_sigma *= psi.A(j);
      psi_sigma *= s;
    }
    return psi_sigma.cplx();
  }
   
  ITensor PauliX(const ITensor &A){
    auto s = findtype(A,Site);
    auto sp = prime(s);
    ITensor pauliX(s,sp);
    pauliX.set(s(1),sp(1),1.0/std::sqrt(2)); 
    pauliX.set(s(1),sp(2),1.0/std::sqrt(2)); 
    pauliX.set(s(2),sp(1),1.0/std::sqrt(2)); 
    pauliX.set(s(2),sp(2),-1.0/std::sqrt(2)); 
    return pauliX;
  }

  ITensor PauliY(const ITensor &A){
    auto s = findtype(A,Site);
    auto sp = prime(s);
    ITensor pauliY(s,sp);
    //pauliY.set(s(1),sp(1),1.0/std::sqrt(2)); 
    //pauliY.set(s(1),sp(2),-I_/std::sqrt(2)); 
    //pauliY.set(s(2),sp(1),1.0/std::sqrt(2)); 
    //pauliY.set(s(2),sp(2),I_/std::sqrt(2)); 
    pauliY.set(s(1),sp(1),1.0/std::sqrt(2)); 
    pauliY.set(s(1),sp(2),1.0/std::sqrt(2)); 
    pauliY.set(s(2),sp(1),I_/std::sqrt(2)); 
    pauliY.set(s(2),sp(2),-I_/std::sqrt(2)); 
    return pauliY;
  }

  MPS RotateMPS(std::vector<std::string> &basis){
 
    auto psi = psi_;
    for(int j=1;j<=N_;j++){
      if (basis[j-1] != "I"){
        psi.position(j);
        auto wf = psi.A(j);
        ITensor rotation;
        if (basis[j-1] == "X")
          rotation = PauliX(wf);
        if (basis[j-1] == std::string("Y"))
          rotation = dag(PauliY(wf));
        //prime(rotation,rotation.index(2),1);
        wf *= rotation;
        wf.noprime();
        psi.setA(j,wf);
      }
    }
    return psi;
  }

  void LoadBases(std::ifstream &fin){
    int num_bases;
    fin >> num_bases;
    std::cout<<num_bases<<std::endl;
    bases_.resize(num_bases,std::vector<std::string>(N_));
    for (int i=0;i<num_bases;i++){
      for(int j=0;j<N_;j++){
        fin >> bases_[i][j];
      }
    }
  }

  void PrintFullWavefunction(MPS &psi){
    if(N_>10){
      std::cout<<"System size too large!"<<std::endl;
      exit(0);
    }
    std::bitset<10> bit;
    std::vector<int> sigma;
    for(int i=0; i<1<<N_;i++){
      bit = i;
      sigma.clear();
      for(int j=0;j<N_;j++){
        sigma.push_back(bit[N_-1-j]);
      }
      auto psi_sigma = CollapseWavefunction(psi,sigma); 
      std::cout<<psi_sigma<<std::endl;
    }
    std::cout<<std::endl;
  }

  void TestRotations(){

    for(int j=0;j<bases_.size();j++){
      std::cout<<"Basis ";
      for (int k=0;k<N_;k++){
        std::cout<<bases_[j][k];
      }
      std::cout<<std::endl<<std::endl;
      MPS psi = RotateMPS(bases_[j]);
      PrintFullWavefunction(psi);
      std::cout<<std::endl;
    }
  }

  void TestSampler(){
      
    int ns = 100000;
    std::bitset<10> bit;
    std::vector<int> conf;
    std::vector<std::complex<double> > psi;
    std::vector<int> hist;
    std::vector<double> prob_approx;
    int ind;
    hist.assign(1<<N_,0);
    for(int i=0; i<1<<N_;i++){
      bit = i;
      conf.clear();
      for (int j=0;j<N_;j++){
        conf.push_back(bit[N_-1-j]);
      }
      psi.push_back(CollapseWavefunction(psi_,conf));
    }
    std::cout <<std::endl<<std::endl; 
    for (int i=0; i<ns;i++){
      OneSample(psi_);
      ind = 0;
      for(int j=0;j<N_;j++){
        ind += state_[N_-j-1]*pow(2,j);    
      }
      hist[ind]++;
    }
    for(int i=0;i<1<<N_;i++){
      bit = i;
      conf.clear();
      for (int j=0;j<N_;j++){
        conf.push_back(bit[N_-1-j]);
        std::cout<<conf[j]<< " ";
      }
      std::cout<< "      |Psi|^2 = " << std::real(psi[i]*std::conj(psi[i]));
      std::cout<< " \tSampled = " <<double(hist[i])/double(ns) << std::endl;
    }
  }
  void PrintState(){
    for (int j=0;j<N_;j++){
      std::cout<<state_[j];
    }
    std::cout<<std::endl<<std::endl;
  }
};

#endif 
