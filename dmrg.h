#ifndef _DMRG_H_
#define _DMRG_H_

#include <string>
#include <iostream>

using namespace itensor;

class DMRG {
    
  int nsites_;
  MPO Hamiltonian_;
  SiteSet sites_;
  std::mt19937 rgen_;
  MPS psi_;
  
public:
    
  DMRG(int nsites):nsites_(nsites){} 
  
  inline MPS GetPsi() {
      return psi_;
  }
  
  //******** HAMILTONIANS *********//
  void TransverseFieldIsing(double h,std::string bc){
    sites_ = SpinHalf(nsites_);
    auto ampo = AutoMPO(sites_);
    for(int j = 1; j < nsites_; ++j) {
      ampo += -4.0,"Sz",j,"Sz",j+1;
      ampo += -2.0*h,"Sx",j;
    }
    ampo += -2.0*h,"Sx",nsites_;
    if (bc=="pbc"){
      ampo += -4.0,"Sz",nsites_,"Sz",0;
    }
    Hamiltonian_ = MPO(ampo);
  }

  void Heisenberg(std::string bc){
    sites_ = SpinHalf(nsites_);
    auto ampo = AutoMPO(sites_);
    for(int j = 1; j < nsites_; ++j) {
      ampo += 0.5,"S+",j,"S-",j+1;
      ampo += 0.5,"S-",j,"S+",j+1;
      ampo += 1.0,"Sz",j,"Sz",j+1;
    }
    if (bc=="pbc"){
      ampo += 0.5,"S+",nsites_,"S-",1;
      ampo += 0.5,"S-",nsites_,"S+",1;
      ampo += 1.0,"Sz",nsites_,"Sz",1;
    }
    Hamiltonian_ = MPO(ampo);
  }

  void InitializeRandom(){
    auto state = InitState(sites_);
    std::uniform_int_distribution<int> distribution(0,1);
    for(int i = 1; i <= nsites_; ++i) {
      auto ran  = distribution(rgen_);
      if (ran == 1)
        state.set(i,"Up");
      else
        state.set(i,"Dn");
    }
    psi_ = MPS(state);
  }

  void InitializeNeel(){
    auto state = InitState(sites_);
    std::uniform_int_distribution<int> distribution(0,1);
    for(int i = 1; i <= nsites_; ++i) {
      if(i%2 == 1)
        state.set(i,"Up");
      else
        state.set(i,"Dn");
      }
    psi_ = MPS(state);
  }


  void Run(){
      auto sweeps = Sweeps(10);
      sweeps.maxm() = 10,20,100,100,200;
      sweeps.cutoff() = 1E-10;
      sweeps.niter() = 2;
      sweeps.noise() = 1E-7,1E-8,0.0;
  
      // Begin the DMRG calculation
      auto energy = dmrg(psi_,Hamiltonian_,sweeps,"Quiet");
      //gs_energy_ = energy; 
      printfln("\nGround State Energy = %.10f",energy/float(nsites_));
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

  double CollapseWavefunction(std::vector<int> & sigma){
    ITensor psi_sigma;
    ITensor s;
    s = BasisState(sigma[0],psi_.A(1));
    psi_sigma = psi_.A(1);
    psi_sigma *= s;

    for (int j=2;j<=nsites_;j++){
      s = BasisState(sigma[j-1],psi_.A(j));
      psi_sigma *= psi_.A(j);
      psi_sigma *= s;
    }
    return psi_sigma.real();
  }

  void PrintFullWavefunction(){
    if(nsites_>10){
      std::cout<<"System size too large!"<<std::endl;
      exit(0);
    }
    std::bitset<10> bit;
    std::vector<int> sigma;
    for(int i=0; i<1<<nsites_;i++){
      bit = i;
      sigma.clear();
      for(int j=0;j<nsites_;j++){
        sigma.push_back(bit[nsites_-1-j]);
      }
      auto psi_sigma = CollapseWavefunction(sigma); 
      std::cout<<psi_sigma<<std::endl;
      //std::cout<<psi_sigma.cplx().real()<<"  "<<psi_sigma.cplx().imag()<<std::endl;
    }
    std::cout<<std::endl<<std::endl;
  }



};

#endif 
