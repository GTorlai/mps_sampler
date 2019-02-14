#include "itensor/all.h"
#include <vector>
#include <fstream>
#include <random>
#include <string>
#include <iomanip>
#include "sampler.h"
#include "dmrg.h"

using namespace itensor;
    

int main(int argc, char* argv[]) {

  int nsites   = 4;
  int nsamples = 10000;
  std::string bc = "obc";
  std::ifstream fin("bases.txt"); 

  DMRG dmrg(nsites);
  //dmrg.TransverseFieldIsing(1.0,bc);
  dmrg.Heisenberg(bc);
  dmrg.InitializeRandom();
  dmrg.Run();
  dmrg.PrintFullWavefunction();
  MPS psi = dmrg.GetWavefunction();

  Sampler sampler(nsites,nsamples,psi);
  sampler.LoadBases(fin);
  sampler.TestRotations(); 
  sampler.TestSampler();  
} 
