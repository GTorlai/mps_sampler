#include "itensor/all.h"
#include <vector>
#include <fstream>
#include <random>
#include <string>
#include <iomanip>
//#include "sampler.h"
#include "dmrg.h"
//#include "parameters.h"

using namespace itensor;
    

int main(int argc, char* argv[]) {

  int nsites   = 4;
  int nsamples = 10000;
  std::string bc = "obc";

  DMRG dmrg(nsites);
  dmrg.TransverseFieldIsing(1.0,bc);
  dmrg.InitializeRandom();
  dmrg.Run();
  dmrg.PrintFullWavefunction();
  



    //Parameters par;
    //par.ReadParameters(argc,argv);
    //par.PrintParameters(); 
    //int N = par.N_;
    //double h = par.h_;

    //DMRG solver(N,h); 
    //solver.BuildHamiltonian();
    //solver.InitializeMPS();
    //solver.run_dmrg();
    //MPS psi = solver.GetPsi();
    //printfln("\nGround State Energy = %.10f",solver.gs_energy_/float(N));

    //Sampler sampler(N,h,psi);
    //
    //
    //sampler.get_partial_tensors();
    //for(int k=0; k<par.ns_; k++){
    //    std::cout<< k << std::endl;
    //    sampler.sample();
    //}
    
    return 0;
} 
