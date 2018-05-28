#include "itensor/all.h"
#include <vector>
#include <fstream>
#include <random>
#include <string>
#include <iomanip>
#include "sampler.h"
#include "dmrg.h"
#include "parameters.h"

using namespace itensor;
    

int main(int argc, char* argv[]) {

    Parameters par;
    par.ReadParameters(argc,argv);
    par.PrintParameters(); 
    int N = par.N_;
    double h = par.h_;

    DMRG solver(N,h); 
    solver.BuildHamiltonian();
    solver.InitializeMPS();
    solver.run_dmrg();
    MPS psi = solver.GetPsi();
    printfln("\nGround State Energy = %.10f",solver.gs_energy_/float(N));

    Sampler sampler(N,h,psi);
    
    
    sampler.get_partial_tensors();
    for(int k=0; k<ns_; k++){
        std::cout<< k << std::endl;
        sampler.sample();
    }
    
    return 0;
} 
