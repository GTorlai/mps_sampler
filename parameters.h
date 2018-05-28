#ifndef __PARAMETERS_H
#define __PARAMETERS_H
#include <stdio.h>
#include <stdlib.h>


// Parameter Class
class Parameters{

public:
   
    int N_;             // Number of hidden units 
    double h_;             // Number of visible units
    int ns_;   // Width of normal distribution for initial weights
   
    // Constructor
    Parameters() {
        // Default values
        N_ = 4;
        h_ = 1.0;
        ns_ = 1000;
    }
    
    // Read parameters from the command line
    void ReadParameters(int argc,char** argv){
        std::string flag;
        
        flag = "-N";
        for(int i=1;i<argc;i++){
            if(flag==argv[i]) N_=atoi(argv[i+1]);
        }
        flag = "-ns";
        for(int i=1;i<argc;i++){
            if(flag==argv[i]) ns_=atoi(argv[i+1]);
        }
        flag = "-h";
        for(int i=1;i<argc;i++){
            if(flag==argv[i]) h_=double(atof(argv[i+1]));
        }
    }
    
    // Print the parameters
    void PrintParameters(){
        std::cout << "\nExact sampling of Matrix Product States\n\n";
        std::cout << " Number of spins: " << N_ << std::endl;
        std::cout << " Magnetic field: " << h_<< std::endl;
        std::cout << std::endl<<std::endl;
    }
};

#endif
