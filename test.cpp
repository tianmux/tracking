#include <vector>
#include <iostream>
#include <ctime>
#include <fstream>
#include "basic_classes.h"
#include <omp.h>

int main(){
    unsigned int N = 1e6;
    std::vector<double> kick{1,2,4};
    std::cout.precision(17);

// Initialize a bunch
    std::cout<<"Start initializing the bunch."<<std::endl;
    std::cout<<"Number of particles:"<<N<<std::endl;
    double start_time = omp_get_wtime();
    std::vector<double> spc0{0,0,0};
    std::vector<double> spcSig{2,2,2};
    std::vector<double> mmt0{0,0,0};
    std::vector<double> mmtSig{2,2,2};
    unsigned int dist_type = 1;
    bunch b1=bunch(dist_type,spc0,mmt0,spcSig,mmtSig,N);
    double time = omp_get_wtime() - start_time;
    std::cout<<"Time spend on bunch initialzation (openMP):"<< time*1000<<" ms"<<std::endl;
    
// Try initialize a beam
    start_time = omp_get_wtime();
    beam bm1 = beam(N);
    time = omp_get_wtime() - start_time;
    std::cout<<"Time spend on initialize a beam:"<< time*1000<<" ms"<<std::endl;
    
    start_time = omp_get_wtime();
    std::string path = "data1.txt";
    bm1.bnch.dump_to_file(path);
    time = omp_get_wtime() - start_time;
    std::cout<<"Time spend on dumping the bunch to file (serial):"<< time*1000<<" ms"<<std::endl;
    
    
// Try initialize a cavity
    cavity cvt = cavity();
    ring rng = ring();
    
// Try to iterate all particles in a bunch to update and apply the kicks
    start_time = omp_get_wtime();
    unsigned int N_turns = 1000;
    for(unsigned int i = 0;i<N_turns;++i){
        cvt.update_coord(bm1.bnch);
        rng.update_coord(bm1.bnch);
    }
    time = omp_get_wtime() - start_time;
    std::cout<<"Time spend on finish one round of kick:"<< time*1000/N_turns<<" ms"<<std::endl;
    
    start_time = omp_get_wtime();
    path = "data2.txt";
    bm1.bnch.dump_to_file(path);
    time = omp_get_wtime() - start_time;
    std::cout<<"Time spend on dumping the bunch to file (serial):"<< time*1000<<" ms"<<std::endl;
    
    
    return 1;
}



























