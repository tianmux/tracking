#include <vector>
#include <iostream>
#include <ctime>
#include <fstream>
#include "basic_class.h"
#include <omp.h>

int main() {
	unsigned int N = 1e3;   // particles per bunch.
	unsigned int N_bnches = 1; // bunches per beam.
	std::vector<double> kick{ 1,2,4 };
	std::cout.precision(17);

	// Try Initialize a bunch
	std::cout << "Start initializing the bunch." << std::endl;
	std::cout << "Number of particles:" << N << std::endl;
	double start_time = omp_get_wtime();
	std::vector<double> spc0{ 0,0,0 };
	std::vector<double> spcSig{ 2e-3,2e-3,2e-3 };// beam size;
	std::vector<double> mmt0{ 0,0,0 };
	std::vector<double> mmtSig{ 2e-3,2e-3,2e-3 };// beam angle spread
	unsigned int dist_type = 1;
	bunch b1 = bunch(dist_type, spc0, mmt0, spcSig, mmtSig, N);
	double time = omp_get_wtime() - start_time;
	std::cout << "Time spend on bunch initialzation (openMP):" << time * 1000 << " ms" << std::endl;

	// Try initialize a beam
	start_time = omp_get_wtime();
	beam bm1 = beam(N,N_bnches);
	time = omp_get_wtime() - start_time;
	std::cout << "Time spend on initialize a beam:" << time * 1000 << " ms" << std::endl;
	/*
	start_time = omp_get_wtime();
	std::string path = "data1.txt";
	bm1.bnch.dump_to_file(path);
	time = omp_get_wtime() - start_time;
	std::cout << "Time spend on dumping the bunch to file (serial):" << time * 1000 << " ms" << std::endl;
	*/

	// Try initialize a cavity
	cavity cvt = cavity();
	cvt.frq[0] = 704e6;
	cvt.V0xR[0] = 1;
	ring rng = ring();
	drift_space dft = drift_space(3);
    std::cout<<"Initialize cavity and ring successfully."<<std::endl;
    
	// Try to iterate all particles in a bunch to update and apply the kicks
	start_time = omp_get_wtime();
	std::string path = "data";
	unsigned int N_turns = 1;
	for (unsigned int i = 0; i<N_turns; ++i) {

	    for (unsigned int j = 0;j< N_bnches;++j){
            bm1.bnches[j].dump_to_file(path+std::to_string(i)+std::to_string(j)+"before");
	        bm1.bnches[j].sort();
    //	    std::cout<<"Calculating wake..."<<std::endl;
            cvt.wake_Naive(bm1,bm1.bnches[j]);
    //	    std::cout<<"Finished with wake calculation."<<std::endl;        
		    cvt.update_coord(bm1.bnches[j]);
	//	    rng.update_coord(bm1.bnches[j]);
	        dft.update_coord(bm1.bnches[j]);
		    bm1.bnches[j].dump_to_file(path+std::to_string(i)+std::to_string(j)+"after");
		}
	}
	
	time = omp_get_wtime() - start_time;
	std::cout << "Total time spend on tracking: "<< time*1000<<" ms."<<std::endl;
	std::cout << "Time spend on finish one round of kick:" << time * 1000 / N_turns << " ms." << std::endl;
	/*
	start_time = omp_get_wtime();
	path = "data2.txt";
	bm1.bnch.dump_to_file(path);
	time = omp_get_wtime() - start_time;
	std::cout << "Time spend on dumping the bunch to file (serial):" << time * 1000 << " ms" << std::endl;
	*/

	return 0;
}
