#include <vector>
#include <iostream>
#include <ctime>
#include <fstream>
#include <omp.h>
#include "ult.h"
#include "basic_class.h"



int main() {
    
	unsigned int N = 1e6;   // particles per bunch.
	unsigned int N_bnches_p_train = 1; // bunches per train.
	unsigned int N_trns = 1; // number of trains
	double dly = 1/704e6;
	double gp = 0/704e6;
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
	bunch b1 = bunch(dist_type, spc0, mmt0, spcSig, mmtSig, N,0.0);
	double time = omp_get_wtime() - start_time;
	std::cout << "Time spend on bunch initialzation (openMP):" << time * 1000 << " ms" << std::endl;
    /*
	// Try initialize a beam with PARMELA output
	start_time = omp_get_wtime();
	std::string parmelaData = "Dist_in_defl_cavity.txt";
	beam bm1 = beam(parmelaData,704e6,N_bnches);
	time = omp_get_wtime() - start_time;
	std::cout << "Time spend on initialize a beam with PARMELA output file:" << time * 1000 << " ms" << std::endl;
	*/
	
	// Try initialize a beam

	start_time = omp_get_wtime();

	beam bm1 = beam(N,N_bnches_p_train, N_trns,dly,gp);

	time = omp_get_wtime() - start_time;
	std::cout << "Time spend on initialize a beam: " << time * 1000 << " ms" << std::endl;
	/*
	// Try initialize a bunch from a PARMELA ouput file
	start_time = omp_get_wtime();
	std::string parmelaData = "Dist00.txt";
	bunch bnch2 = bunch(parmelaData,704e6);
	time = omp_get_wtime() - start_time;
	std::cout << "Time spend on initialize a bunch from PARMELA file:" << time * 1000 << " ms" << std::endl;
	*/
	
	/*
	start_time = omp_get_wtime();
	std::string path = "data1.txt";
	bm1.bnch.dump_to_file(path);
	time = omp_get_wtime() - start_time;
	std::cout << "Time spend on dumping the bunch to file (serial):" << time * 1000 << " ms" << std::endl;
	*/

	// Try initialize cavity, ring and drift space.
	cavity cvt = cavity();
	cavity cvt2 = cavity();
	cvt.frq[0] = 704e6;//78239.554893733875*10;
	cvt.V0zR[0] = 0;
	cvt.V0zI[0] = 0;//-1e5;
	cvt.phi[0] = 0;
	cvt.kz[0] = 1.10584061406e11*2;// 1/4*(w*R/Q)
	cvt.tau_invert[0] = 1/4.521447e-5; // 2 Q/w
	
	cvt2.frq[0] = 78239.554893733875*20;
	cvt2.V0zR[0] = 0;
	cvt2.V0zI[0] = 0;
	cvt2.phi[0] = 0;
	
	ring rng = ring();
	drift_space dft = drift_space(1.73);
    std::cout<<"Initialize cavity and ring successfully."<<std::endl;
    
	// Try to iterate all particles in a bunch to update and apply the kicks

	std::string path = "data";
	unsigned int N_turns = 1;
	bm1.bnches[0].dump_to_file(path+std::to_string(0)+std::to_string(0)+"before");
	std::vector<double> temp_t;
	std::vector<double> temp_pz;
	std::vector<double> temp_vc;
	//std::cout<<bm1.bnches[0].p0<<std::endl;
	start_time = omp_get_wtime();
	for (unsigned int i = 0; i<N_turns; ++i) {
	    rng.update_f0(bm1.bnches[0]);
	//    cvt.frq[0]=rng.f0*360;
	//    cvt2.frq[0]=rng.f0*720;
    //    std::cout<<"Beam energy: "<<(bm1.bnches[0].gamma0-60.0)*bm1.bnches[0].me*c*c/bm1.bnches[0].qe<<std::endl;
    //    cvt.V0zI[0] += 1e5/N_turns;
    //    cvt2.V0zI[0] +=2e5/N_turns;
	    for (unsigned int j = 0;j< N_bnches_p_train*N_trns;++j){
	        bm1.bnches[j].sort();
    //	    std::cout<<"Calculating wake..."<<std::endl;
            cvt.wake_Naive(bm1,bm1.bnches[j]);
	//	    cvt.update_coord(bm1.bnches[j]);
	//	    cvt2.wake_Naive(bm1,bm1.bnches[j]);
	//	    cvt2.update_coord(bm1.bnches[j]);
		    rng.update_coord(bm1.bnches[j]);
	//        dft.update_coord(bm1.bnches[j]);
	//        temp_t.push_back(bm1.bnches[0].t[0]*2*pi*cvt.frq[0]);
	        temp_t.push_back(bm1.bnches[j].t[N-1]*1e9);
            temp_pz.push_back((bm1.bnches[0].px[0]/bm1.bnches[0].p0-1));
            temp_vc.push_back(sqrt(cvt.VzTauR[0][bm1.bnches[0].x.size()]*cvt.VzTauR[0][bm1.bnches[0].x.size()]+cvt.VzTauI[0][bm1.bnches[0].x.size()]*cvt.VzTauI[0][bm1.bnches[0].x.size()]));
            if ((i*N_bnches_p_train*N_trns+j)%(N_bnches_p_train*N_trns*N_turns/10)==0){
	            std::cout<<int(float(i*N_bnches_p_train*N_trns+j)/float(N_bnches_p_train*N_trns*N_turns)*100)<<"%..."<<std::endl;
	        }
	    }
	}
	time = omp_get_wtime() - start_time;
	std::cout << "Total time spend on tracking: "<< time*1000<<" ms."<<std::endl;
	std::cout << "Time spend on finish one round of kick:" << time * 1000 / N_turns << " ms." << std::endl;
	std::string cavity_voltage = "cavity";
	cvt.dump_voltage(cavity_voltage);
	std::string t = "tempT";
	std::string pz = "tempPz";
	std::string vc = "tempVc";
	dump_to_file (t,temp_t);
	dump_to_file(pz,temp_pz);
	dump_to_file (vc, temp_vc);
    bm1.bnches[0].dump_to_file(path+std::to_string(0)+std::to_string(0)+"after");
	
	/*
	start_time = omp_get_wtime();
	path = "data2.txt";
	bm1.bnch.dump_to_file(path);
	time = omp_get_wtime() - start_time;
	std::cout << "Time spend on dumping the bunch to file (serial):" << time * 1000 << " ms" << std::endl;
	*/

	return 0;
}
