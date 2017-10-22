#include <vector>
#include <iostream>
#include <ctime>
#include <fstream>
#include "basic_classes.h"
#include <omp.h>
bool testopenMP(){
    unsigned int N = 1e8;
    clock_t begin = clock();
    std::vector<double> p1;
    for(int i = 0;i<N;++i){
        p1.push_back(i/2.5);
    }
    clock_t end = clock(); 
    double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC*1000.0;
    std::cout<<"Elapsed time with seriel initialzation:"<< elapsed_secs<<" ms"<<std::endl;

    begin = clock();
    std::vector<double> p2(N);
#pragma omp parallel for
    for(int i = 0;i<N;++i){
        p2[i]=(i/2.5);
    }
    end = clock(); 
    elapsed_secs = double(end-begin)/CLOCKS_PER_SEC*1000.0;
    std::cout<<"Elapsed time with openMP initialzation:"<< elapsed_secs<<" ms"<<std::endl;
};

int main(){
    unsigned int N = 1e6;
    std::vector<double> kick{1,2,4};
    std::cout.precision(17);
// Debug
/*
    std::vector<particle> p1;

    clock_t begin = clock();
    for( int i = 0;i<N;++i){
        p1.push_back(particle(1,2,3,100,200,300));
    }
    for( int i = 0;i<N;++i){
        //std::cout<<i<<std::endl;
        //p1[i].display();
        p1[i].move(kick);
        //p1[i].display();
    }
    clock_t end = clock();
    double elapsed_secs = double(end-begin)/CLOCKS_PER_SEC*1000.0;
    std::cout<<"Elapsed time with seriel initialzation:"<< elapsed_secs<<" ms"<<std::endl;
    p1[0].display();

    std::vector<particle> p2(N,particle());

    begin = clock();
#pragma omp parallel for
    for( int i = 0;i<N;++i){
        p2[i] = particle(2,3,4,200,300,400);
    }
#pragma omp parallel for
    for( int i = 0;i<N;++i){
        //std::cout<<i<<std::endl;
        //p1[i].display();
        p2[i].move(kick);
        //p1[i].display();
    }
    end = clock();
    elapsed_secs = double(end-begin)/CLOCKS_PER_SEC*1000.0;
    std::cout<<"Elapsed time with openMP initialzation:"<< elapsed_secs<<" ms"<<std::endl;
    p2[0].display();
*/
    //testopenMP();
    
// Initialize a bunch
    std::cout<<"Start initializing the bunch."<<std::endl;
    std::cout<<"Number of particles:"<<N<<std::endl;
    double start_time = omp_get_wtime();
    std::vector<double> spc0{0,0,0};
    std::vector<double> spcSig{2,2,2};
    std::vector<double> mmt0{0,0,0};
    std::vector<double> mmtSig{2,2,2};
    unsigned int dist_type = 1;
    bunch b1=bunch(dist_type,spc0,mmt0,spcSig,mmtSig,N);;
    double time = omp_get_wtime() - start_time;
    std::cout<<"Time spend on bunch initialzation (openMP):"<< time*1000<<" ms"<<std::endl;

    /*
    start_time = omp_get_wtime();
    std::string path = "data.txt";
    b1.dump_to_file(path);
    time = omp_get_wtime() - start_time;
    std::cout<<"Time spend on dumping the bunch to file (serial):"<< time*1000<<" ms"<<std::endl;
    */
    
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
 
#pragma omp parallel for
    for (int i = 0;i<N;++i){
        cvt.update_coord(bm1.bnch.ptcs[i]);
        rng.update_coord(bm1.bnch.ptcs[i]);
    }
    time = omp_get_wtime() - start_time;
    std::cout<<"Time spend on finish one round of kick:"<< time*1000<<" ms"<<std::endl;
    
    start_time = omp_get_wtime();
    path = "data2.txt";
    bm1.bnch.dump_to_file(path);
    time = omp_get_wtime() - start_time;
    std::cout<<"Time spend on dumping the bunch to file (serial):"<< time*1000<<" ms"<<std::endl;
    
    
    return 1;
}



























