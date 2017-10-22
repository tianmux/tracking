#include <vector>
#include <iostream>
#include <random>

/*Global constants*/

double c = 3e8;
double pi = 3.1415926;

/* Class of particle*/
class particle{

    double M = 1;    /* mass number*/
    double N = 1;    /* charge number*/
    double q = N*1.60217662e-19;
    double m = M*9.10938356e-31;
    
    public:
        std::vector<double> spc;    /* The vector to store 3D spacial coordinates of the particle*/
        std::vector<double> mmt;    /* The vector to store 3D momentum coordiantes of the particle*/
        double time = 0;    /* timer for each particle*/
        
        particle(double x, double y,double z, double px, double py, double pz){
            spc.resize(3);
            mmt.resize(3);
            spc[0] = x;         
            spc[1] = y;
            spc[2] = z;         
            mmt[0] = px;         
            mmt[1] = py;         
            mmt[2] = pz;
        };
        particle(){
            spc.resize(3);
            mmt.resize(3);
            spc[0] = 0;         
            spc[1] = 0;
            spc[2] = 0;         
            mmt[0] = 0;         
            mmt[1] = 0;         
            mmt[2] = 0; 
        };


        void display(){
            std::cout<<spc[0]<<','<<spc[1]<<','<<spc[2]<<std::endl;
            std::cout<<mmt[0]<<','<<mmt[1]<<','<<mmt[2]<<std::endl;
        };
        
                
};


// Class of Bunch
class bunch{
    unsigned int Np;



    public:
        double gamma; // reference gamma of the bunch.
        std::vector<particle> ptcs;
        // Initialize with given normal distribution parameters
        bunch(unsigned int disType, 
        std::vector<double> spc0, // vector that contain centers of distribution in spacial coordinates
        std::vector<double> mmt0, // vector that contain centers of distribution in momentum coordinates
        std::vector<double> spcSig, // vector that contain sigma of distribution in spacial coordinates
        std::vector<double> mmtSig, // vector that contain sigma of distribution in momentum coordinates
        unsigned int N){
            Np = N;
            ptcs.resize(Np);
            /* So far only the Gaussian distribution is supported */
            std::default_random_engine generator; // Normal distribution generator
            std::normal_distribution<double> dist0(spc0[0],spcSig[0]); // normal distribution in x
            std::normal_distribution<double> dist1(spc0[1],spcSig[1]); // normal distribution in y
            std::normal_distribution<double> dist2(spc0[2],spcSig[2]); // normal distribution in z
            std::normal_distribution<double> dist3(mmt0[0],mmtSig[0]); // normal distribution in px
            std::normal_distribution<double> dist4(mmt0[1],mmtSig[1]); // normal distribution in py
            std::normal_distribution<double> dist5(mmt0[2],mmtSig[2]); // normal distribution in pz
            // initialize the particles in one bunch in parallel.
#pragma omp parallel for
            for(int j = 0;j<Np;++j){
            ptcs[j] = particle(dist0(generator),dist1(generator),dist2(generator),dist3(generator),dist4(generator),dist5(generator));            
            }
        };
        
        // Initialize with defualt normal distribution.
        bunch(unsigned int N){
            Np = N;
            ptcs.resize(Np);
            /* So far only the Gaussian distribution is supported */
            std::default_random_engine generator; // Normal distribution generator
            std::normal_distribution<double> dist0(0,1e-3); // normal distribution in x
            std::normal_distribution<double> dist1(0,1e-3); // normal distribution in y
            std::normal_distribution<double> dist2(0,1e-3); // normal distribution in z
            std::normal_distribution<double> dist3(0,1e-3); // normal distribution in px
            std::normal_distribution<double> dist4(0,1e-3); // normal distribution in py
            std::normal_distribution<double> dist5(0,1e-3); // normal distribution in pz
            // initialize the particles in one bunch in parallel.
#pragma omp parallel for
            for(int j = 0;j<Np;++j){
            ptcs[j] = particle(dist0(generator),dist1(generator),dist2(generator),dist3(generator),dist4(generator),dist5(generator));            
            }
        };
        
        // print the 6D coordiantes of a couple of particle, mainly for debug
        void display(){
            ptcs[0].display();
            ptcs[Np/2].display();
            ptcs[Np-1].display();
        };
        // dump the whole bunch to a file. file name is specified by 'path'.
        void dump_to_file(std::string path){
            std::filebuf fb;
            fb.open (path,std::ios::out);
            std::ostream os(&fb);
            for(int i= 0;i<Np;++i){
                os<<ptcs[i].spc[0]<<","<<ptcs[i].spc[1]<<","<<ptcs[i].spc[2]<<';'<<ptcs[i].mmt[0]<<","<<ptcs[i].mmt[1]<<","<<ptcs[i].mmt[2]<<'\n';
            }
        };
        
};

// Class of beam
class beam{
    public:
        double delay; // delay between adjustant bunches, [ns]
        double gap; // gap between bunch trains, [ns]
        
        bunch bnch = bunch(1);
        
        // Defualt initializer.
        beam(unsigned int N){
            delay = 1;
            gap = 20;
            bnch = bunch(N);
        };
        
        // Initialize with given normal distribution parameters
        beam(
            double dly, 
            double gp,
            std::vector<double> spc0, // vector that contain centers of distribution in spacial coordinates
            std::vector<double> mmt0, // vector that contain centers of distribution in momentum coordinates
            std::vector<double> spcSig, // vector that contain sigma of distribution in spacial coordinates
            std::vector<double> mmtSig, // vector that contain sigma of distribution in momentum coordinates
            unsigned int N){
            
            delay = dly;
            gap = gp;
            bnch = bunch(1,spc0,mmt0, spcSig, mmtSig,N);
        };
        

        // 
};
// Class of cavity
class cavity{
    public:
        double frq; // Frequency of the cavity, [MHz]
        double V0; // Voltage of the cavity, [MV]

        cavity(){
            frq = 704;
            V0 = 20;
        };
        void update_coord(particle& ptc){
            ptc.spc[0] += 0;
            ptc.spc[1] += 0;
            ptc.spc[2] += 0;
            ptc.mmt[0] += 0;
            ptc.mmt[1] += 0;
            double temp = 2*pi*frq*ptc.spc[2];
            ptc.mmt[2] += V0*sin(temp);
        };
};

// Class of ring
class ring{
    public:
        double alpx, alpy, betax,betay, gammax,gammay;
        double phix, phiy;
        double R;   // The radius of the ring.
        double cosphix;
        double sinphix;
        double cosphiy;
        double sinphiy;
        ring(){
            R = 610.1754;
            alpx = 0;
            alpy = 0;
            betax = 1;
            betay = 1;
            gammax = 0;
            gammay = 0;
            phix = 1.1;
            phiy = 1.2;
            cosphix = cos(phix);
            sinphix = sin(phix);
            cosphiy = cos(phiy);
            sinphiy = sin(phiy);
        };
        
        void update_coord(particle ptc){
            ptc.spc[0] = ptc.spc[0]*betax*(cosphix+alpx*sinphix)+betax*sinphix*ptc.mmt[0];
            ptc.spc[1] = ptc.spc[1]*betax*(cosphiy+alpx*sinphiy)+betax*sinphiy*ptc.mmt[1];
            ptc.spc[2] += 0;
            ptc.mmt[0] = -ptc.spc[0]*(1+alpx*alpx)/betax*sinphix+cosphix*ptc.mmt[0];
            ptc.mmt[1] = -ptc.spc[1]*(1+alpx*alpy)/betax*sinphiy+cosphiy*ptc.mmt[1];
            ptc.mmt[2] += 0;
        };
};

// Class of lattice
class lattice{
    public:
};
