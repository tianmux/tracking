#include <vector>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <random>

/*Global constants*/

double c = 3e8;
double pi = 3.1415926;


// Class of Bunch
class bunch {
public:
	double M = 1e4;    /* mass number*/
	double N = 1e4;    /* charge number*/
	double qe = 1.60217662e-19;
	double q = N*qe;
	double me = 9.10938356e-31;
	double m = M*me;
    double gamma0 = 60; // everage gamma of the bunch. 
    double p0 = sqrt(gamma0*gamma0-1)*m/M*c;
    
	unsigned int Np;


	std::vector<double> x;
	std::vector<double> y;
	std::vector<double> t;
	std::vector<double> px;
	std::vector<double> py;
	std::vector<double> pz;


	std::vector<int> index;// holds the sorted indices of the particles, sort is based time coordinates of the particle. 

	// Initialize with given normal distribution parameters
	bunch(unsigned int disType,
		std::vector<double> spc0, // vector that contain centers of distribution in spacial coordinates
		std::vector<double> mmt0, // vector that contain centers of distribution in momentum coordinates
		std::vector<double> spcSig, // vector that contain sigma of distribution in spacial coordinates
		std::vector<double> mmtSig, // vector that contain sigma of distribution in momentum coordinates
		unsigned int Npar) {
		Np = Npar;
        
		x.resize(Np);
		y.resize(Np);
		t.resize(Np);
		px.resize(Np);
		py.resize(Np);
		pz.resize(Np);
        index.resize(Np);
		
		/* So far only the Gaussian distribution is supported */
		std::default_random_engine generator; // Normal distribution generator
		std::normal_distribution<double> dist0(spc0[0], spcSig[0]); // normal distribution in x
		std::normal_distribution<double> dist1(spc0[1], spcSig[1]); // normal distribution in y
		std::normal_distribution<double> dist2(spc0[2], spcSig[2]); // normal distribution in t
		std::normal_distribution<double> dist3(mmt0[0], mmtSig[0]); // normal distribution in px
		std::normal_distribution<double> dist4(mmt0[1], mmtSig[1]); // normal distribution in py
		std::normal_distribution<double> dist5(mmt0[2], mmtSig[2]); // normal distribution in gamma
	    // initialize the particles in one bunch in parallel.
#pragma omp parallel for
		for (int j = 0; j<Np; ++j) {
			x[j] = dist0(generator);
			y[j] = dist1(generator);
			t[j] = dist2(generator);
			px[j] = dist3(generator)*p0;
			py[j] = dist4(generator)*p0;
			pz[j] = (dist5(generator)+1)*p0;
			index[j] = j;
		}
	};

	// Initialize with defualt normal distribution.
	bunch(unsigned int N) {
		Np = N;
		x.resize(Np);
		y.resize(Np);
		t.resize(Np);
		px.resize(Np);
		py.resize(Np);
		pz.resize(Np);
        index.resize(Np);
		/* So far only the Gaussian distribution is supported */
		std::default_random_engine generator; // Normal distribution generator
		std::normal_distribution<double> dist0(0, 1e-3); // normal distribution in x
		std::normal_distribution<double> dist1(0, 1e-3); // normal distribution in y
		std::normal_distribution<double> dist2(0, 1e-10); // normal distribution in t
		std::normal_distribution<double> dist3(0, 1e-3); // normal distribution in px
		std::normal_distribution<double> dist4(0, 1e-3); // normal distribution in py
		std::normal_distribution<double> dist5(0, 1e-3); // normal distribution in gamma
		// initialize the particles in one bunch in parallel.
#pragma omp parallel for
		for (int j = 0; j<Np; ++j) {
			x[j] = dist0(generator);
			y[j] = dist1(generator);
			t[j] = dist2(generator);
			px[j] = dist3(generator)*p0;
			py[j] = dist4(generator)*p0;
			pz[j] = (dist5(generator)+1)*p0;
			index[j] = j;
		}
	};
// Initialize with defualt normal distribution.
	bunch() {
		Np = 1;
		x.resize(Np);
		y.resize(Np);
		t.resize(Np);
		px.resize(Np);
		py.resize(Np);
		pz.resize(Np);
        index.resize(Np);
		/* So far only the Gaussian distribution is supported */
		std::default_random_engine generator; // Normal distribution generator
		std::normal_distribution<double> dist0(0, 1e-3); // normal distribution in x
		std::normal_distribution<double> dist1(0, 1e-3); // normal distribution in y
		std::normal_distribution<double> dist2(0, 1e-10); // normal distribution in t
		std::normal_distribution<double> dist3(0, 1e-3); // normal distribution in px
		std::normal_distribution<double> dist4(0, 1e-3); // normal distribution in py
		std::normal_distribution<double> dist5(0, 1e-3); // normal distribution in gamma
		// initialize the particles in one bunch in parallel.
#pragma omp parallel for
		for (int j = 0; j<Np; ++j) {
			x[j] = dist0(generator);
			y[j] = dist1(generator);
			t[j] = dist2(generator);
			px[j] = dist3(generator)*p0;
			py[j] = dist4(generator)*p0;
			pz[j] = (dist5(generator)+1)*p0;
			index[j] = j;
		}
	};
	
	// Initialize the bunch with PARMELA output file.
	bunch(std::string path,double frq){
	    std::ifstream in;
	    std::string val;
	    std::vector<double> values;
	    in.open(path);
	    getline(in,val);
	    int i = 0;
	    while(in>>val){
	        values.push_back(stod(val));
	    }
	    
	    Np = values.size()/6;
	    x.resize(Np);
		y.resize(Np);
		t.resize(Np);
		px.resize(Np);
		py.resize(Np);
		pz.resize(Np);
        index.resize(Np);

	    for (int i = 0;i<Np;++i){
	        x[i] = values[6*i];// meter
	        px[i] = values[6*i+1];//rad, will be converted to real momentum later.
			y[i] = values[6*i+2];
			py[i] = values[6*i+3];
			t[i] = values[6*i+4]/180/(2*pi*frq);// raw is degree, convert to time
			pz[i] = values[6*i+5]*1e6*qe;// raw data is kinetic energy in MeV for single micro particle.
			index[i] = i;
	    }
	    std::cout<<"Number of particles per bunch: " <<Np<<std::endl;
	    std::cout<<"Kinetic energy per particle (MeV): "<<std::accumulate(pz.begin(),pz.end(),0.0)/Np/(1e6*qe)<<std::endl;
		for (int i = 0;i<Np;++i){
	        pz[i] = sqrt(pz[i]*pz[i]+2*pz[i]*me*c*c)/c; // convert the kinetic energy to momentum for micro particle.
	        px[i] = px[i]*pz[i]; // convert the transverse momentum to real momentum.
			py[i] = py[i]*pz[i];
	    }

	    p0 = std::accumulate(pz.begin(),pz.end(),0.0)/Np;
	    std::cout<<p0<<std::endl;
		gamma0 = sqrt(p0*p0/(me*c*c*me)+1);
		std::cout<<gamma0<<std::endl;
	};
	
    // get the sorted index based on the time coordinates of the particles.
    void sort(){
        std::sort(index.begin(),index.end(),[&](const double& a,const double& b){return (t[a]<t[b]);}
        );
    };
    
	// dump the whole bunch to a file. file name is specified by 'path'.
	void dump_to_file(std::string path) {
		std::filebuf fb;
		fb.open(path, std::ios::out);
		std::ostream os(&fb);
		for (int i = 0; i<Np; ++i) {
			os << x[i] << "," << y[i] << "," << t[i] << "," << px[i] << "," << py[i] << "," << pz[i] << '\n';
		}
		fb.close();
	};
};

// Class of beam
class beam {
public:
	double delay; // delay between adjustant bunches, [ns]
	double gap; // gap between bunch trains, [ns]
    unsigned int N_bnch; // total number of bunches.
	std::vector<bunch> bnches;

	// Defualt initializer.
	beam(unsigned int N,unsigned int N_bnches) {
        bnches.resize(N_bnches);
		delay = 1;
		gap = 20;
        for (int i = 0;i<N_bnches;++i){
		    bnches[i] = bunch(N);
		}
	};
    // Initialize with PARMELA output file
    beam(std::string path,double frq,unsigned int N_bnches) {
        bnches.resize(N_bnches);
		delay = 1;
		gap = 20;
        for (int i = 0;i<N_bnches;++i){
		    bnches[i] = bunch(path,frq);
		}
	};
	// Initialize with given normal distribution parameters
	beam(
		double dly,
		double gp,
		std::vector<double> spc0, // vector that contain centers of distribution in spacial coordinates
		std::vector<double> mmt0, // vector that contain centers of distribution in momentum coordinates
		std::vector<double> spcSig, // vector that contain sigma of distribution in spacial coordinates
		std::vector<double> mmtSig, // vector that contain sigma of distribution in momentum coordinates
		unsigned int N, unsigned int N_bnches) // particles per bunch and bunches per beam.
		{
        bnches.resize(N_bnches);
		delay = dly;
		gap = gp;
		for (int i = 0;i<N_bnches;++i){
		    bnches[i] = bunch(1, spc0, mmt0, spcSig, mmtSig, N);
		}
	};
};


// Class of cavity
class cavity {
public:
	std::vector<double> frq; // Frequency of the cavity, [MHz]
	std::vector<double> V0xR,V0yR,V0zR; // Real parts of the Voltage of the cavity, [MV]
	std::vector<double> V0xI,V0yI,V0zI; // Imaginary parts of the Voltage of the cavity, [MV]
	
	
    std::vector<std::vector<double>> VxTauR;
	std::vector<std::vector<double>> VyTauR;
	std::vector<std::vector<double>> VzTauR;
	std::vector<std::vector<double>> VxTauI;
	std::vector<std::vector<double>> VyTauI;
	std::vector<std::vector<double>> VzTauI;// vectors to hold the cumulated wake info at each time point (corresponding to each particle), one vector for each mode.

    std::vector<double> kx; 
    std::vector<double> ky;
    std::vector<double> kz;// loss factors for each mode, use this to calculate wake of each mode.
    std::vector<double> tau_invert; // one over tau, save some calculation for each particle.
    std::vector<double> phi; // phase of each mode.
    int N_mod=1;// number of modes to consider
    int bnch_count = 0; // The number of bunches that has already passed the cavity, used to calculate wake field.
    
	cavity() {
		frq.resize(N_mod,0.0);
		phi.resize(N_mod,0.0);
		V0xR.resize(N_mod,0.0);
		V0yR.resize(N_mod,0.0);
		V0zR.resize(N_mod,0.0);
		V0xI.resize(N_mod,0.0);
		V0yI.resize(N_mod,0.0);
		V0zI.resize(N_mod,0.0);
		VxTauR.resize(N_mod);
		VyTauR.resize(N_mod);
		VzTauR.resize(N_mod);
		VxTauI.resize(N_mod);
		VyTauI.resize(N_mod);
		VzTauI.resize(N_mod);
		kx.resize(N_mod,0.0);
		ky.resize(N_mod,0.0);
		kz.resize(N_mod,0.0);
		tau_invert.resize(N_mod,0.0);
	};
	
	// Most naive way of calculating the wake, basically sum all the wake from each macro particles for each mode, 
	// Equavlent to convolute the bunch with wake function, slow as hell...
	
	void wake_Naive(beam& bm, bunch& bnch){

	    for (int i = 0;i<N_mod;++i){
	        VxTauR[i].resize(bnch.Np+1,0); // extra one element at the end, use to store the wake in this mode when bunch leaves the cavity.
	        VyTauR[i].resize(bnch.Np+1,0);
	        VzTauR[i].resize(bnch.Np+1,0);
	        VxTauI[i].resize(bnch.Np+1,0); 
	        VyTauI[i].resize(bnch.Np+1,0);
	        VzTauI[i].resize(bnch.Np+1,0);
	    }
        int j = 0;
#pragma omp parallel for private(j)
	    for(int i = 0;i<N_mod;++i){ // iterate over number of modes.
	        double Vbx = kx[i]*bnch.N;// always negative(?) imaginary.
	        double Vby = ky[i]*bnch.N;// always negative(?) imaginary.
	        double Vbz = kz[i]*bnch.N;// always negative real.
	        double decay = exp(-bm.delay*tau_invert[i]);// decay of wake from last bunch
	        double dphi = bm.delay*frq[i]*2.0*pi; // phase shift of the wake from last bunch.
	        // rotated wake from last bunch:
	        VxTauR[i][0] = VxTauR[i][bnch.Np]*decay*cos(dphi)-VxTauI[i][bnch.Np]*decay*sin(dphi);
	        VyTauR[i][0] = VyTauR[i][bnch.Np]*decay*cos(dphi)-VyTauI[i][bnch.Np]*decay*sin(dphi);
	        VzTauR[i][0] = VzTauR[i][bnch.Np]*decay*cos(dphi)-VzTauI[i][bnch.Np]*decay*sin(dphi);
	        VxTauI[i][0] = VxTauR[i][bnch.Np]*decay*sin(dphi)+VxTauI[i][bnch.Np]*decay*cos(dphi);
	        VyTauI[i][0] = VyTauR[i][bnch.Np]*decay*sin(dphi)+VyTauI[i][bnch.Np]*decay*cos(dphi);
	        VzTauI[i][0] = VzTauR[i][bnch.Np]*decay*sin(dphi)+VzTauI[i][bnch.Np]*decay*cos(dphi);
	        
	        
	        for (j = 1;j<bnch.Np;++j){ // iterate over every particles. 
	                                   // ith particle sees the wake of all
	                                   // previous particles with a phase shift
	                                   // w(t(i)-t(i-1))
	                                   // ignore the decay of the wake between
	                                   // two macro particles. 
	            double dt = bnch.t[bnch.index[j]]-bnch.t[bnch.index[j-1]];
	            double cosin = cos(2*pi*frq[i]*dt);
	            double sine = sin(2*pi*frq[i]*dt);
	            VxTauR[i][j] = VxTauR[i][j-1]*cosin-(VxTauI[i][j-1]+Vbx)*sine;
	            VyTauR[i][j] = VyTauR[i][j-1]*cosin-(VyTauI[i][j-1]+Vby)*sine;
	            VzTauR[i][j] = (VzTauR[i][j-1]+Vbz)*cosin-VzTauI[i][j-1]*sine;
	            VxTauI[i][j] = VxTauR[i][j-1]*sine+(VxTauI[i][j-1]+Vbx)*cosin;
	            VyTauI[i][j] = VyTauR[i][j-1]*sine+(VyTauI[i][j-1]+Vby)*cosin;
	            VzTauI[i][j] = (VzTauR[i][j-1]+Vbz)*sine+VzTauI[i][j-1]*cosin;
	        }
	        
	        VxTauR[i][bnch.Np] = VxTauR[i][j-1];
	        VyTauR[i][bnch.Np] = VyTauR[i][j-1];
	        VzTauR[i][bnch.Np] = VzTauR[i][j-1]+Vbz;
	        VxTauI[i][bnch.Np] = VxTauI[i][j-1]+Vbx;
	        VyTauI[i][bnch.Np] = VyTauI[i][j-1]+Vby;
	        VzTauI[i][bnch.Np] = VzTauI[i][j-1];
	    }
	};
	
	void update_coord(bunch& bnch) {
	    double qoc = bnch.qe/c;
	    //wake_Naive(bnch);
#pragma omp parallel	    
	    for (int i = 0;i<N_mod;++i){ // iterate over all modes
	        double fi=frq[i];
	        double phiN = phi[i]/180*pi;
	        double taui=tau_invert[i];
	        double Vbz = kz[i]*bnch.N;// always negative real.
#pragma omp for
		    for (int j = 0; j<bnch.Np; ++j) { // iterate over all particles

			    int tempID = bnch.index[j];
			    double cosphi = cos(2 * pi*fi*bnch.t[tempID]+phiN);
			    double sinphi = sin(2 * pi*fi*bnch.t[tempID]+phiN);
			    double expi = exp(-(bnch.t[tempID])*taui);
			    // The voltage each particle sees is the real part of the complex voltage in cavity at that time, 
			    // It should be the combination of the cavity voltage (including the previous wake_voltage), and half of the self field(z). 
			    bnch.px[tempID] += qoc*(V0xR[i]*cosphi-V0xI[i]*sinphi+VxTauR[i][j]); // kick in x direction. cannot see self field.
			    bnch.py[tempID] += qoc*(V0yR[i]*cosphi-V0yI[i]*sinphi+VyTauR[i][j]); // kick in y direction.
			    bnch.pz[tempID] += qoc*(V0zR[i]*cosphi-V0zI[i]*sinphi+VzTauR[i][j]+Vbz*0.5); // kick in z direction.
		    }
		}
		// Update the bunch momentum and energy.
		bnch.p0 = std::accumulate(bnch.pz.begin(),bnch.pz.end(),0.0)/bnch.pz.size();
		bnch.gamma0 = sqrt(bnch.p0*bnch.p0/(bnch.me*c)+1);
	};
};

// Class of ring
class ring {
public:
	double alpx, alpy, betax, betay, gammax, gammay; // ring lattice parameters.
	double phix, phiy;
	double R;   // The radius of the ring.
	double GAMTSQ; // Transition gamma.
	double cosphix;
	double sinphix;
	double cosphiy;
	double sinphiy;
	double *TM;//Transfer matrix.
	ring() {
		R = 610.1754; // RHIC ring.
		GAMTSQ = 691.69;
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
		TM = new double[6 * 6];
		for (unsigned int i = 0; i < 36; ++i) {
			TM[i] = 0;
		}
		TM[0] = betax*(cosphix + alpx*sinphix);
		TM[1] = betax*sinphix;
		TM[6] = -(1 + alpx*alpx) / betax*sinphix;
		TM[7] = cosphix;
		TM[14] = betay*(cosphiy + alpy*sinphiy);
		TM[15] = betay*sinphiy;
		TM[20] = -(1 + alpy*alpy) / betay*sinphiy;
		TM[21] = cosphiy;
	};
	void updt_coord_blas(double* bnch) {

	};
	void update_coord(bunch& bnch) {
	    double alpha0 = 1/GAMTSQ;
	    double GM0SQ_inv = 1/(bnch.gamma0*bnch.gamma0);
	    double p0_inv = 1/bnch.p0;
	    double ita = alpha0-GM0SQ_inv;
	    double f0_inv = (c*sqrt(1-GM0SQ_inv))/(2*pi*R);
#pragma omp parallel for        
        // zeroth order transport
		for (int i = 0; i<bnch.Np; ++i) {
			bnch.x[i] = bnch.x[i] * betax*(cosphix + alpx*sinphix) + betax*sinphix*bnch.px[i];
			bnch.y[i] = bnch.y[i] * betay*(cosphiy + alpy*sinphiy) + betay*sinphiy*bnch.py[i];
			bnch.t[i] += f0_inv/(1-ita*(bnch.pz[i]-bnch.p0)/bnch.p0);
			bnch.px[i] = -bnch.x[i] * (1 + alpx*alpx) / betax*sinphix + cosphix*bnch.px[i];
			bnch.py[i] = -bnch.y[i] * (1 + alpy*alpy) / betay*sinphiy + cosphiy*bnch.py[i];
			bnch.pz[i] += 0;

		}
	};
};

// Class of drift space
class drift_space{
public:
    double L;// length of the drift space.
    drift_space(double dft_L = 1.73){
        L = dft_L;// Default length is 1 m.
    }
    void update_coord(bunch& bnch){
        double p0_invert = 1/bnch.p0;

#pragma omp paprallel for 
        for (int i = 0;i<bnch.Np;++i){
            bnch.x[i] += L*bnch.px[i]*p0_invert;
            bnch.y[i] += L*bnch.py[i]*p0_invert;
        }
    }
};
// Class of lattice
class lattice {
public:
};

// Class of input parameters
class inputPara{
public:
    double* frqs; // frequencies of each mode in a cavity.
    double* V0s; // Existing voltage of each mode in a cavity at the beginning.
    double* RoQs; // R over Qs of each mode in a cavity.
    double* Qs; // Loaded Qs of each mode in a cavity.
    double* ks; // loss factors of each mode in a cabity.
    double R=610.1754; // radius of the ring.
    double gammaT; // Transision energy of the ring.
    
};

