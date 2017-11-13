#include <vector>
#include <iostream>
#include <random>

/*Global constants*/

double c = 3e8;
double pi = 3.1415926;


// Class of Bunch
class bunch {


public:
	double M = 1;    /* mass number*/
	double N = 1;    /* charge number*/
	double q = N*1.60217662e-19;
	double m = M*9.10938356e-31;

	unsigned int Np;
	double gamma; // reference gamma of the bunch.
	double *x, *y, *t, *px, *py, *ga;
	double *sixD;
	// Initialize with given normal distribution parameters
	bunch(unsigned int disType,
		std::vector<double> spc0, // vector that contain centers of distribution in spacial coordinates
		std::vector<double> mmt0, // vector that contain centers of distribution in momentum coordinates
		std::vector<double> spcSig, // vector that contain sigma of distribution in spacial coordinates
		std::vector<double> mmtSig, // vector that contain sigma of distribution in momentum coordinates
		unsigned int N) {
		Np = N;
		x = new double[Np];
		y = new double[Np];
		t = new double[Np];
		px = new double[Np];
		py = new double[Np];
		ga = new double[Np];
		sixD = new double[Np * 6];
		/* So far only the Gaussian distribution is supported */
		std::default_random_engine generator; // Normal distribution generator
		std::normal_distribution<double> dist0(spc0[0], spcSig[0]); // normal distribution in x
		std::normal_distribution<double> dist1(spc0[1], spcSig[1]); // normal distribution in y
		std::normal_distribution<double> dist2(spc0[2], spcSig[2]); // normal distribution in z
		std::normal_distribution<double> dist3(mmt0[0], mmtSig[0]); // normal distribution in px
		std::normal_distribution<double> dist4(mmt0[1], mmtSig[1]); // normal distribution in py
		std::normal_distribution<double> dist5(mmt0[2], mmtSig[2]); // normal distribution in pz
																	// initialize the particles in one bunch in parallel.
#pragma omp parallel for
		for (int j = 0; j<Np; ++j) {
			x[j] = dist0(generator);
			y[j] = dist1(generator);
			t[j] = dist2(generator);
			px[j] = dist3(generator);
			py[j] = dist4(generator);
			ga[j] = dist5(generator);
			sixD[j] = x[j];
			sixD[j+1] = y[j];
			sixD[j+2] = t[j];
			sixD[j+3] = px[j];
			sixD[j+4] = py[j];
			sixD[j+5] = ga[j];
		}
	};

	// Initialize with defualt normal distribution.
	bunch(unsigned int N) {
		Np = N;
		x = new double[Np];
		y = new double[Np];
		t = new double[Np];
		px = new double[Np];
		py = new double[Np];
		ga = new double[Np];
		sixD = new double[Np * 6];

		/* So far only the Gaussian distribution is supported */
		std::default_random_engine generator; // Normal distribution generator
		std::normal_distribution<double> dist0(0, 1e-3); // normal distribution in x
		std::normal_distribution<double> dist1(0, 1e-3); // normal distribution in y
		std::normal_distribution<double> dist2(0, 1e-3); // normal distribution in z
		std::normal_distribution<double> dist3(0, 1e-3); // normal distribution in px
		std::normal_distribution<double> dist4(0, 1e-3); // normal distribution in py
		std::normal_distribution<double> dist5(0, 1e-3); // normal distribution in pz
														 // initialize the particles in one bunch in parallel.
#pragma omp parallel for
		for (int j = 0; j<Np; ++j) {
			x[j] = dist0(generator);
			y[j] = dist1(generator);
			t[j] = dist2(generator);
			px[j] = dist3(generator);
			py[j] = dist4(generator);
			ga[j] = dist5(generator);
			sixD[j] = x[j];
			sixD[j + 1] = y[j];
			sixD[j + 2] = t[j];
			sixD[j + 3] = px[j];
			sixD[j + 4] = py[j];
			sixD[j + 5] = ga[j];
		}
	};

	// print the 6D coordiantes of a couple of particle, mainly for debug
	void display() {

	};
	// dump the whole bunch to a file. file name is specified by 'path'.
	void dump_to_file(std::string path) {
		std::filebuf fb;
		fb.open(path, std::ios::out);
		std::ostream os(&fb);
		for (int i = 0; i<Np; ++i) {
			os << x[i] << "," << y[i] << "," << t[i] << ";" << px[i] << "," << py[i] << "," << ga[i] << '\n';
		}
	};

};

// Class of beam
class beam {
public:
	double delay; // delay between adjustant bunches, [ns]
	double gap; // gap between bunch trains, [ns]

	bunch bnch = bunch(1);

	// Defualt initializer.
	beam(unsigned int N) {
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
		unsigned int N) {

		delay = dly;
		gap = gp;
		bnch = bunch(1, spc0, mmt0, spcSig, mmtSig, N);
	};


	// 
};


// Class of cavity
class cavity {
public:
	double *frq; // Frequency of the cavity, [MHz]
	double *V0r; // Real parts of the Voltage of the cavity, [MV],x,y,z direction in the fashion of [Vx1,Vy1,Vz1,Vx2,Vy2,Vz2,...,VxN,VyN,VzN]
	double *V0i; // Imaginary parts of the Voltage of the cavity,
    double *k; // loss factors for each mode, use this to calculate wake of each mode.
    double *tau_invert; // decay factor for each mode 
    double *phi; // phase of each mode.
    int N=2;// number of modes to consider
    int bnch_cout = 0; // The number of bunches that has already passed the cavity, used to calculate wake field.
    
	cavity() {
		frq = new double[N];
		phi = new double[N];
		V0r = new double[N*3];
		V0i = new double[N*3];
		k = new double[N*3];
		tau_invert = new double[N];
	};
	// Most naive way of calculating the wake, basically sum all the wake from each macro particles for each mode, 
	// Equavlent to convolute the bunch with wake function, slow as hell...
	void wake_Naive(bunch& bnch){
	    for(int i = 0;i<N;++i){ // iterate over number of modes.
	        for (int j = 0;j<bnch.Np;++j){ // iterate over every particles.
	            V0r[i+2] +=k[i]*bnch.N*cos(2*pi*frq[i]*bnch.t[i]);
	        }
	    }
	};
	
	void update_coord(bunch& bnch) {
	    for (int j = 0;j<N;++j){ // iterate over all modes
	        double fj=frq[j];
	        double tau=tau_invert[j];
#pragma omp parallel for        
		    for (int i = 0; i<bnch.Np; ++i) { // iterate over all particles
			    //bnch.x[i] += 0;
			    //bnch.y[i] += 0;
			    //bnch.t[i] += 0;
			    double cosphi = cos(2 * pi*fj*bnch.t[i]+phi[j]);
			    bnch.px[i] += V0r[j]* cosphi*exp(-(bnch.t[i])*tau); // kick in x direction.
			    bnch.py[i] += V0r[j+1]* cosphi*exp(-(bnch.t[i])*tau); // kick in y direction.
			    bnch.ga[i] += V0r[j+2]* cosphi*exp(-(bnch.t[i])*tau); //kick in z direction.
		    }
		}
	};
};

// Class of ring
class ring {
public:
	double alpx, alpy, betax, betay, gammax, gammay; // ring lattice parameters.
	double phix, phiy;
	double R;   // The radius of the ring.
	double cosphix;
	double sinphix;
	double cosphiy;
	double sinphiy;
	double *TM;//Transfer matrix.
	ring() {
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
#pragma omp parallel for        
		for (int i = 0; i<bnch.Np; ++i) {
			bnch.x[i] = bnch.x[i] * betax*(cosphix + alpx*sinphix) + betax*sinphix*bnch.px[i];
			bnch.y[i] = bnch.y[i] * betay*(cosphiy + alpy*sinphiy) + betay*sinphiy*bnch.py[i];
			bnch.t[i] += 0;
			bnch.px[i] = -bnch.x[i] * (1 + alpx*alpx) / betax*sinphix + cosphix*bnch.px[i];
			bnch.py[i] = -bnch.y[i] * (1 + alpy*alpy) / betay*sinphiy + cosphiy*bnch.py[i];
			bnch.ga[i] += 0;

		}
	};
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

