#include <vector>
#include <iostream>
#include <random>

/*Global constants*/

double c = 3e8;
double pi = 3.1415926;


// Class of Bunch
class bunch {

	double M = 1;    /* mass number*/
	double N = 1;    /* charge number*/
	double q = N*1.60217662e-19;
	double m = M*9.10938356e-31;

public:
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
	double frq; // Frequency of the cavity, [MHz]
	double V0; // Voltage of the cavity, [MV]

	cavity() {
		frq = 704;
		V0 = 20;
	};
	void update_coord(bunch& bnch) {
#pragma omp parallel for        
		for (int i = 0; i<bnch.Np; ++i) {
			bnch.x[i] += 0;
			bnch.y[i] += 0;
			bnch.t[i] += 0;
			bnch.px[i] += 0;
			bnch.py[i] += 0;
			bnch.ga[i] += V0*sin(2 * pi*frq*bnch.t[i]);
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
