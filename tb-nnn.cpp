# include <cstdlib>
# include <iostream>
# include <iomanip>
# include <cmath>
# include <vector>
# include <complex>
# include <eigen3/Eigen/Eigenvalues>
# include <fstream>

using namespace std;
using namespace Eigen;


///////////////// Edit residual function so that it takes in tb parameters and 
///////////////// spits out residual without fourier transform yet. Then look
///////////////// at the nelder mead algorithm and change accordingly

vector<double> dot_product(vector<int> vec1, vector<int> vec2, int typeint, vector<vector<double> > neighbors, vector<vector<double> > next_neighbors){
	vector<double> p;
	int i;
	if(typeint == 1){
		for(i = 0; i<neighbors.size(); i++){
			p.push_back((vec1[0]*neighbors[i][0]+vec1[1]*neighbors[i][1]+vec1[2]*neighbors[i][2])
		        /(pow(neighbors[i][0]*neighbors[i][0]+neighbors[i][1]*neighbors[i][1]
		        +neighbors[i][2]*neighbors[i][2],0.5)));
		}
	}
	if(typeint == 2){
		for(i = 0; i<neighbors.size(); i++){
			p.push_back(((vec1[0]*neighbors[i][0]+vec1[1]*neighbors[i][1]+vec1[2]*neighbors[i][2])
	            /(pow(neighbors[i][0]*neighbors[i][0]+neighbors[i][1]*neighbors[i][1]
	            +neighbors[i][2]*neighbors[i][2],0.5)))*((vec2[0]*neighbors[i][0]+vec2[1]*neighbors[i][1]
			    +vec2[2]*neighbors[i][2])/(pow(neighbors[i][0]*neighbors[i][0]+neighbors[i][1]*neighbors[i][1]
			    +neighbors[i][2]*neighbors[i][2],0.5))));
		}
	}
	if(typeint == 3){
		for(i = 0; i<next_neighbors.size(); i++){
			p.push_back((vec1[0]*next_neighbors[i][0]+vec1[1]*next_neighbors[i][1]+vec1[2]*next_neighbors[i][2])
				/(pow(next_neighbors[i][0]*next_neighbors[i][0]+next_neighbors[i][1]*next_neighbors[i][1]
			    +next_neighbors[i][2]*next_neighbors[i][2],0.5)));
		}
	}
	if(typeint == 4){
		for(i = 0; i<next_neighbors.size(); i++){
			p.push_back(((vec1[0]*next_neighbors[i][0]+vec1[1]*next_neighbors[i][1]+vec1[2]*next_neighbors[i][2])
	            /(pow(next_neighbors[i][0]*next_neighbors[i][0]+next_neighbors[i][1]*next_neighbors[i][1]
	            +next_neighbors[i][2]*next_neighbors[i][2],0.5)))*((vec2[0]*next_neighbors[i][0]+vec2[1]*next_neighbors[i][1]
			    +vec2[2]*next_neighbors[i][2])/(pow(next_neighbors[i][0]*next_neighbors[i][0]+next_neighbors[i][1]*next_neighbors[i][1]
			    +next_neighbors[i][2]*next_neighbors[i][2],0.5))));
		}
	}			
	return p;
}

vector<complex<double> > phase(vector<double> k, vector<vector<double> > neighbors, vector<vector<double> > next_neighbors){
	int i;
	vector<int> x{1,0,0};
	vector<int> y{0,1,0};
	vector<int> z{0,0,1}; 
	vector<complex<double> > nnp;
	vector<vector<double> > nnd;
	vector<complex<double> > nnnp;
	vector<vector<double> >nnnd;
	vector<vector<double> > empty{};
	vector<int> zero_vec;
	complex<double> j(0.0,1.0);

	for(i = 0; i<neighbors.size();i++){
		nnp.push_back(exp(j*(k[0]*neighbors[i][0]+k[1]*neighbors[i][1]+k[2]*neighbors[i][2])));
	}
	nnd.push_back(dot_product(x,zero_vec,1,neighbors,empty));
	nnd.push_back(dot_product(y,zero_vec,1,neighbors,empty));
	nnd.push_back(dot_product(z,zero_vec,1,neighbors,empty));
	nnd.push_back(dot_product(x,x,2,neighbors,empty));
	nnd.push_back(dot_product(y,y,2,neighbors,empty));
	nnd.push_back(dot_product(z,z,2,neighbors,empty));
	nnd.push_back(dot_product(x,y,2,neighbors,empty));
	nnd.push_back(dot_product(x,z,2,neighbors,empty));
	nnd.push_back(dot_product(y,z,2,neighbors,empty));

	for(i = 0; i<next_neighbors.size();i++){
		nnnp.push_back(exp(j*(k[0]*next_neighbors[i][0]+k[1]*next_neighbors[i][1]+k[2]*next_neighbors[i][2])));
	}
	nnnd.push_back(dot_product(x,zero_vec,3,empty,next_neighbors));
	nnnd.push_back(dot_product(y,zero_vec,3,empty,next_neighbors));
	nnnd.push_back(dot_product(z,zero_vec,3,empty,next_neighbors));
	nnnd.push_back(dot_product(x,x,4,empty,next_neighbors));
	nnnd.push_back(dot_product(y,y,4,empty,next_neighbors));
	nnnd.push_back(dot_product(z,z,4,empty,next_neighbors));
	nnnd.push_back(dot_product(x,y,4,empty,next_neighbors));
	nnnd.push_back(dot_product(x,z,4,empty,next_neighbors));
	nnnd.push_back(dot_product(y,z,4,empty,next_neighbors));

	vector<complex<double> > factors;
	complex<double> sumfac(0.0);
	for(int l = 0; l < nnd.size(); l++){
		sumfac = 0.0;
		for(int m = 0; m < nnp.size(); m++){
			sumfac += nnp[m]*nnd[l][m];
		}
		factors.push_back(sumfac);
	}
	for(int l = 0; l < nnd.size(); l++){
		sumfac = 0.0;
		for(int m = 0; m < nnp.size(); m++){
			sumfac += nnnp[m]*nnnd[l][m];
		}
		factors.push_back(sumfac);
	}
	return factors;
}

struct eigen{
	vector<complex<double> > evalues;
	vector<vector<complex<double> > > evectors;
};

eigen band_eigenvalues_and_eigenvectors(complex<double> l1,   complex<double> m1,   complex<double> n1, 
					 complex<double> l1l1, complex<double> m1m1, complex<double> n1n1, 
					 complex<double> l1m1, complex<double> l1n1, complex<double> m1n1, 
					 complex<double> l2,   complex<double> m2,   complex<double> n2, 
					 complex<double> l2l2, complex<double> m2m2, complex<double> n2n2, 
					 complex<double> l2m2, complex<double> l2n2, complex<double> m2n2,
					 double es, double ep, double vssS, double vspS, double vppS, double vppP, 
					 double vssS2, double vspS2, double vppS2, double vppP2){
	complex<double> one(1.0,0.0);
	complex<double> l1c = conj(l1);
	complex<double> m1c = conj(m1);
	complex<double> n1c = conj(n1);
	complex<double> l1l1c = conj(l1l1);
	complex<double> m1m1c = conj(m1m1);
	complex<double> n1n1c = conj(n1n1);
	complex<double> l1m1c = conj(l1m1);
	complex<double> l1n1c = conj(l1n1);
	complex<double> m1n1c = conj(m1n1);
	complex<double> l2c = conj(l2);
	complex<double> m2c = conj(m2);
	complex<double> n2c = conj(n2);
	complex<double> l2l2c = conj(l2l2);
	complex<double> m2m2c = conj(m2m2);
	complex<double> n2n2c = conj(n2n2);
	complex<double> l2m2c = conj(l2m2);
	complex<double> l2n2c = conj(l2n2);
	complex<double> m2n2c = conj(m2n2);

	Matrix<complex<double>, 12, 12> hamiltonian;
	hamiltonian <<          es, l1l1*vssS, l2l2c*vssS2,                           0,                           0,                           0,                     l1*vspS,                     m1*vspS,                     n1*vspS,                    -l2c*vspS2,                    -m2c*vspS2,                    -n2c*vspS2,
                    l1l1c*vssS,         es,  l1l1*vssS,                   -l1c*vspS,                   -m1c*vspS,                   -n1c*vspS,                           0,                           0,                           0,                       l1*vspS,                       m1*vspS,                       n1*vspS,
        			l2l2*vssS2, l1l1c*vssS,         es,                    l2*vspS2,                    m2*vspS2,                    n2*vspS2,                   -l1c*vspS,                   -m1c*vspS,                   -n1c*vspS,                             0,                             0,                             0,
        			         0,    l1*vspS, -l2c*vspS2,                          ep,                           0,                           0,   l1l1*vppS+(one-l1l1)*vppP,            l1m1*(vppS-vppP),            l1n1*(vppS-vppP), l2l2c*vppS2+(one-l2l2c)*vppP2,           l2m2c*(vppS2-vppP2),           l2n2c*(vppS2-vppP2),
        			         0,    m1*vspS, -m2c*vspS2,                           0,                          ep,                           0,            l1m1*(vppS-vppP),   m1m1*vppS+(one-m1m1)*vppP,            m1n1*(vppS-vppP),           l2m2c*(vppS2-vppP2), m2m2c*vppS2+(one-m2m2c)*vppP2,           m2n2c*(vppS2-vppP2),
        			         0,    n1*vspS, -n2c*vspS2,                           0,                           0,                          ep,            l1n1*(vppS-vppP),            m1n1*(vppS-vppP),   n1n1*vppS+(one-n1n1)*vppP,           l2n2c*(vppS2-vppP2),           m2n2c*(vppS2-vppP2), n2n2c*vppS2+(one-n2n2c)*vppP2,
        			 -l1c*vspS,          0,    l1*vspS, l1l1c*vppS+(one-l1l1c)*vppP,           l1m1c*(vppS-vppP),           l1n1c*(vppS-vppP),                          ep,                           0,                           0,     l1l1*vppS+(one-l1l1)*vppP,              l1m1*(vppS-vppP),              l1n1*(vppS-vppP),
        			 -m1c*vspS,          0,    m1*vspS,           l1m1c*(vppS-vppP), m1m1c*vppS+(one-m1m1c)*vppP,           m1n1c*(vppS-vppP),                           0,                          ep,                           0,              l1m1*(vppS-vppP),     m1m1*vppS+(one-m1m1)*vppP,              m1n1*(vppS-vppP),
        			 -n1c*vspS,          0,    n1*vspS,           l1n1c*(vppS-vppP),           m1n1c*(vppS-vppP), n1n1c*vppS+(one-n1n1c)*vppP,                           0,                           0,                          ep,              l1n1*(vppS-vppP),              m1n1*(vppS-vppP),     n1n1*vppS+(one-n1n1)*vppP,
        			  l2*vspS2,  -l1c*vspS,          0, l2l2*vppS2+(one-l2l2)*vppP2,          l2m2*(vppS2-vppP2),          l2n2*(vppS2-vppP2), l1l1c*vppS+(one-l1l1c)*vppP,           l1m1c*(vppS-vppP),           l1n1c*(vppS-vppP),                            ep,                             0,                             0,
        			  m2*vspS2,  -m1c*vspS,          0,            l2m2*(vppS-vppP),  m2m2*vppS+(one-m2m2)*vppP2,          m2n2*(vppS2-vppP2),           l1m1c*(vppS-vppP), m1m1c*vppS+(one-m1m1c)*vppP,           m1n1c*(vppS-vppP),                             0,                            ep,                             0,
        			  m2*vspS2,  -n1c*vspS,          0,            l2n2*(vppS-vppP),          m2n2*(vppS2-vppP2), n2n2*vppS2+(one-n2n2)*vppP2,           l1n1c*(vppS-vppP),           m1n1c*(vppS-vppP), n1n1c*vppS+(one-n1n1c)*vppP,                             0,                             0,                            ep;
    SelfAdjointEigenSolver<Matrix<complex<double>, 12, 12> > s(hamiltonian);
   
    vector<complex<double> > evals;
    vector<vector<complex<double> > > evecs;
   
    for(int n = 0; n < 12; n++){
    	vector<complex<double> > evec;
    	evals.push_back(s.eigenvalues()(n));
    	for(int p = 0; p < 12; p++){
    		evec.push_back(s.eigenvectors()(n,p));
    	}
    	evecs.push_back(evec);
    }

    eigen e = {};
    e.evalues = evals;
    e.evectors = evecs;

	return e;}

struct band_struct{
	vector<vector<complex<double> > > bandvalues;
	vector<vector<vector<complex<double> > > > bandvectors;
};

band_struct band_structure(vector<double> params, vector<vector<double> > neighbors, 
	vector<vector<double> > next_neighbors, vector<vector<double> > path){
	vector<vector<complex<double> > > bands;
	vector<vector<vector<complex<double> > > > bandvecs;
	vector<complex<double> > g;
	eigen e;

	for(int k = 0; k < path.size(); k++){
		g = phase(path[k], neighbors, next_neighbors);
		e = band_eigenvalues_and_eigenvectors(g[0], g[1],  g[2],  g[3], g[4], g[5],  g[6], g[7], 
			                                  g[8], g[9], g[10], g[11], g[12], g[13], g[14], g[15],
			                                  g[16], g[17], params[0], params[1], params[2], params[3],
			                                  params[4], params[5], params[6],params[7],params[8],params[9]);
		bands.push_back(e.evalues);
		bandvecs.push_back(e.evectors);
	}
	band_struct b = {};
	b.bandvalues = bands;
	b.bandvectors = bandvecs;

	return b;}

double OVERLAP(vector<complex<double>> dft_band_eigenvector_complex, vector<complex<double>> tb_band_eigenvector){
	double overlap = 0;
	complex<double> normal = 0;
	for(int element = 0; element < dft_band_eigenvector_complex.size(); element++){
		normal += conj(dft_band_eigenvector_complex[element])*tb_band_eigenvector[element];
	}
	overlap = double(norm(normal));
	return overlap;
}

double residual(band_struct band_structure(vector<double> params, vector<vector<double> > neighbors, 
	vector<vector<double> > next_neighbors, vector<vector<double> > dft_path), vector<double> params, 
	vector<vector<double> > neighbors, vector<vector<double> > next_neighbors,
	vector<vector<double> > dft_path, vector<vector<double> > dft_bands,
	vector<vector<vector<complex<double> > > > dft_band_eigenvectors_complex){
	
	band_struct BAND_STRUCTURE = band_structure(params, neighbors, next_neighbors, dft_path);
	vector<vector<complex<double> > > BAND_VALUES = BAND_STRUCTURE.bandvalues;
	vector<vector<vector<complex<double> > > > tb_band_eigenvectors = BAND_STRUCTURE.bandvectors;

	vector<double> band;
	vector<vector<double> > bands(BAND_VALUES.size(),vector<double> (BAND_VALUES[0].size()));
	vector<vector<double> > tb_bands;
	double normalization = 0;
	for(int i = 0; i < BAND_VALUES.size(); i++){
		for(int j = 0; j < BAND_VALUES[0].size(); j++){
			bands[i][j] = real(BAND_VALUES[i][j]);
			normalization = OVERLAP(tb_band_eigenvectors[i][j],tb_band_eigenvectors[i][j]);
			for(int element = 0; element < tb_band_eigenvectors[i][j].size(); element++){
				tb_band_eigenvectors[i][j][element] /= normalization;
			}
		}
	}

	for(int i = 0; i < dft_bands.size(); i++){
		for(int j = 0; j < dft_bands[0].size(); j++){
			normalization = OVERLAP(dft_band_eigenvectors_complex[i][j],dft_band_eigenvectors_complex[i][j]);
			for(int element = 0; element < dft_band_eigenvectors_complex[i][j].size(); element++){
				dft_band_eigenvectors_complex[i][j][element] /= normalization;
			}
		}
	}

	for(int i = 0; i < bands.size(); i++){
		band = bands[i];
		// sort(band.begin(),band.end());
		tb_bands.push_back(band);
	}

	double diff_energies_squared;
	double product;
	double res = 0;
	double overlap = 0;
	for(int k = 0; k < dft_path.size(); k++){
		for(int i = 0; i < dft_band_eigenvectors_complex[0].size(); i++){
			for(int j = 0; j < tb_band_eigenvectors[0].size(); j++){
				diff_energies_squared = pow(dft_bands[k][i]-tb_bands[k][j],2);
				overlap = OVERLAP(dft_band_eigenvectors_complex[k][i],tb_band_eigenvectors[k][j]);
				product = diff_energies_squared*overlap;
				res += product;
			}
		}
	}
	return res;}

void nelmin (double residual(band_struct band_structure(vector<double> params, vector<vector<double> > neighbors, vector<vector<double> > next_neighbors, vector<vector<double> > dft_path),
	vector<double> params, vector<vector<double> > neighbors, vector<vector<double> > next_neighbors, vector<vector<double> > dft_path, vector<vector<double> > dft_bands, 
	vector<vector<vector<complex<double> > > > dft_band_eigenvectors_complex), 
	
	band_struct band_structure(vector<double> params, vector<vector<double> > neighbors, vector<vector<double> > next_neighbors, vector<vector<double> > dft_path), 

	int n, vector<double> params, vector<vector<double> > neighbors, 
	vector<vector<double> > next_neighbors, vector<vector<double> > dft_path, vector<vector<double> > dft_bands,
	vector<vector<vector<complex<double> > > > dft_band_eigenvectors_complex, vector<double> xmin, 
    double *ynewlo, double reqmin, vector<double> step, int konvge, int kcount, 
    int *icount, int *numres, int *ifault );
void timestamp ( void );

int main(){

	////////////////////////////////////////////////////////////////////////////////////////
	//
	//  Get data from SCAN DFT
	//
	////////////////////////////////////////////////////////////////////////////////////////

	ifstream PROCAR;
	PROCAR.open("PROCAR");
	int num_bands;
	int num_kpoints;
	int num_ions;
	double dft_energy;
	vector<double> dft_energies;
	vector<vector<double> > dft_bands;
	double kx, ky, kz;
	vector<vector<double> > dft_path;
	vector<vector<double> > dft_eigenvectors;
	vector<vector<vector<double> > >dft_band_eigenvectors;
	double s, py, pz, px, dxy, dyz, dz2, dxz, dx2;
	string throwaway;
	PROCAR >> throwaway >> throwaway >> throwaway;
	PROCAR >> throwaway >> throwaway >> throwaway;
	PROCAR >> num_kpoints >> throwaway >> throwaway >> throwaway >> num_bands >> throwaway >> throwaway >> throwaway >> num_ions;
	PROCAR >> throwaway >> throwaway >> throwaway >> kx >> ky >> kz >> throwaway >> throwaway >> throwaway;
	int counter = 4;

	while(PROCAR){	
		if((counter-4)%265==0){
			//k-points
			if(counter-4 != 0){
				dft_bands.push_back(dft_energies);
				dft_band_eigenvectors.push_back(dft_eigenvectors);
				dft_energies.clear();
				dft_eigenvectors.clear();
				dft_path.push_back(vector<double > {kx,ky,kz});
				PROCAR >> throwaway >> throwaway >> throwaway >> kx >> ky >> kz >> throwaway >> throwaway >> throwaway;
			}
		}
		if(((counter-4)%265)%(num_ions+3)==1){
			//bands
			PROCAR >> throwaway >> throwaway >> throwaway >> throwaway >> dft_energy >> throwaway >> throwaway >> throwaway;
			// cout << dft_energy << endl;
			dft_energies.push_back(dft_energy);

		}
		if((((counter-4)%265)%(num_ions+3)==0) && (counter-4)%265!=0){
			//projection
			PROCAR >> throwaway >> s >> py >> pz >> px >> dxy >> dyz >> dz2 >> dxz >> dx2 >> throwaway;
			// cout << s << " " << py << " " << pz << " " << px << " " << dxy << " " << dyz << " " << dz2 << " " << dxz << " " << dx2 << endl;
			dft_eigenvectors.push_back(vector<double > {s,py,pz,px,dxy,dyz,dz2,dxz,dx2});
		}
		if(!((counter-4)%265==0) && !(((counter-4)%265)%(num_ions+3)==1) && !((((counter-4)%265)%(num_ions+3)==0) && (counter-4)%265!=0)){
			PROCAR >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway;
		}
		counter+=1;
	}
	for(int i = 0; i < num_bands; i++){
		cout << dft_bands[0][i] << endl;
	}
	dft_path.erase(dft_path.begin()+100);
	dft_path.erase(dft_path.begin()+200);
	dft_path.erase(dft_path.begin()+300);
	dft_path.erase(dft_path.begin()+400);
	dft_path.erase(dft_path.begin()+500);
	dft_path.erase(dft_path.begin()+600);
	dft_path.erase(dft_path.begin()+700);
	dft_path.erase(dft_path.begin()+800);
	
	vector<double> ticks;
	ticks.push_back(0.0);
	double tick = 0;
	ticks.push_back(tick);
	for(int i = 1; i <dft_path.size(); i++){
		if(pow(pow(dft_path[i][0]-dft_path[i-1][0],2)+pow(dft_path[i][1]-dft_path[i-1][1],2)+pow(dft_path[i][2]-dft_path[i-1][2],2),0.5) < 1.0){
			tick += pow(pow(dft_path[i][0]-dft_path[i-1][0],2)+pow(dft_path[i][1]-dft_path[i-1][1],2)+pow(dft_path[i][2]-dft_path[i-1][2],2),0.5);
		}
		else{
			tick += 0.000000000000001;
		}

		ticks.push_back(tick);
	}

	dft_band_eigenvectors.erase(dft_band_eigenvectors.begin()+100);
	dft_band_eigenvectors.erase(dft_band_eigenvectors.begin()+200);
	dft_band_eigenvectors.erase(dft_band_eigenvectors.begin()+300);
	dft_band_eigenvectors.erase(dft_band_eigenvectors.begin()+400);
	dft_band_eigenvectors.erase(dft_band_eigenvectors.begin()+500);
	dft_band_eigenvectors.erase(dft_band_eigenvectors.begin()+600);
	dft_band_eigenvectors.erase(dft_band_eigenvectors.begin()+700);
	dft_band_eigenvectors.erase(dft_band_eigenvectors.begin()+800);

	dft_bands.erase(dft_bands.begin()+100);
	dft_bands.erase(dft_bands.begin()+200);
	dft_bands.erase(dft_bands.begin()+300);
	dft_bands.erase(dft_bands.begin()+400);
	dft_bands.erase(dft_bands.begin()+500);
	dft_bands.erase(dft_bands.begin()+600);
	dft_bands.erase(dft_bands.begin()+700);
	dft_bands.erase(dft_bands.begin()+800);

	vector<vector<double> > dft_bands_T(dft_bands[0].size(),vector<double> (dft_bands.size()));

	for(int i = 0; i < dft_bands.size(); i++){
		for(int j = 0; j < dft_bands[0].size(); j++){
			dft_bands_T[j][i] = dft_bands[i][j];
		}
	}

	vector<vector<vector<complex<double> > > > dft_band_eigenvectors_complex;
	vector<vector<complex<double> > > dft_band_eigenvector_complexs;
	vector<complex<double> > dft_band_eigenvector_complex;
	for(int i = 0; i < dft_band_eigenvectors.size(); i++){
		dft_band_eigenvector_complexs.clear();
		for(int j = 0; j < dft_band_eigenvectors[0].size(); j++){
			dft_band_eigenvector_complex.clear();
			for(int k = 0; k < dft_band_eigenvectors[0][0].size(); k++){
				dft_band_eigenvector_complex.push_back(complex<double>(dft_band_eigenvectors[i][j][k],0.0));
			}
			dft_band_eigenvector_complexs.push_back(dft_band_eigenvector_complex);
		}
		dft_band_eigenvectors_complex.push_back(dft_band_eigenvector_complexs);
	}

	///////////////////////////////////////////////////////////////////////////////////////
	//// Plot DFT Bandstructure
	///////////////////////////////////////////////////////////////////////////////////////
	// ofstream band_structure_data_dft;
	// band_structure_data_dft.open("BAND_STRUCTURE_DFT.dat");
	// counter = 0;
	// for(auto i : dft_bands_T){
	// 	for(auto j : i){
	// 		band_structure_data_dft << ticks.at(counter) << "\t" << j << endl;
	// 		counter +=1;
	// 	}
	// 	counter = 0;
	// 	band_structure_data_dft << endl;
	// }
	// band_structure_data_dft.close();
	// string command = "python plot_bands_c++.py";
	// const char* com = command.c_str();
	// system(com);

	////////////////////////////////////////////////////////////////////////////////////////
	//**************************************************************************************
	////////////////////////////////////////////////////////////////////////////////////////

	double a = 5.431; 
	vector<double> n1{0.25*a,-0.25*a,0.25*a};
	vector<double> n2{0.25*a,0.25*a,-0.25*a};
	vector<double> n3{-0.25*a,0.25*a,0.25*a};
	vector<double> n4{-0.25*a,-0.25*a,-0.25*a};
	vector<double> nn1{0.5*a,0.0,-0.5*a};
	vector<double> nn2{0.0,0.5*a,-0.5*a};
	vector<double> nn3{-0.5*a,0,-0.5*a};
	vector<double> nn4{0.0,-0.5*a,-0.5*a};
	vector<double> nn5{0.5*a,0.5*a,0.0};
	vector<double> nn6{-0.5*a,0.5*a,0.0};
	vector<double> nn7{-0.5*a,-0.5*a,0.0};
	vector<double> nn8{0.5*a,-0.5*a,0.0};
	vector<double> nn9{0.5*a,0,0.5*a};
	vector<double> nn10{0.0,0.5*a,0.5*a};
	vector<double> nn11{-0.5*a,0,0.5*a};
	vector<double> nn12{0.0,-0.5*a,0.5*a};

	vector<vector<double> > neighbors{n1,n2,n3,n4};
	vector<vector<double> > next_neighbors{nn1,nn2,nn3,nn4,nn5,nn6,nn7,nn8,nn9,nn10,nn11,nn12};


	// vector<double> params{-4.03,2.9,-2.035, 1.47, -1.45, 4.1825, -1.425657743, 1.212435565, 2.045116134,-1.204159458};
	vector<double> params{1,1,1, 1, 1, 1,1, 1, 1,1};
	
	// int i;
	// int icount;
	// int ifault;
	// int kcount;
	// int konvge;
	// int n;
	// int numres;
	// double reqmin;
	// double ynewlo;

	// n = 10;
	// vector<double> step(n);
 // 	vector<double> xmin(n);
	// konvge = 1;
 // 	kcount = 1;
 // 	reqmin = 1.0E-08;
	// cout << "\n";
	// cout << "  Starting point X:\n";
	// cout << "\n";
	// for ( i = 0; i < n; i++ )
	// {
	// cout << "  " << setw(14) << params[i] << "\n";
	// }
	// ynewlo = residual(band_structure, params, neighbors, next_neighbors, dft_path, dft_bands, 
	// 	dft_band_eigenvectors_complex);

	// cout << "\n";
	// cout << "  F(X) = " << ynewlo << "\n";

	// nelmin ( residual, band_structure, n, params, neighbors, next_neighbors, dft_path, dft_bands, dft_band_eigenvectors_complex,xmin, &ynewlo, reqmin, step,
	// konvge, kcount, &icount, &numres, &ifault );

	// cout << "\n";
	// cout << "  Return code IFAULT = " << ifault << "\n";
	// cout << "\n";
	// cout << "  Estimate of minimizing value X*:\n";
	// cout << "\n";
	// for ( i = 0; i < n; i++ )
	// {
	// cout << "  " << setw(14) << xmin[i] << "\n";
	// }

	// cout << "\n";
	// cout << "  F(X*) = " << ynewlo << "\n";

	// cout << "\n";
	// cout << "  Number of iterations = " << icount << "\n";
	// cout << "  Number of restarts =   " << numres << "\n";

	///////////////////////////////////////////////////////////////////////////////////////
	//// Plot TB Bandstructure
	///////////////////////////////////////////////////////////////////////////////////////

	// band_struct BAND_STRUCTURE = band_structure(xmin, neighbors, next_neighbors, dft_path);
	// vector<vector<complex<double> > > BAND_VALUES = BAND_STRUCTURE.bandvalues;
	// vector<vector<vector<complex<double> > > > tb_band_eigenvectors = BAND_STRUCTURE.bandvectors;

	// vector<double> band;
	// vector<vector<double> > bands(BAND_VALUES.size(),vector<double> (BAND_VALUES[0].size()));
	// vector<vector<double> > bands_T(BAND_VALUES[0].size(),vector<double> (BAND_VALUES.size()));
	// vector<vector<double> > tb_bands_T(BAND_VALUES[0].size(),vector<double> (BAND_VALUES.size()));
	// vector<vector<double> > tb_bands;

	// for(int i = 0; i < BAND_VALUES.size(); i++){
	// 	for(int j = 0; j < BAND_VALUES[0].size(); j++){
	// 		bands[i][j] = real(BAND_VALUES[i][j]);
	// 		bands_T[j][i] = bands[i][j];
	// 	}
	// }

	// for(int i = 0; i < bands.size(); i++){
	// 	band = bands[i];
	// 	sort(band.begin(),band.end());
	// 	tb_bands.push_back(band);
	// }
	// for(int i = 0; i < tb_bands.size(); i++){
	// 	for(int j = 0; j < tb_bands[0].size(); j++){
	// 		tb_bands_T[j][i] = tb_bands[i][j];
	// 	}
	// }
	// ofstream band_structure_data_tb;
	// band_structure_data_tb.open("BAND_STRUCTURE_TB.dat");
	// counter = 0;
	// for(auto i : tb_bands_T){
	// 	for(auto j : i){
	// 		band_structure_data_tb << ticks.at(counter) << "\t" << j << endl;
	// 		counter +=1;
	// 	}
	// 	counter = 0;
	// 	band_structure_data_tb << endl;
	// }
	// band_structure_data_tb.close();
	// string command = "python plot_bands_c++.py";
	// const char* com = command.c_str();
	// system(com);

	////////////////////////////////////////////////////////////////////////////////////////
	//**************************************************************************************
	////////////////////////////////////////////////////////////////////////////////////////

	return 0;
}

//****************************************************************************

void nelmin (double residual(band_struct band_structure(vector<double> params, vector<vector<double> > neighbors, vector<vector<double> > next_neighbors, vector<vector<double> > dft_path),
	vector<double> params, vector<vector<double> > neighbors, vector<vector<double> > next_neighbors, vector<vector<double> > dft_path, vector<vector<double> > dft_bands, 
	vector<vector<vector<complex<double> > > > dft_band_eigenvectors_complex), 
	
	band_struct band_structure(vector<double> params, vector<vector<double> > neighbors, vector<vector<double> > next_neighbors, vector<vector<double> > dft_path), 

	int n, vector<double> params, vector<vector<double> > neighbors, 
	vector<vector<double> > next_neighbors, vector<vector<double> > dft_path, vector<vector<double> > dft_bands,
	vector<vector<vector<complex<double> > > > dft_band_eigenvectors_complex, vector<double> xmin, 
    double *ynewlo, double reqmin, vector<double> step, int konvge, int kcount, 
    int *icount, int *numres, int *ifault )

//****************************************************************************
//
//  Purpose:
//
//    NELMIN minimizes a function using the Nelder-Mead algorithm.
//
//  Discussion:
//
//    This routine seeks the minimum value of a user-specified function.
//
//    Simplex function minimisation procedure due to Nelder+Mead(1965),
//    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
//    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
//    25, 97) and Hill(1978, 27, 380-2)
//
//    The function to be minimized must be defined by a function of
//    the form
//
//      function residual ( x, f )
//      double residual
//      double x(*)
//
//    and the name of this subroutine must be declared EXTERNAL in the
//    calling routine and passed as the argument FN.
//
//    This routine does not include a termination test using the
//    fitting of a quadratic surface.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 February 2008
//
//  Author:
//
//    Original FORTRAN77 version by R ONeill.
//    C++ version by John Burkardt.
//
//  Reference:
//
//    John Nelder, Roger Mead,
//    A simplex method for function minimization,
//    Computer Journal,
//    Volume 7, 1965, pages 308-313.
//
//    R ONeill,
//    Algorithm AS 47:
//    Function Minimization Using a Simplex Procedure,
//    Applied Statistics,
//    Volume 20, Number 3, 1971, pages 338-345.
//
//  Parameters:
//
//    Input, double FN ( double x[] ), the name of the routine which evaluates
//    the function to be minimized.
//
//    Input, int N, the number of variables.
//
//    Input/output, double START[N].  On input, a paramsing point
//    for the iteration.  On output, this data may have been overwritten.
//
//    Output, double XMIN[N], the coordinates of the point which
//    is estimated to minimize the function.
//
//    Output, double YNEWLO, the minimum value of the function.
//
//    Input, double REQMIN, the terminating limit for the variance
//    of function values.
//
//    Input, double STEP[N], determines the size and shape of the
//    initial simplex.  The relative magnitudes of its elements should reflect
//    the units of the variables.
//
//    Input, int KONVGE, the convergence check is carried out 
//    every KONVGE iterations.
//
//    Input, int KCOUNT, the maximum number of function 
//    evaluations.
//
//    Output, int *ICOUNT, the number of function evaluations 
//    used.
//
//    Output, int *NUMRES, the number of restarts.
//
//    Output, int *IFAULT, error indicator.
//    0, no errors detected.
//    1, REQMIN, N, or KONVGE has an illegal value.
//    2, iteration terminated because KCOUNT was exceeded without convergence.
//
{
  double ccoeff = 0.5;
  double del;
  double dn;
  double dnn;
  double ecoeff = 2.0;
  double eps = 0.001;
  int i;
  int ihi;
  int ilo;
  int j;
  int jcount;
  int l;
  int nn;
  double rcoeff = 1.0;
  double rq;
  double x;
  double y2star;
  double ylo;
  double ystar;
  double z;
//
//  Check the input parameters.
//
  if ( reqmin <= 0.0 )
  {
    *ifault = 1;
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    return;
  }

  vector<double> p(n*(n+1));
  vector<double> pstar(n);
  vector<double> p2star(n);
  vector<double> pbar(n);
  vector<double> y(n+1);

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;
//
//  Initial or restarted loop.
//
  for ( ; ; )
  {
    for ( i = 0; i < n; i++ )
    { 
      p[i+n*n] = params[i];
    }
    y[n] = residual (band_structure, params, neighbors, next_neighbors,
	dft_path, dft_bands, dft_band_eigenvectors_complex);
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = params[j];
      params[j] = params[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = params[i];
      }
      y[j] = residual (band_structure, params, neighbors, next_neighbors,
	  dft_path, dft_bands, dft_band_eigenvectors_complex);
      *icount = *icount + 1;
      params[j] = x;
    }
//                    
//  The simplex construction is complete.
//                    
//  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
//  the vertex of the simplex to be replaced.
//                
    ylo = y[0];
    ilo = 0;

    for ( i = 1; i < nn; i++ )
    {
      if ( y[i] < ylo )
      {
        ylo = y[i];
        ilo = i;
      }
    }
//
//  Inner loop.
//
    for ( ; ; )
    {
      if ( kcount <= *icount )
      {
        break;
      }
      *ynewlo = y[0];
      ihi = 0;

      for ( i = 1; i < nn; i++ )
      {
        if ( *ynewlo < y[i] )
        {
          *ynewlo = y[i];
          ihi = i;
        }
      }
//
//  Calculate PBAR, the centroid of the simplex vertices
//  excepting the vertex with Y value YNEWLO.
//
      for ( i = 0; i < n; i++ )
      {
        z = 0.0;
        for ( j = 0; j < nn; j++ )
        { 
          z = z + p[i+j*n];
        }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
//
//  Reflection through the centroid.
//
      for ( i = 0; i < n; i++ )
      {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }
      ystar = residual (band_structure, pstar, neighbors, next_neighbors,
		dft_path, dft_bands, dft_band_eigenvectors_complex);
      *icount = *icount + 1;
//
//  Successful reflection, so extension.
//
      if ( ystar < ylo )
      {
        for ( i = 0; i < n; i++ )
        {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
        y2star =  residual (band_structure, p2star, neighbors, next_neighbors,
		dft_path, dft_bands, dft_band_eigenvectors_complex);
        *icount = *icount + 1;
//
//  Check extension.
//
        if ( ystar < y2star )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Retain extension or contraction.
//
        else
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = p2star[i];
          }
          y[ihi] = y2star;
        }
      }
//
//  No extension.
//
      else
      {
        l = 0;
        for ( i = 0; i < nn; i++ )
        {
          if ( ystar < y[i] )
          {
            l = l + 1;
          }
        }

        if ( 1 < l )
        {
          for ( i = 0; i < n; i++ )
          {
            p[i+ihi*n] = pstar[i];
          }
          y[ihi] = ystar;
        }
//
//  Contraction on the Y(IHI) side of the centroid.
//
        else if ( l == 0 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
          y2star = residual (band_structure, p2star, neighbors, next_neighbors,
		dft_path, dft_bands, dft_band_eigenvectors_complex);
          *icount = *icount + 1;
//
//  Contract the whole simplex.
//
          if ( y[ihi] < y2star )
          {
            for ( j = 0; j < nn; j++ )
            {
              for ( i = 0; i < n; i++ )
              {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
              y[j] = residual (band_structure, xmin, neighbors, next_neighbors,
		dft_path, dft_bands, dft_band_eigenvectors_complex);
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;

            for ( i = 1; i < nn; i++ )
            {
              if ( y[i] < ylo )
              {
                ylo = y[i];
                ilo = i;
              }
            }
            continue;
          }
//
//  Retain contraction.
//
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
//
//  Contraction on the reflection side of the centroid.
//
        else if ( l == 1 )
        {
          for ( i = 0; i < n; i++ )
          {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
          y2star = residual (band_structure, p2star, neighbors, next_neighbors,
		dft_path, dft_bands, dft_band_eigenvectors_complex);
          *icount = *icount + 1;
//
//  Retain reflection?
//
          if ( y2star <= ystar )
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
          else
          {
            for ( i = 0; i < n; i++ )
            {
              p[i+ihi*n] = pstar[i];
            }
            y[ihi] = ystar;
          }
        }
      }
//
//  Check if YLO improved.
//
      if ( y[ihi] < ylo )
      {
        ylo = y[ihi];
        ilo = ihi;
      }
      jcount = jcount - 1;

      if ( 0 < jcount )
      {
        continue;
      }
//
//  Check to see if minimum reached.
//
      if ( *icount <= kcount )
      {
        jcount = konvge;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + y[i];
        }
        x = z / dnn;

        z = 0.0;
        for ( i = 0; i < nn; i++ )
        {
          z = z + pow ( y[i] - x, 2 );
        }

        if ( z <= rq )
        {
          break;
        }
      }
    }
//
//  Factorial tests to check that YNEWLO is a local minimum.
//
    for ( i = 0; i < n; i++ )
    {
      xmin[i] = p[i+ilo*n];
    }
    *ynewlo = y[ilo];

    if ( kcount < *icount )
    {
      *ifault = 2;
      break;
    }

    *ifault = 0;

    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
      z = residual (band_structure, xmin, neighbors, next_neighbors,
		dft_path, dft_bands, dft_band_eigenvectors_complex);
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
      z = residual (band_structure, xmin, neighbors, next_neighbors,
		dft_path, dft_bands, dft_band_eigenvectors_complex);
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
      break;
    }
//
//  Reparams the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      params[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }

  return;
}
//****************************************************************************

void timestamp ( void )

//****************************************************************************
//
//  Purpose:
//
//    TIMESTAMP prints the current YMDHMS date as a time stamp.
//
//  Example:
//
//    31 May 2001 09:45:54 AM
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    24 September 2003
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    None
//
{
# define TIME_SIZE 40

  static char time_buffer[TIME_SIZE];
  const struct tm *tm;
  size_t len;
  time_t now;

  now = time ( NULL );
  tm = localtime ( &now );

  len = strftime ( time_buffer, TIME_SIZE, "%d %B %Y %I:%M:%S %p", tm );

  cout << time_buffer << "\n";

  return;
# undef TIME_SIZE
}