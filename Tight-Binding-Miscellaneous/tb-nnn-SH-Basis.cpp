
# include <cstdlib>
# include <iostream>
# include <string>
# include <iomanip>
# include <cmath>
# include <fstream>
# include <vector>
# include <cstdlib>
# include <omp.h>
using namespace std;
const double a = 5.431;
const int N = 12;
const int P = 300;
// P = 894
struct comp{double real = 0; double imag = 0; };
struct dot{double p1[4] = {0}; double p2[12] = {0}; };
struct phi{comp factors[18] = {}; };
struct proj_bandvecs_and_vals{double bandvalues[N] = {0}; comp projected_vec[N][N] = {}; };
struct eigen{double evalues[N] = {0}; double evectors[N][N] = {}; };
struct bandvecs_and_vals{double bandvalues[N] = {0}; double bandvectors[N][N] = {0}; };

comp conj(comp number){
	// cout << "conj" << endl;
	comp conjugate = {};
	conjugate.real = number.real;
	conjugate.imag = -1.0*number.imag;
	return conjugate;
}

comp mult_comp(comp c_number, double number){
	// cout << "mult_comp" << endl;
	comp mult = {};
	mult.real = c_number.real*number;
	mult.imag = c_number.imag*number;
	return mult;
}

comp mult_two_comps(comp number1, comp number2){
	// cout << "mult_two_comps" << endl;
	comp new_number = {};
	new_number.real = number1.real*number2.real - number1.imag*number2.imag;
	new_number.imag = number1.real*number2.imag + number1.imag*number2.real;
	return new_number;
}

comp add_comp(comp number1, comp number2){
	// cout << "add_comp" << endl;
	comp add = {};
	add.real = number1.real + number2.real;
	add.imag = number1.imag + number2.imag;
	return add;
}

double norm(comp number){
	// cout << "norm" << endl;
	return pow(pow(number.real,2) + pow(number.imag,2),0.5);
}

double s_orb(double x, double y, double z){
	return 0.5*pow(M_PI,-0.5);
}
double py_orb(double x, double y, double z){
	return pow(3/(4*M_PI),0.5)*y/pow(x*x+y*y+z*z,0.5);
}
double pz_orb(double x, double y, double z){
	return pow(3/(4*M_PI),0.5)*z/pow(x*x+y*y+z*z,0.5);
}
double px_orb(double x, double y, double z){
	return pow(3/(4*M_PI),0.5)*x/pow(x*x+y*y+z*z,0.5);
}
double dxy_orb(double x, double y, double z){
	return 0.5*pow(15/M_PI,0.5)*x*y/(x*x+y*y+z*z);
}
double dyz_orb(double x, double y, double z){
	return 0.5*pow(15/M_PI,0.5)*y*z/(x*x+y*y+z*z);
}
double dz_2_orb(double x, double y, double z){
	return 0.25*pow(5/M_PI,0.5)*(-x*x-y*y+2*z*z)/(x*x+y*y+z*z);
}
double dxz_orb(double x, double y, double z){
	return 0.5*pow(15/M_PI,0.5)*z*x/(x*x+y*y+z*z);
}
double dx_2_y_2_orb(double x, double y, double z){
	return 0.25*pow(15/M_PI,0.5)*(x*x-y*y)/(x*x+y*y+z*z);
}
double fyz_2_orb(double x, double y, double z){
	return 0.25*pow(21/(2*M_PI),0.5)*(y*(4*z*z-x*x-y*y))/pow(x*x+y*y+z*z,1.5);
}
double fz_3_orb(double x, double y, double z){
	return 0.25*pow(7/M_PI,0.5)*(z*(2*z*z-3*x*x-3*y*y))/pow(x*x+y*y+z*z,1.5);
}
double fxz_2_orb(double x, double y, double z){
	return 0.25*pow(21/(2*M_PI),0.5)*(x*(4*z*z-x*x-y*y))/pow(x*x+y*y+z*z,1.5);
}

typedef double (*Functions) (double x, double y, double z);

proj_bandvecs_and_vals projection_gauss_10_pt_quad(bandvecs_and_vals bandvecs_and_values(double params[10], double neighbors[4][3], 
	double next_neighbors[12][3], double kx, double ky, double kz, bool get_vecs),
	double f0(double x, double y, double z), 
	double f1(double x, double y, double z),
	double f2(double x, double y, double z),
	double f3(double x, double y, double z),
	double f4(double x, double y, double z),
	double f5(double x, double y, double z),
	double f6(double x, double y, double z),
	double f7(double x, double y, double z),
	double f8(double x, double y, double z),
	double f9(double x, double y, double z),
	double f10(double x, double y, double z),
	double f11(double x, double y, double z), 
	double params[10], double neighbors[4][3], double next_neighbors[12][3],
	double ax, double bx, double ay, double by, double az, double bz, double akx, double bkx, double aky, double bky, double akz, double bkz, 
	double kx, double ky, double kz){
	proj_bandvecs_and_vals TB_Proj = {};
    Functions functions[] = 
    {
    	f0,f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11
    };
	double points[10]  = {-0.1488743389816312108848,
			0.1488743389816312108848,
		   -0.4333953941292471907993,
		    0.4333953941292471907993,
		   -0.6794095682990244062343,
		    0.6794095682990244062343,
		   -0.8650633666889845107321,
		    0.8650633666889845107321,
		   -0.973906528517171720078,
		    0.973906528517171720078};
	double weights[10] = {0.2955242247147528701739,
		   0.2955242247147528701739,
		   0.269266719309996355091,
		   0.269266719309996355091,
		   0.2190863625159820439955,
		   0.2190863625159820439955,
		   0.149451349150580593146,
		   0.149451349150580593146,
		   0.0666713443086881375936,
		   0.0666713443086881375936};
	double half_diff_x = (bx-ax)*0.5;
	double half_diff_y = (by-ay)*0.5;
	double half_diff_z = (bz-az)*0.5;
	double half_diff_kx = (bkx-akx)*0.5;
	double half_diff_ky = (bky-aky)*0.5;
	double half_diff_kz = (bkz-akz)*0.5;
	double half_sum_x = (bx+ax)*0.5;
	double half_sum_y = (by+ay)*0.5;
	double half_sum_z = (bz+az)*0.5;
	double half_sum_kx = (bkx+akx)*0.5;
	double half_sum_ky = (bky+aky)*0.5;
	double half_sum_kz = (bkz+akz)*0.5;
	double partial_integral_real = 0;
	double partial_integral_imag = 0;

	double real_part[N][N] = {0};
	double imag_part[N][N] = {0};

	/////
	///// read in values of the fourier transformed spherical harmonics
	///// do 11 point gaussian quadrature instead of 10, because you want to include high symmetry k points that have zeros
	///// 
	#pragma omp parallel
	{					
	#pragma omp for collapse(3), reduction(+: real_part[:12][:12], imag_part[:12][:12])
	for(int l = 0; l < 10; l++){
		for(int m = 0; m < 10; m++){
			for(int n = 0; n < 10; n++){
				double end = half_diff_x*half_diff_y*half_diff_z*half_diff_kx*half_diff_ky*half_diff_kz*weights[l]*weights[m]*weights[n];
				bandvecs_and_vals TB = bandvecs_and_values(params, neighbors, next_neighbors, half_diff_kx*points[l]+half_sum_kx, half_diff_ky*points[m]+half_sum_ky, half_diff_kz*points[n]+half_sum_kz,true);		
				for(int i = 0; i < 10; i++){
					for(int j = 0; j < 10; j++){
						for(int k = 0; k < 10; k++){
							for(int s = 0; s < N; s++){
								for(int e = 0; e < N; e++){	
									real_part[s][e] += TB.bandvectors[s][e]*functions[e](half_diff_x*points[i]+half_sum_x,half_diff_y*points[j]+half_sum_y,half_diff_z*points[k]+half_sum_z)
									*weights[i]*weights[j]*weights[k]*
									cos((half_diff_kx*points[l]+half_sum_kx)*(half_diff_x*points[i]+half_sum_x)+
										(half_diff_ky*points[m]+half_sum_ky)*(half_diff_y*points[j]+half_sum_y)+
										(half_diff_kz*points[n]+half_sum_kz)*(half_diff_z*points[k]+half_sum_z))*end;

									imag_part[s][e] += TB.bandvectors[s][e]*functions[e](half_diff_x*points[i]+half_sum_x,half_diff_y*points[j]+half_sum_y,half_diff_z*points[k]+half_sum_z)
									*weights[i]*weights[j]*weights[k]*
									sin((half_diff_kx*points[l]+half_sum_kx)*(half_diff_x*points[i]+half_sum_x)+
										(half_diff_ky*points[m]+half_sum_ky)*(half_diff_y*points[j]+half_sum_y)+
										(half_diff_kz*points[n]+half_sum_kz)*(half_diff_z*points[k]+half_sum_z))*end;
									}
								}
							}
						}
					}
				}
			}	
		}
	}
	bandvecs_and_vals TB = bandvecs_and_values(params, neighbors, next_neighbors, kx, ky, kz, false);
	double normalization = 0;
	for(int s = 0; s < N; s++){
		TB_Proj.bandvalues[s] = TB.bandvalues[s];
		normalization = 0;
			for(int e = 0; e < N; e++){	
				normalization += pow(real_part[s][e],2) + pow(imag_part[s][e],2);
			}
			for(int e = 0; e < N; e++){	
				TB_Proj.projected_vec[s][e].real = real_part[s][e]/pow(normalization,0.5);
				TB_Proj.projected_vec[s][e].imag = imag_part[s][e]/pow(normalization,0.5);	
			}
		}
	return TB_Proj;
}

// double Fourier_real(double gauss(double f(double x, double y, double z), double ax, double bx, double ay, double by, double az, double bz), double Y(double x, double y, double z), double K[3]){
// 	double fourier = 0;
// 	double mult_phase(double )
// 	fourier = pow(2*M_PI,-0.5)*gauss(,);
// 	return fourier;
// }

dot dot_product(int vec1[3], int vec2[3], int typeint, double neighbors[4][3], double next_neighbors[12][3]){
	// cout << "dot_product" << endl;
	dot p = {};
	if(typeint == 1){
		for(int i = 0; i< 4; i++){
			p.p1[i] = (vec1[0]*neighbors[i][0]+vec1[1]*neighbors[i][1]+vec1[2]*neighbors[i][2])
		        /(pow(neighbors[i][0]*neighbors[i][0]+neighbors[i][1]*neighbors[i][1]
		        +neighbors[i][2]*neighbors[i][2],0.5)); 
		    p.p2[i] =(vec1[0]*next_neighbors[i][0]+vec1[1]*next_neighbors[i][1]+vec1[2]*next_neighbors[i][2])
				/(pow(next_neighbors[i][0]*next_neighbors[i][0]+next_neighbors[i][1]*next_neighbors[i][1]
			    +next_neighbors[i][2]*next_neighbors[i][2],0.5));
		}
	}
	else{
		for(int i = 0; i< 4; i++){
			p.p1[i] = ((vec1[0]*neighbors[i][0]+vec1[1]*neighbors[i][1]+vec1[2]*neighbors[i][2])
	            /(pow(neighbors[i][0]*neighbors[i][0]+neighbors[i][1]*neighbors[i][1]
	            +neighbors[i][2]*neighbors[i][2],0.5)))*((vec2[0]*neighbors[i][0]+vec2[1]*neighbors[i][1]
			    +vec2[2]*neighbors[i][2])/(pow(neighbors[i][0]*neighbors[i][0]+neighbors[i][1]*neighbors[i][1]
			    +neighbors[i][2]*neighbors[i][2],0.5)));
			
			p.p2[i] = ((vec1[0]*next_neighbors[i][0]+vec1[1]*next_neighbors[i][1]+vec1[2]*next_neighbors[i][2])
	            /(pow(next_neighbors[i][0]*next_neighbors[i][0]+next_neighbors[i][1]*next_neighbors[i][1]
	            +next_neighbors[i][2]*next_neighbors[i][2],0.5)))*((vec2[0]*next_neighbors[i][0]+vec2[1]*next_neighbors[i][1]
			    +vec2[2]*next_neighbors[i][2])/(pow(next_neighbors[i][0]*next_neighbors[i][0]+next_neighbors[i][1]*next_neighbors[i][1]
			    +next_neighbors[i][2]*next_neighbors[i][2],0.5)));
		}
	}


	return p;			
}

phi phase(double kx, double ky, double kz, double neighbors[4][3], double next_neighbors[12][3]){
	// cout << "phase" << endl;
	int x[3] = {1,0,0};
	int y[3] = {0,1,0};
	int z[3] = {0,0,1}; 
	comp nnp[4] = {};
	double nnd[9][4] = {0};
	comp nnnp[12] = {};
	double nnnd[9][12] = {0};
	int zero_vec[3] = {0};
	phi p = {};
	for(int i = 0; i<4;i++){
		nnp[i].real = cos(kx*neighbors[i][0]+ky*neighbors[i][1]+kz*neighbors[i][2]);
		nnp[i].imag = sin(kx*neighbors[i][0]+ky*neighbors[i][1]+kz*neighbors[i][2]);
	}
	dot temp0 = {}, temp1 = {}, temp2 = {}, temp3 = {}, temp4 = {}, temp5 = {}, temp6 = {}, temp7 = {}, temp8 = {};
	temp0 = dot_product(x,zero_vec,1,neighbors,next_neighbors);
	temp1 = dot_product(y,zero_vec,1,neighbors,next_neighbors);
	temp2 = dot_product(z,zero_vec,1,neighbors,next_neighbors);
	temp3 = dot_product(x,x,2,neighbors,next_neighbors);
	temp4 = dot_product(y,y,2,neighbors,next_neighbors);
	temp5 = dot_product(z,z,2,neighbors,next_neighbors);
	temp6 = dot_product(x,y,2,neighbors,next_neighbors);
	temp7 = dot_product(x,z,2,neighbors,next_neighbors);
	temp8 = dot_product(y,z,2,neighbors,next_neighbors);
	for(int i = 0; i < 4; i++){
		nnd[0][i] = temp0.p1[i]; 
		nnd[1][i] = temp1.p1[i];
		nnd[2][i] = temp2.p1[i];
		nnd[3][i] = temp3.p1[i];
		nnd[4][i] = temp4.p1[i];
		nnd[5][i] = temp5.p1[i];
		nnd[6][i] = temp6.p1[i];
		nnd[7][i] = temp7.p1[i];
		nnd[8][i] = temp8.p1[i];
	}
	
	for(int i = 0; i<12;i++){
		nnnp[i].real = cos(kx*next_neighbors[i][0]+ky*next_neighbors[i][1]+kz*next_neighbors[i][2]);
		nnnp[i].imag = sin(kx*next_neighbors[i][0]+ky*next_neighbors[i][1]+kz*next_neighbors[i][2]);
	}

	for(int i = 0; i < 12; i++){
		nnnd[0][i] = temp0.p2[i]; 
		nnnd[1][i] = temp1.p2[i];
		nnnd[2][i] = temp2.p2[i];
		nnnd[3][i] = temp3.p2[i];
		nnnd[4][i] = temp4.p2[i];
		nnnd[5][i] = temp5.p2[i];
		nnnd[6][i] = temp6.p2[i];
		nnnd[7][i] = temp7.p2[i];
		nnnd[8][i] = temp8.p2[i];
	}

	for(int l = 0; l < 9; l++){
		p.factors[l].real = 0;
		p.factors[l].imag = 0;
		for(int m = 0; m < 4; m++){
			p.factors[l].real += nnp[m].real*nnd[l][m];
			p.factors[l].imag += nnp[m].imag*nnd[l][m];
		}
	}
	for(int l = 0; l < 9; l++){
		p.factors[l+9].real = 0;
		p.factors[l+9].imag = 0;
		for(int m = 0; m < 12; m++){
			p.factors[l+9].real += nnnp[m].real*nnnd[l][m];
			p.factors[l+9].imag += nnnp[m].imag*nnnd[l][m];
		}
	}
	return p;
}

eigen Householder_with_QR(comp hamiltonian[N][N], bool Vecs){
	// cout << "Householder_with_QR" << endl;
	eigen Eigen = {};
	double Ak[N][N] = {0};
	vector<vector<vector<double> > > Qk;
	vector<vector<double> > Identity(N,vector<double>(N));
	comp b[N] = {};
	comp c[N-1] = {};
	double gamma = 0;
	double phi = 0;
	double alpha = 0;
	comp rho = {};
	comp beta = {};
	comp u[N] = {};
	comp t[N] = {};
	comp r[N] = {};
	double G[N-1][N][N] = {0.0};
	int row = 0;
	int column = 0;
	double temp[N][N] = {0.0};
	double temp2[N][N] = {0.0};
	int sum_diags = 0;
	int counter = 0;
	double temporary = 0;
	int i = 0;
	int j = 0;
	int k = 0;
	int s = 0;
	double sign = 0;
	double delta = 0;
	double mu = 0;
	for(k = 0; k<N-2;k++){
		b[k] = hamiltonian[k][k];
		gamma = 0;
		for(i = k + 1; i < N; i++){
			if(norm(hamiltonian[i][k]) > gamma){
				gamma = norm(hamiltonian[i][k]);
			}
		phi = 0;
		if(gamma < 0.000000001){
			c[k].real = 0;
			c[k].imag = 0;
		}
		else{
			alpha = 0;
			if(norm(hamiltonian[k+1][k]) < 0.000000001){
				beta.real = 1;
				beta.imag = 0;
			}
			else{
				beta = mult_comp(hamiltonian[k+1][k],1/norm(hamiltonian[k+1][k]));
			}
			
			for( i = k + 1; i < N; i++){
				u[i] = mult_comp(hamiltonian[i][k],1/gamma);
				alpha += pow(norm(u[i]),2);
			}

			alpha = pow(alpha, 0.5);
			c[k].real = -1*alpha;
			phi = 1/(alpha*(alpha+norm(u[k+1])));	
			u[k+1] = add_comp(mult_comp(beta, alpha),u[k+1]);
			rho.real = 0;
			rho.imag = 0;

			for( s = k + 1; s < N; s++){
				t[s].real = 0;
				t[s].imag = 0;
				for( j = k+1; j<s; j++){
					t[s] = add_comp(t[s],mult_two_comps(u[j],hamiltonian[s][j]));
				}
				for( j = s+1; j<N; j++){
					t[s] = add_comp(t[s],mult_two_comps(u[j],hamiltonian[j][s]));
				}
				t[s] = mult_comp(t[s],phi);
				rho = add_comp(rho,mult_two_comps(conj(u[s]),t[s]));
			}

			for( s = k + 1; s < N; s++){
				r[s] = mult_two_comps(mult_comp(rho, phi/2),u[s]);
			}
			for( i = k + 1; i <N; i++){
				for( j = k + 1; j < i; j++){
					hamiltonian[i][j] = add_comp(hamiltonian[i][j],add_comp(add_comp(mult_two_comps(mult_comp(t[i],-1),conj(u[j])),mult_two_comps(mult_comp(conj(t[j]),-1),u[i])),
					add_comp(mult_two_comps(u[i],conj(r[j])),mult_two_comps(r[i],conj(u[j])))));
				}
			}
		}
		}
	}
	c[N-2] = hamiltonian[N-1][N-2];
	b[N-2] = hamiltonian[N-2][N-2];
	b[N-1] = hamiltonian[N-1][N-1];

	for( i = 0; i < N; i++){
		for( j = 0; j < N; j++){
			if(i == j){
				Ak[i][j] = b[i].real ;
			}
			else if(j == i + 1){
				if (abs(c[j].imag) < 0.000000000001){
					if (abs(c[j].real) < 0.000000000001){
						Ak[i][j] = 0.0;
					}
					else{
						Ak[i][j] = c[j].real;
					}
				}
				else{
					Ak[i][j] = c[j].real;
				}
			}
			else if(i == j + 1){
				if (abs(c[i].imag) < 0.000000000001){
					if (abs(c[i].real) < 0.000000000001){
						Ak[i][j] = 0.0;
					}
					else{
						Ak[i][j] = c[i].real;
					}
				}
				else{
					Ak[i][j] = c[i].real;
				}
			}
			else{
				Ak[i][j] = 0.0;
			}
		}
	}

	while(sum_diags < N){
		delta = (Ak[N-2][N-2]-Ak[N-1][N-1])/2;
		if(delta >= 0){
			sign = 1;
		}
		else{
			sign = -1;
		}

		mu = Ak[N-1][N-1] - ((sign*pow(Ak[N-1][N-2],2))/(abs(delta)+pow(pow(delta,2)+pow(Ak[N-1][N-2],2),0.5)));

		for(i = 0; i < N; i++){
			Ak[i][i] -= mu;
		}

		counter +=1;
		for( k = 1; k < N; k++){
			row = k;
			column = k - 1;
			for( i = 0; i<N;i++){
				for( j = 0; j<N;j++){
					G[k-1][i][j] = 0.0;
					if (i == j && i != row - 1 && j != column && i != row  && j != column + 1 ){
						G[k-1][i][j] = 1;
					}
					else if((i == row -1 && j == column) || (i == row && j == column + 1) ){
						G[k-1][i][j] = Ak[row - 1][column]/pow(pow(Ak[row - 1][column],2)+pow(Ak[row][column],2),0.5);
					}
					else if((i == row && j == column) || (i == row - 1 && j == column + 1) ){
						G[k-1][i][j] = (j-i)*Ak[row][column]/pow(pow(Ak[row - 1][column],2)+pow(Ak[row][column],2),0.5);
					}
				}
			}

			for( i =0;i<N;i++){
				for( j = 0;j < N; j++){
					temp[i][j] = 0.0;
					for( s = 0; s < N; s++){
						temp[i][j] += G[k-1][i][s]*Ak[s][j];

					}
				} 
			}
			copy(&temp[0][0],&temp[0][0]+N*N,&Ak[0][0]);

			for(i =0;i<N;i++){
				for(j = 0;j < N; j++){
					temp[i][j] = 0.0;
					for(int p = 0; p < N; p++){
						temp[i][j] += Ak[i][p]*G[k-1][j][p];
					}
				} 
			}
			copy(&temp[0][0],&temp[0][0]+N*N,&Ak[0][0]);

		}
		for(i = 0; i < N; i++){
			Ak[i][i] += mu;
			Identity[i][i] = 1;
		}
		if(Vecs){
			Qk.push_back(Identity);
			for( k = 1; k < N;k++){
				for( i =0;i<N;i++){
					for( j = 0;j < N; j++){
						temp2[i][j] = 0.0;
						for( s = 0; s < N; s++){
							temp2[i][j] += Qk[counter-1][i][s]*G[k-1][j][s];
						}
					}
				}

				for(i = 0; i < N; i ++){
					for(j = 0; j < N; j ++){
						Qk[counter-1][i][j] = temp2[i][j];
					}
				}
			}
		}
		for( k = 1; k < N;k++){
			sum_diags = 1;
			for( i =0;i<N;i++){
				for( j = 0;j < N; j++){
					if(i == j + 1){
						if(abs(Ak[i][j]) < 0.000000001){
							sum_diags += 1;
						}
					}
				}
			}
		}
	}
	if(Vecs){
		for(k = 1; k<counter;k++){
			for(i = 0; i < N; i++){
				for(j = 0; j < N; j++){
					Identity[i][j] = 0;
					for( s = 0; s < N; s++){
						Identity[i][j] += Qk[0][i][s]*Qk[k][s][j];
					}
				}
			}
			for(i = 0; i < N; i++){
				for(j = 0; j < N; j++){
					Qk[0][i][j] = Identity[i][j];
				}
			}
		}
	}
	for( i=0;i<N;i++)
	{		
		for( j=i+1;j<N;j++)
		{
			if(Ak[i][i]>Ak[j][j])
			{
				temporary  = Ak[i][i];
				Ak[i][i] = Ak[j][j];
				Ak[j][j] = temporary;
				if(Vecs){
					for(s = 0; s < N; s++){
						temporary = Qk[0][s][i];
						Qk[0][s][i] = Qk[0][s][j];
						Qk[0][s][j] = temporary;
					}
				}
			}
		}
	}

	for( i = 0;i <N;i++){
		Eigen.evalues[i] = Ak[i][i];
		if(Vecs){
			for(j = 0;j < N; j++){
				Eigen.evectors[i][j] = Qk[0][j][i];
			}
		}
	}
	return Eigen;
}

eigen band_eigenvalues_and_eigenvectors(comp l1,   comp m1,   comp n1, 
					 comp l1l1, comp m1m1, comp n1n1, 
					 comp l1m1, comp l1n1, comp m1n1, 
					 comp l2,   comp m2,   comp n2, 
					 comp l2l2, comp m2m2, comp n2n2, 
					 comp l2m2, comp l2n2, comp m2n2,
					 double es, double ep, double vssS, double vspS, double vppS, double vppP, 
					 double vssS2, double vspS2, double vppS2, double vppP2, bool get_vecs){
	// cout << "band_eigenvalues_and_eigenvectors" << endl;
	comp l1c = conj(l1);
	comp m1c = conj(m1);
	comp n1c = conj(n1);
	comp l1l1c = conj(l1l1);
	comp m1m1c = conj(m1m1);
	comp n1n1c = conj(n1n1);
	comp l1m1c = conj(l1m1);
	comp l1n1c = conj(l1n1);
	comp m1n1c = conj(m1n1);
	comp l2c = conj(l2);
	comp m2c = conj(m2);
	comp n2c = conj(n2);
	comp l2l2c = conj(l2l2);
	comp m2m2c = conj(m2m2);
	comp n2n2c = conj(n2n2);
	comp l2m2c = conj(l2m2);
	comp l2n2c = conj(l2n2);
	comp m2n2c = conj(m2n2);
	comp zero_comp = {};
	zero_comp.real = 0.0;
	zero_comp.imag = 0.0;
	comp es_comp = {};
	es_comp.real = es;
	es_comp.imag = 0.0; 
	comp ep_comp = {};
	ep_comp.real = ep;
	ep_comp.imag = 0.0; 
	comp comp_vppP = {};
	comp_vppP.real = vppP;
	comp_vppP.imag = 0.0;
	comp comp_vppP2 = {};
	comp_vppP2.real = vppP2;
	comp_vppP2.imag = 0.0;
	comp hamiltonian[N][N] = {{es_comp,				mult_comp(l1l1,vssS),	mult_comp(l2l2c,vssS2),		zero_comp,										zero_comp,										zero_comp,										mult_comp(l1,vspS),							mult_comp(m1,vspS),							mult_comp(n1,vspS),							mult_comp(l2c,-1*vspS2),						mult_comp(m2c,-1*vspS2),						mult_comp(n2c,-1*vspS2)						},
                    {mult_comp(l1l1c,vssS),	es_comp,				mult_comp(l1l1,vssS),		mult_comp(l1c,-1*vspS),							mult_comp(m1c,-1*vspS),							mult_comp(n1c,-1*vspS),							zero_comp,									zero_comp,									zero_comp,									mult_comp(l1,vspS),								mult_comp(m1,vspS),								mult_comp(n1,vspS)							},
        			{mult_comp(l2l2,vssS2),	mult_comp(l1l1c,vssS),	es_comp,					mult_comp(l2,vspS2),							mult_comp(m2,vspS2),							mult_comp(n2,vspS2),							mult_comp(l1c,-1*vspS),						mult_comp(m1c,-1*vspS),						mult_comp(n1c,-1*vspS),						zero_comp,										zero_comp,										zero_comp									},
        			{zero_comp,				mult_comp(l1,vspS),		mult_comp(l2c,-1*vspS2),	ep_comp,										zero_comp,										zero_comp,										add_comp(mult_comp(l1l1c,vppS-vppP),comp_vppP),	mult_comp(l1m1,vppS-vppP),					mult_comp(l1n1,vppS-vppP),					add_comp(mult_comp(l2l2c,vppS2-vppP2),comp_vppP2),	mult_comp(l2m2c,vppS2-vppP2),					mult_comp(l2n2c,vppS2-vppP2)				},
        			{zero_comp,				mult_comp(m1,vspS),		mult_comp(m2c,-1*vspS2),	zero_comp,										ep_comp,										zero_comp,										mult_comp(l1m1,vppS-vppP),					add_comp(mult_comp(m1m1,vppS-vppP),comp_vppP),	mult_comp(m1n1,vppS-vppP),					mult_comp(l2m2c,vppS2-vppP2),					add_comp(mult_comp(m2m2c,vppS2-vppP2),comp_vppP2),	mult_comp(m2n2c,vppS2-vppP2)				},
        			{zero_comp,				mult_comp(n1,vspS),		mult_comp(n2c,-1*vspS2),	zero_comp,										zero_comp,										ep_comp,										mult_comp(l1n1,vppS-vppP),					mult_comp(m1n1,vppS-vppP),					add_comp(mult_comp(n1n1,vppS-vppP),comp_vppP),	mult_comp(l2n2c,vppS2-vppP2),					mult_comp(m2n2c,vppS2-vppP2),					add_comp(mult_comp(n2n2c,vppS2-vppP2),comp_vppP2) },
        			{mult_comp(l1c,-1*vspS),zero_comp,				mult_comp(l1,vspS),			add_comp(mult_comp(l1l1c,vppS-vppP),comp_vppP),		mult_comp(l1m1c,vppS-vppP),						mult_comp(l1n1c,vppS-vppP),						ep_comp,									zero_comp,									zero_comp,									add_comp(mult_comp(l1l1,vppS-vppP),comp_vppP),		mult_comp(l1m1,vppS-vppP),						mult_comp(l1n1,vppS-vppP)					},
        			{mult_comp(m1c,-1*vspS),zero_comp,				mult_comp(m1,vspS),			mult_comp(l1m1c,vppS-vppP),						add_comp(mult_comp(m1m1c,vppS-vppP),comp_vppP),		mult_comp(m1n1c,vppS-vppP),						zero_comp,									ep_comp,									zero_comp,									mult_comp(l1m1,vppS-vppP),						add_comp(mult_comp(m1m1,vppS-vppP),comp_vppP),		mult_comp(m1n1,vppS-vppP)					},
        			{mult_comp(n1c,-1*vspS),zero_comp,				mult_comp(n1,vspS),			mult_comp(l1n1c,vppS-vppP),						mult_comp(m1n1c,vppS-vppP),						add_comp(mult_comp(n1n1c,vppS-vppP),comp_vppP),		zero_comp,									zero_comp,									ep_comp,									mult_comp(l1n1,vppS-vppP),						mult_comp(m1n1,vppS-vppP),						add_comp(mult_comp(n1n1,vppS-vppP),comp_vppP)	},
        			{mult_comp(l2,vspS2),	mult_comp(l1c,-1*vspS),	zero_comp,					add_comp(mult_comp(l2l2,vppS2-vppP2),comp_vppP2),	mult_comp(l2m2,vppS2-vppP2),					mult_comp(l2n2,vppS2-vppP2),					add_comp(mult_comp(l1l1c,vppS-vppP),comp_vppP),	mult_comp(l1m1c,vppS-vppP),					mult_comp(l1n1c,vppS-vppP),					ep_comp,										zero_comp,										zero_comp									},
        			{mult_comp(m2,vspS2),	mult_comp(m1c,-1*vspS),	zero_comp,					mult_comp(l2m2,vppS-vppP),						add_comp(mult_comp(m2m2,vppS2-vppP2),comp_vppP2),	mult_comp(m2n2,vppS2-vppP2),					mult_comp(l1m1c,vppS-vppP),					add_comp(mult_comp(m1m1c,vppS-vppP),comp_vppP),	mult_comp(m1n1c,vppS-vppP),					zero_comp,										ep_comp,										zero_comp									},
        			{mult_comp(m2,vspS2),	mult_comp(n1c,-1*vspS),	zero_comp,					mult_comp(l2n2,vppS-vppP),						mult_comp(m2n2,vppS2-vppP2),					add_comp(mult_comp(n2n2,vppS2-vppP2),comp_vppP2),	mult_comp(l1n1c,vppS-vppP),					mult_comp(m1n1c,vppS-vppP),					add_comp(mult_comp(n1n1c,vppS-vppP),comp_vppP),	zero_comp,										zero_comp,										ep_comp										}
        		};			
    eigen e = Householder_with_QR(hamiltonian, get_vecs);

	return e;
}


bandvecs_and_vals bandvectors_and_values(double params[10], double neighbors[4][3], 
	double next_neighbors[12][3], double kx, double ky, double kz, bool get_vecs){
	// cout << "bandvecs_and_values" << endl;
	bandvecs_and_vals b = {};
	int k = 0;
	int x = 0;
	int i = 0;
	int j = 0;

	phi g = phase(kx, ky, kz, neighbors, next_neighbors);
	eigen e = band_eigenvalues_and_eigenvectors(g.factors[0], g.factors[1],  g.factors[2],  g.factors[3], g.factors[4], g.factors[5],  g.factors[6], g.factors[7], 
		                                  g.factors[8], g.factors[9], g.factors[10], g.factors[11], g.factors[12], g.factors[13], g.factors[14], g.factors[15],
		                                  g.factors[16], g.factors[17], params[0], params[1], params[2], params[3],
		                                  params[4], params[5], params[6],params[7],params[8],params[9],get_vecs);

	for( i = 0; i < N; i++){
		b.bandvalues[i] = e.evalues[i];
		for(j = 0; j < N; j++){
			b.bandvectors[i][j] = e.evectors[i][j];
		}
	}
	

	return b;
}

double residual(double params[10], double neighbors[4][3], double next_neighbors[12][3], 
	double kx, double ky, double kz, double dft_bands[24], double dft_band_eigenvectors[24][N]){
	
	bandvecs_and_vals BAND_STRUCTURE = bandvectors_and_values(params, neighbors, next_neighbors, kx, ky, kz,true);

	double bands[N] = {0};
	double bandvecs[N][N] = {0};
	double diff_energies_squared = 0;
	double product = 0;
	double res = 0;
	double overlap = 0;

	
	for(int i = 0; i<N;i++){
		bands[i] = BAND_STRUCTURE.bandvalues[i];
		for(int j = 0; j<N;j++){
			bandvecs[i][j] = BAND_STRUCTURE.bandvectors[i][j];										
		}
	}


	for(int i = 0; i < 24; i++){
		for(int j = 0; j < N; j++){
			diff_energies_squared = pow(dft_bands[i]-bands[j],2);
			overlap = 0;
			for(int element = 0; element < N; element++){
				// cout << dft_band_eigenvectors[k][i][element] << "\t" << bandvecs[k][j][element] << endl;
				overlap += dft_band_eigenvectors[i][element]*bandvecs[j][element];
			}
			product = diff_energies_squared*pow(abs(overlap),0.5);
			res += product;
		}
	}

	return res;
}

void nelmin ( double fn (double x[10], double neighbors[4][3], double next_neighbors[12][3], 
	double dft_path[3], double dft_bands[24], double dft_band_eigenvectors[24][N]), int n, double start[10], double neighbors[4][3], double next_neighbors[12][3], 
	double dft_path[3], double dft_bands[24], double dft_band_eigenvectors[24][N], double xmin[], 
  double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
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
//      function fn ( x, f )
//      double fn
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
//    Input/output, double START[N].  On input, a starting point
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
  double *p;
  double *p2star;
  double *pbar;
  double *pstar;
  double rcoeff = 1.0;
  double rq;
  double x;
  double *y;
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

  p = new double[n*(n+1)];
  pstar = new double[n];
  p2star = new double[n];
  pbar = new double[n];
  y = new double[n+1];

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
      p[i+n*n] = start[i];
    }
    y[n] = fn ( start , neighbors, next_neighbors, dft_path, dft_bands, dft_band_eigenvectors);
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      y[j] = fn ( start , neighbors, next_neighbors, dft_path, dft_bands, dft_band_eigenvectors);
      *icount = *icount + 1;
      start[j] = x;
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
      ystar = fn ( pstar , neighbors, next_neighbors, dft_path, dft_bands, dft_band_eigenvectors);
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
        y2star = fn ( p2star , neighbors, next_neighbors, dft_path, dft_bands, dft_band_eigenvectors);
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
          y2star = fn ( p2star , neighbors, next_neighbors, dft_path, dft_bands, dft_band_eigenvectors);
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
              y[j] = fn ( xmin , neighbors, next_neighbors, dft_path, dft_bands, dft_band_eigenvectors);
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
          y2star = fn ( p2star , neighbors, next_neighbors, dft_path, dft_bands, dft_band_eigenvectors);
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
      z = fn ( xmin , neighbors, next_neighbors, dft_path, dft_bands, dft_band_eigenvectors);
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
        break;
      }
      xmin[i] = xmin[i] - del - del;
      z = fn ( xmin , neighbors, next_neighbors, dft_path, dft_bands, dft_band_eigenvectors);
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
//  Restart the procedure.
//
    for ( i = 0; i < n; i++ )
    {
      start[i] = xmin[i];
    }
    del = eps;
    *numres = *numres + 1;
  }
  delete [] p;
  delete [] pstar;
  delete [] p2star;
  delete [] pbar;
  delete [] y;

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

int main(){

	////////////////////////////////////////////////////////////////////////////////////////
	//
	//  Get data from SCAN DFT
	//
	////////////////////////////////////////////////////////////////////////////////////////

	ifstream PROCAR;
	PROCAR.open("PROCAR");
	double kx = 0, ky = 0, kz = 0;
	double tmp_dft_path[900][3] = {0};
	double dft_energy = 0;
	double tmp_dft_bands[900][24] = {0};
	double s = 0, py = 0, pz = 0, px = 0, dxy = 0, dyz = 0, dz2 = 0, dxz = 0, dx2 = 0;
	double tmp_dft_band_eigenvectors[900][24][N] = {0};

	string throwaway = "";

	PROCAR >> throwaway >> throwaway >> throwaway;
	PROCAR >> throwaway >> throwaway >> throwaway;
	PROCAR >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway;
	PROCAR >> throwaway >> throwaway >> throwaway >> kx >> ky >> kz >> throwaway >> throwaway >> throwaway;
	tmp_dft_path[0][0] = kx;
	tmp_dft_path[0][1] = ky;
	tmp_dft_path[0][2] = kz;

	int counter = 4;
	int element = 0;
	int e_element = 0;
	int ev_element = 0;
	while(PROCAR){	

		if((counter-4)%265==0){
			//k-points
			if(counter-4 != 0){
				kx = 0, ky = 0, kz = 0;
				PROCAR >> throwaway >> throwaway >> throwaway >> kx >> ky >> kz >> throwaway >> throwaway >> throwaway;
				element += 1;
				tmp_dft_path[element][0] = kx;
				tmp_dft_path[element][1] = ky;
				tmp_dft_path[element][2] = kz;
				e_element = 0;
				ev_element = 0;
			}
		}
		if(((counter-4)%265)%(8+3)==1){
			dft_energy = 0;
			PROCAR >> throwaway >> throwaway >> throwaway >> throwaway >> dft_energy >> throwaway >> throwaway >> throwaway;
			tmp_dft_bands[element][e_element] = dft_energy;
			e_element += 1;
		}
		if((((counter-4)%265)%(8+3)==0) && (counter-4)%265!=0){
			s = 0, py = 0, pz = 0, px = 0, dxy = 0, dyz = 0, dz2 = 0, dxz = 0, dx2 = 0;
			PROCAR >> throwaway >> s >> py >> pz >> px >> dxy >> dyz >> dz2 >> dxz >> dx2>>throwaway;
			tmp_dft_band_eigenvectors[element][ev_element][0] = s;
			tmp_dft_band_eigenvectors[element][ev_element][1] = py;
			tmp_dft_band_eigenvectors[element][ev_element][2] = pz;
			tmp_dft_band_eigenvectors[element][ev_element][3] = px;
			tmp_dft_band_eigenvectors[element][ev_element][4] = dxy;
			tmp_dft_band_eigenvectors[element][ev_element][5] = dyz;
			tmp_dft_band_eigenvectors[element][ev_element][6] = dz2;
			tmp_dft_band_eigenvectors[element][ev_element][7] = dxz;
			tmp_dft_band_eigenvectors[element][ev_element][8] = dx2;
			tmp_dft_band_eigenvectors[element][ev_element][9] = 0;
			tmp_dft_band_eigenvectors[element][ev_element][10] = 0;
			tmp_dft_band_eigenvectors[element][ev_element][11] = 0;
			ev_element += 1;
		}
		if(!((counter-4)%265==0) && !(((counter-4)%265)%(8+3)==1) && !((((counter-4)%265)%(8+3)==0) && (counter-4)%265!=0)){
			PROCAR >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway;
		}
		counter+=1;
	}

	double dft_path[P][3] = {0};
	double dft_band_eigenvectors[P][24][N] = {0};
	double dft_bands[P][24] = {0};
	double kx_ = 0, ky_ = 0, kz_ = 0;
	double epsilon = 0.000000001;
	double normalization = 0;
	counter = 0;
	for(int i = 1; i < P; i++){
		kx_ = tmp_dft_path[i-1][0];
		ky_ = tmp_dft_path[i-1][1];	
		kz_ = tmp_dft_path[i-1][2];
		kx = tmp_dft_path[i][0];
		ky = tmp_dft_path[i][1];
		kz = tmp_dft_path[i][2];
		if(!(abs(kx - kx_) < epsilon && abs(ky - ky_) < epsilon  && abs(kz - kz_) < epsilon)){
			dft_path[counter][0] = kx_;
			dft_path[counter][1] = ky_;
			dft_path[counter][2] = kz_;
			for(int j = 0; j < 24; j++){
				dft_bands[counter][j] = tmp_dft_bands[i-1][j];
				for(int k = 0; k < N; k++){
					normalization = 0;
					for(int element = 0; element < N; element++){
						normalization += pow(tmp_dft_band_eigenvectors[i-1][j][element],2);
					}
					dft_band_eigenvectors[counter][j][k] = tmp_dft_band_eigenvectors[i-1][j][k]/pow(normalization,0.5);
				}
			}

			counter += 1;
		}
	}

	if(P == 894){
	dft_path[893][0] = 0;
	dft_path[893][1] = 0;
	dft_path[893][2] = 0;		
	}


	double ticks[P];
	ticks[0] = 0;
	double tick = 0;
	double max_mom_x = 0;
	double max_mom_y = 0;
	double max_mom_z = 0;

	for(int i = 1; i < P; i++){
		if(dft_path[i-1][0]>max_mom_x){
			max_mom_x = dft_path[i-1][0];

		}
		if(dft_path[i-1][1]>max_mom_y){
			max_mom_y = dft_path[i-1][1];
		}
		if(dft_path[i-1][2]>max_mom_z){
			max_mom_z = dft_path[i-1][2];

		}
		if(pow(pow(dft_path[i][0]-dft_path[i-1][0],2)+pow(dft_path[i][1]-dft_path[i-1][1],2)+pow(dft_path[i][2]-dft_path[i-1][2],2),0.5) < 1.0){
			tick += pow(pow(dft_path[i][0]-dft_path[i-1][0],2)+pow(dft_path[i][1]-dft_path[i-1][1],2)+pow(dft_path[i][2]-dft_path[i-1][2],2),0.5);
		}
		else{
			tick += 0.0000000000001;
		}

		ticks[i] = tick;
	}

	max_mom_x *= M_PI/a;
	max_mom_y *= M_PI/a;
	max_mom_z *= M_PI/a;


	///////////////////////////////////////////////////////////////////////////////////////
	//// Plot DFT Bandstructure
	///////////////////////////////////////////////////////////////////////////////////////
	// ofstream bandvecs_and_values_data_dft;
	// bandvecs_and_values_data_dft.open("BAND_STRUCTURE_DFT.dat");
	// counter = 0;

	// for(int j = 0; j < 24; j++){
	// 	for(int k = 0; k < P; k++){
	// 		bandvecs_and_values_data_dft << ticks[k] << "\t" << dft_bands[k][j] << endl;
	// 	}
	// 	bandvecs_and_values_data_dft << endl;
	// }
	// bandvecs_and_values_data_dft.close();
	// string command = "gnuplot -e \"set yrange [-10:10]; plot \\\"BAND_STRUCTURE_DFT.dat\\\" u 1:2 w l; pause mouse \"";
	// const char* com = command.c_str();
	// system(com);

	////////////////////////////////////////////////////////////////////////////////////////
	//**************************************************************************************
	////////////////////////////////////////////////////////////////////////////////////////

	double neighbors[4][3] = {
	{0.25*a,-0.25*a,0.25*a},
	{0.25*a,0.25*a,-0.25*a},
	{-0.25*a,0.25*a,0.25*a},
	{-0.25*a,-0.25*a,-0.25*a}};

	double next_neighbors[12][3] = {
	{0.5*a,0.0,-0.5*a},
	{0.0,0.5*a,-0.5*a},
	{-0.5*a,0,-0.5*a},
	{0.0,-0.5*a,-0.5*a},
	{0.5*a,0.5*a,0.0},
	{-0.5*a,0.5*a,0.0},
	{-0.5*a,-0.5*a,0.0},
	{0.5*a,-0.5*a,0.0},
	{0.5*a,0,0.5*a},
	{0.0,0.5*a,0.5*a},
	{-0.5*a,0,0.5*a},
	{0.0,-0.5*a,0.5*a}};


	// double params[10] = {-4.03,2.9,-2.035, 1.47, -1.45, 4.1825, -1.425657743, 1.212435565, 2.045116134,-1.204159458};
	// double params[10] = {1,-3.40277,-1.77416,-0.84730,1.53114,-0.716431,-0.0891011,-0.105639,0.388664,-0.126518};
	// double params[10] ={1.39458,	2.59735,	-4.0652,	-3.65098,	0.202101,	-4.21768,	-4.30094,	-2.95345,	-0.385795,	3.19677};
	double params[10] = {0.567518,	3.28803,	12.9791,	-0.821661,	2.32315,	-11.4802,	-7.48998,	0.730091,	-0.500956,	3.66033	};
	// double new_params[10];
	// double min_res = 1000000000;
	// double re = 0;
	// for(int i =0; i <1000; i++){
	// 	for(int j = 0; j<10;j++){
	// 		new_params[j] = params[j] += fRand(-1,1);
	// 	}
	// 	cout << i << endl;
	// 	re = residual(new_params,neighbors,next_neighbors,dft_path,dft_bands,dft_band_eigenvectors);		
	// 	if(re < min_res){
	// 		min_res = re;
	// 		cout << re << endl;
	// 		cout << "{";
	// 		for(int l = 0; l<10;l++){
	// 			cout << params[l] << ",\t";
	// 			params[l] = new_params[l];
	// 		}
	// 		cout << "}"<<endl;
	// 	}
	// }
	// cout << residual(params,neighbors,next_neighbors,dft_path,dft_bands,dft_band_eigenvectors) << endl;

	// int icount;
	// int ifault;
	// int kcount;
	// int konvge;
	// int n;
	// int numres;
	// double reqmin;
	// double *start;
	// double *step;
	// double *xmin;
	// double ynewlo;
	// n = 10;
	// start = new double[n];	
	// step = new double[n];
	// xmin = new double[n];
	// konvge = 10;
	// kcount = 100;
	// cout << "\n";
	// cout << "  Starting point X:\n";
	// cout << "\n";
	// for (int  i = 0; i < n; i++ )
	// {	
	// 	step[i] = 1;
	// 	start[i] = params[i];
	// 	cout << "  " << setw(14) << start[i] << "\n";
	// }
	// ynewlo = residual(start,neighbors,next_neighbors,dft_path,dft_bands,dft_band_eigenvectors);
	// cout << "\n";
	// cout << "  F(X) = " << ynewlo << "\n";
	// nelmin ( residual, n, start, neighbors, next_neighbors, dft_path, dft_bands, dft_band_eigenvectors, xmin, &ynewlo, reqmin, step,
	// konvge, kcount, &icount, &numres, &ifault );
	// cout << "\n";
	// cout << "  Return code IFAULT = " << ifault << "\n";
	// cout << "\n";
	// cout << "  Estimate of minimizing value X*:\n";
	// cout << "\n";
	// for (int i = 0; i < n; i++ )
	// {
	// cout << "  " << setw(14) << xmin[i] << "\n";
	// }

	// cout << "\n";
	// cout << "  F(X*) = " << ynewlo << "\n";

	// cout << "\n";
	// cout << "  Number of iterations = " << icount << "\n";
	// cout << "  Number of restarts =   " << numres << "\n";

	// delete [] start;
	// delete [] step;
	// delete [] xmin;
	// ///////////////////////////////////////////////////////////////////////////////////////
	// //// Plot TB Bandstructure
	// ///////////////////////////////////////////////////////////////////////////////////////

	kx = 0.5;
	ky = 0.5;
	kz = 0.5;

	proj_bandvecs_and_vals BANDS = projection_gauss_10_pt_quad(bandvectors_and_values,s_orb, py_orb, pz_orb,px_orb,dxy_orb,dyz_orb,dz_2_orb,dxz_orb,dx_2_y_2_orb,fyz_2_orb,fz_3_orb,fxz_2_orb,
	params, neighbors, next_neighbors, -a, a, -a, a, -a, a, -M_PI/a,M_PI/a,-M_PI/a,M_PI/a,-M_PI/a,M_PI/a,kx,ky,kz);

	for(int i = 0; i < N; i++){
		cout << BANDS.bandvalues[i] << ": " << endl << endl;
		for(int j = 0; j < N; j++){
			cout << BANDS.projected_vec[i][j].real << "\t" << BANDS.projected_vec[i][j].imag << endl;
		}
		cout << endl;
	}


	// double bands[P][N] = {0};

	// for(int i = 0; i < P; i++){
	// 	for(int j = 0; j < N; j++){
	// 		bands[i][j] = BAND_STRUCTURE.bandvalues[i][j];
	// 	}
	// }
	// ofstream bandvecs_and_values_data_tb;
	// bandvecs_and_values_data_tb.open("BAND_STRUCTURE_TB.dat");

	// for(int j = 0; j < N; j++){
	// 	for(int k = 0; k < P; k++){
	// 		bandvecs_and_values_data_tb << ticks[k] << "\t" << bands[k][j] << endl;
	// 	}
	// 	bandvecs_and_values_data_tb << endl;
	// }
	// bandvecs_and_values_data_tb.close();
	// string command2= "gnuplot -e \"set xrange [0:1]; set yrange [-20:20]; plot \\\"BAND_STRUCTURE_TB.dat\\\" u 1:2 w l; pause mouse \"";
	// const char* com2 = command2.c_str();
	// system(com2);

	////////////////////////////////////////////////////////////////////////////////////////
	//**************************************************************************************
	////////////////////////////////////////////////////////////////////////////////////////

	return 0;
}