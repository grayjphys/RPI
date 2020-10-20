
# include <cstdlib>
# include <iostream>
# include <string>
# include <iomanip>
# include <cmath>
# include <fstream>
# include <vector>
using namespace std;
const int N = 12;
const int P = 1;
// P = 894
struct comp{double real = 0; double imag = 0; };
struct dot{double p1[4] = {0}; double p2[12] = {0}; };
struct phi{comp factors[18] = {}; };
struct eigen{double evalues[N] = {0}; double evectors[N][N] = {}; };
struct band_struct{double bandvalues[P][N] = {0}; double bandvectors[P][N][N] = {0}; };

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

phi phase(double k[3], double neighbors[4][3], double next_neighbors[12][3]){
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
		nnp[i].real = cos(k[0]*neighbors[i][0]+k[1]*neighbors[i][1]+k[2]*neighbors[i][2]);
		nnp[i].imag = sin(k[0]*neighbors[i][0]+k[1]*neighbors[i][1]+k[2]*neighbors[i][2]);
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
		nnnp[i].real = cos(k[0]*next_neighbors[i][0]+k[1]*next_neighbors[i][1]+k[2]*next_neighbors[i][2]);
		nnnp[i].imag = sin(k[0]*next_neighbors[i][0]+k[1]*next_neighbors[i][1]+k[2]*next_neighbors[i][2]);
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

eigen Householder_with_QR(comp hamiltonian[N][N]){
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

	// double sign = 0;
	// double delta = (b[N-2].real-b[N-1].real)/2;
	// double mu = b[N-1].real -  sign*pow(c[N-2].real,2)/(abs(delta)+pow(pow(delta,2)+pow(c[N-2].real,2),0.5));
	for( i = 0; i < N; i++){
		for( j = 0; j < N; j++){
			if(i == j){
				// sign = (delta > 0) - (delta < 0);
				Ak[i][j] = b[i].real ;
			}
			else if(j == i + 1){
				if (abs(c[j].imag) < 0.000000001){
					if (abs(c[j].real) < 0.000000001){
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
				if (abs(c[i].imag) < 0.000000001){
					if (abs(c[i].real) < 0.000000001){
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
			for( i =0;i<N;i++){
				for( j = 0;j < N; j++){
					Ak[i][j] = temp[i][j];
				}
			}

		}
		Qk.push_back(Identity);
		for( k = 1; k < N;k++){
			for( i =0;i<N;i++){
				for( j = 0;j < N; j++){
					temp[i][j] = 0.0;
					temp2[i][j] = 0.0;
					for( s = 0; s < N; s++){
						temp[i][j] += Ak[i][s]*G[k-1][j][s];
						temp2[i][j] += Qk[counter-1][i][s]*G[k-1][j][s];
					}
				}
			}
			sum_diags = 1;

			for(i = 0; i < N; i ++){
				for(j = 0; j < N; j ++){
					Ak[i][j] = temp[i][j];
					Qk[counter-1][i][j] = temp2[i][j];
				}
			}
			for( i =0;i<N;i++){
				for( j = 0;j < N; j++){
					if(i == j + 1){
						if(abs(Ak[i][j]) < 0.000000001){
							sum_diags += 1;
						};
					}
				}
			}
		}
	}

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

	for( i=0;i<N;i++)
	{		
		for( j=i+1;j<N;j++)
		{
			if(Ak[i][i]>Ak[j][j])
			{
				temporary  = Ak[i][i];
				Ak[i][i] = Ak[j][j];
				Ak[j][j] = temporary;
				for(s = 0; s < N; s++){
					temporary = Qk[0][s][i];
					Qk[0][s][i] = Qk[0][s][j];
					Qk[0][s][j] = temporary;
				}
			}
		}
	}

	for( i = 0;i <N;i++){
		Eigen.evalues[i] = Ak[i][i];
		for(j = 0;j < N; j++){
			Eigen.evectors[i][j] = Qk[0][j][i];
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
					 double vssS2, double vspS2, double vppS2, double vppP2){
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
    eigen e = Householder_with_QR(hamiltonian);

	return e;
}

band_struct band_structure(double params[10], double neighbors[4][3], 
	double next_neighbors[12][3], double path[P][3]){
	// cout << "band_structure" << endl;
	band_struct b = {};
	int k = 0;
	int x = 0;
	int i = 0;
	int j = 0;

	for(k = 0; k < P; k++){
		double vec[3];
		for(x = 0; x < 3; x++){
			vec[x] = path[k][x];
		}
		phi g = phase(vec, neighbors, next_neighbors);
		eigen e = band_eigenvalues_and_eigenvectors(g.factors[0], g.factors[1],  g.factors[2],  g.factors[3], g.factors[4], g.factors[5],  g.factors[6], g.factors[7], 
			                                  g.factors[8], g.factors[9], g.factors[10], g.factors[11], g.factors[12], g.factors[13], g.factors[14], g.factors[15],
			                                  g.factors[16], g.factors[17], params[0], params[1], params[2], params[3],
			                                  params[4], params[5], params[6],params[7],params[8],params[9]);
		for( i = 0; i < N; i++){
			b.bandvalues[k][i] = e.evalues[i];
			for(j = 0; j < N; j++){
				b.bandvectors[k][i][j] = e.evectors[i][j];
			}
		}
	}

	return b;
}

double OVERLAP(double dft_band_eigenvector[N], double tb_band_eigenvector[N]){
	double over = 0;
	for(int element = 0; element < N; element++){
		over += dft_band_eigenvector[element]*tb_band_eigenvector[element];
	}

	return over;
}

double residual(double params[10], double neighbors[4][3], double next_neighbors[12][3], 
	double dft_path[P][3], double dft_bands[P][24], double dft_band_eigenvectors[P][24][N]){
	
	band_struct BAND_STRUCTURE = band_structure(params, neighbors, next_neighbors, dft_path);

	double bands[P][N] = {0};
	double normal = 0;
	double bandvecs[P][N][N] = {0};

	for(int k = 0; k<P;k++){
		for(int i = 0; i<N;i++){
			bands[k][i] = BAND_STRUCTURE.bandvalues[k][i];
			normal = 0;
			for(int j = 0; j<24;j++){
				if(j < N){
					bandvecs[k][i][j] = BAND_STRUCTURE.bandvectors[k][i][j];						
				}				
			}
		}
	}

	double diff_energies_squared;
	double product;
	double res = 0;
	double overlap = 0;
	for(int k = 0; k < P; k++){
		for(int i = 0; i < 24; i++){
			for(int j = 0; j < N; j++){
				diff_energies_squared = pow(dft_bands[k][i]-bands[k][j],2);
				product = diff_energies_squared*pow(abs(OVERLAP(dft_band_eigenvectors[k][i],bandvecs[k][j])),0.5);
				res += product;
			}
		}
	}

	return res;
}

