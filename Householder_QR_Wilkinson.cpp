#include <iostream>
#include <iomanip>
#include <cmath>
#include <vector>

using namespace std;

const int N = 4;
// const double 0.000000001 = (-double(1) + double(7)/double(3) - double(4)/double(3));
struct eigen{double evalues[N] = {0}; double evectors[N][N] = {}; };
struct comp{double real; double imag; };

comp conj(comp number);

comp mult_comp(comp c_number, double number);

comp mult_two_comps(comp number1, comp number2);

comp add_comp(comp number1, comp number2);

double norm(comp number);

int main(){

	int i = 0;
	int j = 0;

	comp A[N][N];
	double Ak[N][N]= {0.0};
	vector<vector<vector<double> > > Qk;
	vector<vector<double> > Identity(N,vector<double>(N));
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			A[i][j].real = 2*(static_cast <double> (rand()) / static_cast <double> (RAND_MAX)) - 1; 
			A[i][j].imag = 2*(static_cast <double> (rand()) / static_cast <double> (RAND_MAX)) - 1;
			if(i == j){
				Identity[i][j] = 1.0;
			}
		}
	}

	comp element;
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			element = conj(A[j][i]);
			A[i][j].real = A[i][j].real + element.real;
			A[i][j].imag = A[i][j].imag + element.imag;
		}
	}

	// cout << "A: Hermitian Matrix" << endl;
	// for(i = 0; i < N; i++){
	// 	for(j = 0; j < N; j++){
	// 		cout << "(";
	// 		if(A[i][j].real >= 0){
	// 			cout << " ";
	// 		}
	// 		cout << fixed <<  setprecision(5) << A[i][j].real;
	// 		if(A[i][j].imag >= 0){
	// 			cout << "+";
	// 		} 
	// 		cout << A[i][j].imag << "i)\t";
	// 	}
	// 	cout << endl;
	// }
	// cout << endl;

	int k = 0;
	comp b[N];
	comp c[N-1];
	double gamma = 0;
	double phi = 0;
	double alpha = 0;
	comp rho;
	comp beta;
	comp u[N];
	comp t[N];
	comp r[N];
	for(k = 0; k<N-2;k++){
		b[k] = A[k][k];
		gamma = 0;
		for(int i = k + 1; i < N; i++){
			if(norm(A[i][k]) > gamma){
				gamma = norm(A[i][k]);
			}
		phi = 0;
		if(gamma < 0.00000000000001){
			c[k].real = 0;
			c[k].imag = 0;
		}
		else{
			alpha = 0;
			if(norm(A[k+1][k]) < 0.00000000000001){
				beta.real = 1;
				beta.imag = 0;
			}
			else{
				beta = mult_comp(A[k+1][k],1/norm(A[k+1][k]));
			}
			
			for(int i = k + 1; i < N; i++){
				u[i] = mult_comp(A[i][k],1/gamma);
				alpha += pow(norm(u[i]),2);
			}

			alpha = pow(alpha, 0.5);
			c[k].real = -1*alpha;
			phi = 1/(alpha*(alpha+norm(u[k+1])));	
			u[k+1] = add_comp(mult_comp(beta, alpha),u[k+1]);
			rho.real = 0;
			rho.imag = 0;

			for(int s = k + 1; s < N; s++){
				t[s].real = 0;
				t[s].imag = 0;
				for(int j = k+1; j<s; j++){
					t[s] = add_comp(t[s],mult_two_comps(u[j],A[s][j]));
				}
				for(int j = s+1; j<N; j++){
					t[s] = add_comp(t[s],mult_two_comps(u[j],A[j][s]));
				}
				t[s] = mult_comp(t[s],phi);
				rho = add_comp(rho,mult_two_comps(conj(u[s]),t[s]));
			}

			for(int s = k + 1; s < N; s++){
				r[s] = mult_two_comps(mult_comp(rho, phi/2),u[s]);
			}
			for(int i = k + 1; i <N; i++){
				for(int j = k + 1; j < i; j++){
					A[i][j] = add_comp(A[i][j],add_comp(add_comp(mult_two_comps(mult_comp(t[i],-1),conj(u[j])),mult_two_comps(mult_comp(conj(t[j]),-1),u[i])),
					add_comp(mult_two_comps(u[i],conj(r[j])),mult_two_comps(r[i],conj(u[j])))));
				}
			}
		}
		}
	}
	c[N-2] = A[N-1][N-2];
	b[N-2] = A[N-2][N-2];
	b[N-1] = A[N-1][N-1];

	for( i = 0; i < N; i++){
		for( j = 0; j < N; j++){
			if(i == j){
				Ak[i][j] = b[i].real;
			}
			else if(j == i + 1){
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
			else if(i == j + 1){
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
			else{
				Ak[i][j] = 0.0;
			}
		}
	}


	// cout << "A0: Real Symmetric Tridiagonal Matrix" << endl;
	// for(i = 0; i < N; i++){
	// 	for(j = 0; j < N; j++){
	// 		if(Ak[i][j] >= 0){
	// 			cout << " ";
	// 		}
	// 		cout << fixed <<  setprecision(5) << Ak[i][j] << " ";
	// 	}
	// 	cout << endl;
	// }
	// cout << endl;

	

	double G[N-1][N][N] = {0.0};
	int row = 0;
	int column = 0;
	double temp[N][N] = {0.0};
	double temp2[N][N] = {0.0};
	int sum_diags = 0;
	int counter = 0;
	double sign = 0;
	double delta = 0;
	double mu = 0;
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
		for(k = 1; k < N; k++){
			row = k;
			column = k - 1;
			for(i = 0; i<N;i++){
				for(j = 0; j<N;j++){
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

			for(i =0;i<N;i++){
				for(j = 0;j < N; j++){
					temp[i][j] = 0.0;
					for(int l = 0; l < N; l++){
						temp[i][j] += G[k-1][i][l]*Ak[l][j];						
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
		}
		Qk.push_back(Identity);
		for(k = 1; k < N; k++){
			for(i = 0;i < N; i++){
				for(j = 0;j < N; j++){					
					temp2[i][j] = 0.0;
					for(int l = 0; l < N; l++){
						temp2[i][j] += Qk[counter-1][i][l]*G[k-1][j][l];
					}
				}
			}


			sum_diags = 1;

			for(i = 0; i < N; i ++){
				for(j = 0; j < N; j ++){
					Qk[counter-1][i][j] = temp2[i][j];
				}
			}
			for(i =0;i<N;i++){
				for(j = 0;j < N; j++){
					if(i == j + 1){
						if(abs(Ak[i][j]) < 0.0000000001){
							sum_diags += 1;
						}
					}
				}
			}
		}

		// cout << "A" << counter << endl;
		// for(i = 0; i < N; i++){
		// 	for(j = 0; j < N; j++){
		// 		if(Ak[i][j] >= 0){
		// 			cout << " ";
		// 		}
		// 		cout << fixed <<  setprecision(5) << Ak[i][j] << " ";
		// 	}
		// 	cout << endl;
		// }
		// cout << endl;
	} 

	cout << "Eigenvalues" << endl;
	for(i = 0; i < N; i++){
		if(Ak[i][i] >= 0){
				cout << " ";
		}
		cout << fixed <<  setprecision(16) << Ak[i][i] << endl;
	}
	cout << endl;

	for(k = 1; k<counter;k++){
		for(i = 0; i < N; i++){
			for(j = 0; j < N; j++){
				Identity[i][j] = 0;
				for(int l = 0; l < N; l++){
					Identity[i][j] += Qk[0][i][l]*Qk[k][l][j];
				}
			}
		}
		for(i = 0; i < N; i++){
			for(j = 0; j < N; j++){
				Qk[0][i][j] = Identity[i][j];
			}
		}

	}
	cout << "Eigenvectors" << endl;
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			if(Qk[0][j][i] >= 0){
					cout << " ";
			}
			cout << fixed <<  setprecision(16) << Qk[0][j][i] << endl;
		}
		cout << endl;
	}
	cout << endl;

	double temporary = 0;
	for(int i=0;i<N;i++)
	{		
		for(int j=i+1;j<N;j++)
		{
			if(Ak[i][i]>Ak[j][j])
			{
				temporary  = Ak[i][i];
				Ak[i][i] = Ak[j][j];
				Ak[j][j] = temporary;
				for(int s = 0; s < N; s++){
					temporary = Qk[0][s][i];
					Qk[0][s][i] = Qk[0][s][j];
					Qk[0][s][j] = temporary;
				}
			}
		}
	}

	eigen Eigen = {};
	for( i = 0;i <N;i++){
		Eigen.evalues[i] = Ak[i][i];
		for(j = 0;j < N; j++){
			Eigen.evectors[i][j] = Qk[0][j][i];
		}
	}

	cout << "Eigenvectors" << endl;
	for(i = 0; i < N; i++){
		for(j = 0; j < N; j++){
			if(Eigen.evectors[i][j] >= 0){
					cout << " ";
			}
			cout << fixed <<  setprecision(16) << Eigen.evectors[i][j] << endl;
		}
		cout << endl;
	}
	cout << endl;

	return 0;
}

comp conj(comp number){
	comp conjugate;
	conjugate.real = number.real;
	conjugate.imag = -1.0*number.imag;
	return conjugate;
}

comp mult_comp(comp c_number, double number){
	comp mult;
	mult.real = c_number.real*number;
	mult.imag = c_number.imag*number;
	return mult;
}

comp mult_two_comps(comp number1, comp number2){
	comp new_number;
	new_number.real = number1.real*number2.real - number1.imag*number2.imag;
	new_number.imag = number1.real*number2.imag + number1.imag*number2.real;
	return new_number;
}
comp add_comp(comp number1, comp number2){
	comp add;
	add.real = number1.real + number2.real;
	add.imag = number1.imag + number2.imag;
	return add;
}

double norm(comp number){
	return pow(pow(number.real,2) + pow(number.imag,2),0.5);
}