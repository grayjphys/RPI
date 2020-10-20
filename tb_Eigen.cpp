// Re-do the DFT calculations including NBANDS=50. It looks like the TB conduction band is very high. <- This could be the problem with the sp3 model. Try the sp3d5 model.

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

# define my_sizeof(type) ((char *)(&type+1)-(char*)(&type))

const double a = 5.431;
const int N = 8;
const int P = 396;
// P = 892

struct proj_bandvecs{ complex<double> projected_vec[N][N] = {}; };
struct eigen{double evalues[N] = {0}; complex<double> evectors[N][N] = {}; };

eigen band_eigenvalues_and_eigenvectors(double params[6], double neighbors[4][3], double kx, double ky, double kz){

	eigen e = {};
	complex<double> g[4] = {(cos(0.5*M_PI*kx)*cos(0.5*M_PI*ky)*cos(0.5*M_PI*kz),-sin(0.5*M_PI*kx)*sin(0.5*M_PI*ky)*sin(0.5*M_PI*kz)),
							(-cos(0.5*M_PI*kx)*sin(0.5*M_PI*ky)*sin(0.5*M_PI*kz),sin(0.5*M_PI*kx)*cos(0.5*M_PI*ky)*cos(0.5*M_PI*kz)),
							(-sin(0.5*M_PI*kx)*cos(0.5*M_PI*ky)*sin(0.5*M_PI*kz),cos(0.5*M_PI*kx)*sin(0.5*M_PI*ky)*cos(0.5*M_PI*kz)),
							(-sin(0.5*M_PI*kx)*sin(0.5*M_PI*ky)*cos(0.5*M_PI*kz),cos(0.5*M_PI*kx)*cos(0.5*M_PI*ky)*sin(0.5*M_PI*kz))};
	complex<double> gc[4] = {conj(g[0]),conj(g[1]),conj(g[2]),conj(g[3])};
	double es = params[0];
	double ep = params[1];
	double vss = params[2];
	double vsp = params[3];
	double vxx = params[4];
	double vxy = params[5];

	Matrix<complex<double>, N, N> hamiltonian;
	hamiltonian << es,	      g[0]*vss,	 0,		     0,			 0,			 g[1]*vsp, g[2]*vsp, g[3]* vsp,
                   gc[0]*vss, es,		 -gc[1]*vsp, -gc[2]*vsp, -gc[3]*vsp, 0,		   0,		 0,
    		       0,	      -g[1]*vsp, ep,		 0,			 0,			 g[0]*vxx, g[3]*vxy, g[2]* vxy,
    			   0,		  -g[2]*vsp, 0,	         ep,		 0,			 g[3]*vxy, g[0]*vxx, g[1]* vxy,
    			   0,		  -g[3]*vsp, 0,	         0,	         ep,		 g[2]*vxy, g[1]*vxy, g[0]* vxx,
    			   gc[1]*vsp, 0,		 gc[0]*vxx,	 gc[3]*vxy,  gc[2]*vxy,	 ep,	   0,		 0,
    			   gc[2]*vsp, 0,		 gc[3]*vxy,	 gc[0]*vxx,  gc[1]*vxy,	 0,		   ep,	     0,
    			   gc[3]*vsp, 0,		 gc[2]*vxy,  gc[1]*vxy,  gc[0]*vxx,	 0,		   0,		 ep;             

    SelfAdjointEigenSolver<Matrix<complex<double>, N, N> > s(hamiltonian);
    
    for(int i = 0; i <N; i++){
    	e.evalues[i] = s.eigenvalues()(i);
	    for(int j = 0; j <N; j++){
	    	e.evectors[i][j] = s.eigenvectors()(i,j);
	    }
    }
	return e;
}


proj_bandvecs projection_gauss_10_pt_quad(complex<double> sh[N][1331], double params[6], double neighbors[4][3], double akx, double bkx, double aky, double bky, double akz, double bkz){
	proj_bandvecs TB_Proj = {};

	static double points[11]  = {-0.97822865814605699280393800112285739077142240891978,
                          -0.88706259976809529907515776930392726663167575122531,
                          -0.73015200557404932409341625203115345804964306202613,
                          -0.51909612920681181592572566945860955448022711511993,
                          -0.26954315595234497233153198540086152467962186243905,
                           0,
                           0.26954315595234497233153198540086152467962186243905,
                           0.51909612920681181592572566945860955448022711511993,
                           0.73015200557404932409341625203115345804964306202613,
                           0.88706259976809529907515776930392726663167575122531,
                           0.97822865814605699280393800112285739077142240891978};
	static double weights[11] = {0.055668567116173666482753720442548578728515625696898,
                          0.1255803694649046246346942992239401001976157913954,
                          0.1862902109277342514260976414316558916912847480402,
                          0.23319376459199047991852370484317513943179817231696,
                          0.26280454451024666218068886989050919537276467760315,
                          0.27292508677790063071448352833634218915604196989478,
                          0.26280454451024666218068886989050919537276467760315,
                          0.23319376459199047991852370484317513943179817231696,
                          0.1862902109277342514260976414316558916912847480402,
                          0.1255803694649046246346942992239401001976157913954,
                          0.0556685671161736664827537204425485787285156256969};
	static double half_diff_kx = (bkx-akx)*0.5;
	static double half_diff_ky = (bky-aky)*0.5;
	static double half_diff_kz = (bkz-akz)*0.5;
	static double half_sum_kx = (bkx+akx)*0.5;
	static double half_sum_ky = (bky+aky)*0.5;
	static double half_sum_kz = (bkz+akz)*0.5;
	complex<double> integral[N][N] = {0};
	double end = 0;
	eigen TB = {};
	for(int l = 0; l < 11; l++){
		for(int m = 0; m < 11; m++){
			for(int n = 0; n < 11; n++){
				end = half_diff_kx*half_diff_ky*half_diff_kz*weights[l]*weights[m]*weights[n];
				TB = band_eigenvalues_and_eigenvectors(params, neighbors,  half_diff_kx*points[l]+half_sum_kx, half_diff_ky*points[m]+half_sum_ky, half_diff_kz*points[n]+half_sum_kz);		
				for(int s = 0; s < N; s++){
					integral[s][0] += TB.evectors[s][0]*sh[s][l+11*m+121*n]*end;
					integral[s][1] += TB.evectors[s][1]*sh[s][l+11*m+121*n]*end;
					integral[s][2] += TB.evectors[s][2]*sh[s][l+11*m+121*n]*end;
					integral[s][3] += TB.evectors[s][3]*sh[s][l+11*m+121*n]*end;
					integral[s][4] += TB.evectors[s][4]*sh[s][l+11*m+121*n]*end;
					integral[s][5] += TB.evectors[s][5]*sh[s][l+11*m+121*n]*end;
					integral[s][6] += TB.evectors[s][6]*sh[s][l+11*m+121*n]*end;
					integral[s][7] += TB.evectors[s][7]*sh[s][l+11*m+121*n]*end;
				}
			}
		}
	}

	for(int s = 0; s < N; s++){
		TB_Proj.projected_vec[s][0] = integral[s][0];
		TB_Proj.projected_vec[s][1] = integral[s][1];
		TB_Proj.projected_vec[s][2] = integral[s][2];
		TB_Proj.projected_vec[s][3] = integral[s][3];
		TB_Proj.projected_vec[s][4] = integral[s][4];
		TB_Proj.projected_vec[s][5] = integral[s][5];
		TB_Proj.projected_vec[s][6] = integral[s][6];
		TB_Proj.projected_vec[s][7] = integral[s][7];
		}
	return TB_Proj;
}

double residual(complex<double> sh[N][1331], double params[6], double neighbors[4][3], double fit_points[][3], int points[], double dft_bands[P][56], double dft_band_eigenvectors[P][56][N]){
	
	double diff_energies_squared = 0;
	double product = 0;
	double res = 0;
	complex<double> overlap = (0.0,0.0);
	eigen BAND_STRUCTURE = {};
	proj_bandvecs TB_Proj = projection_gauss_10_pt_quad(sh,	params, neighbors, -M_PI/a, M_PI/a, -M_PI/a, M_PI/a, -M_PI/a, M_PI/a);

	for(int k = 0; k < (my_sizeof(points)/my_sizeof(points[0])); k++){
		BAND_STRUCTURE = band_eigenvalues_and_eigenvectors(params, neighbors,  fit_points[k][0], fit_points[k][1], fit_points[k][2]);
		for(int i = 0; i < 56; i++){
			for(int j = 0; j < N; j++){
				diff_energies_squared = pow(dft_bands[points[k]][i]-BAND_STRUCTURE.evalues[j],2);
				overlap.real(0);
				overlap.imag(0);
				for(int element = 0; element < N; element++){
					overlap += TB_Proj.projected_vec[j][element]*dft_band_eigenvectors[points[k]][i][element];
				}
				product = diff_energies_squared*abs(norm(overlap));
				res += product;
			}
		}
	}

	return res;
}

void nelmin (complex<double> sh[N][1331], double fit_points[][3], int points[],
	int n, double start[10], double neighbors[4][3], double dft_bands[P][56], double dft_band_eigenvectors[P][56][N], double xmin[], 
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
    cout << string("REQMIN, N, or KONVGE has an illegal value.");
    return;
  }

  if ( n < 1 )
  {
    *ifault = 1;
    cout << string("REQMIN, N, or KONVGE has an illegal value.");
    return;
  }

  if ( konvge < 1 )
  {
    *ifault = 1;
    cout << string("REQMIN, N, or KONVGE has an illegal value.");

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

    y[n] = residual(sh, start, neighbors, fit_points, points, dft_bands, dft_band_eigenvectors);
    *icount = *icount + 1;

    for ( j = 0; j < n; j++ )
    {
      x = start[j];
      start[j] = start[j] + step[j] * del;
      for ( i = 0; i < n; i++ )
      {
        p[i+j*n] = start[i];
      }
      y[j] = residual(sh, start, neighbors, fit_points, points, dft_bands, dft_band_eigenvectors);
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
      ystar = residual(sh, pstar, neighbors, fit_points, points, dft_bands, dft_band_eigenvectors);
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
        y2star = residual(sh, p2star, neighbors, fit_points, points, dft_bands, dft_band_eigenvectors);
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
          y2star = residual(sh, p2star, neighbors, fit_points, points, dft_bands, dft_band_eigenvectors);
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
              y[j] = residual(sh, xmin, neighbors, fit_points, points, dft_bands, dft_band_eigenvectors);
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
          y2star = residual(sh, p2star, neighbors, fit_points, points, dft_bands, dft_band_eigenvectors);
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
          cout << string("iteration terminated because KCOUNT was exceeded without convergence.");
      break;
    }

    *ifault = 0;
    for ( i = 0; i < n; i++ )
    {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
      z = residual(sh, xmin, neighbors, fit_points, points, dft_bands, dft_band_eigenvectors);
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
            cout << string("iteration terminated because KCOUNT was exceeded without convergence.");
        break;
      }
      xmin[i] = xmin[i] - del - del;
      z = residual(sh, xmin, neighbors, fit_points, points, dft_bands, dft_band_eigenvectors);
      *icount = *icount + 1;
      if ( z < *ynewlo )
      {
        *ifault = 2;
            cout << string("iteration terminated because KCOUNT was exceeded without convergence.");
        break;
      }
      xmin[i] = xmin[i] + del;
    }

    if ( *ifault == 0 )
    {
        cout << string("no errors detected.");
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

double fRand(double fMin, double fMax)
{
    double f = (double)rand() / RAND_MAX;
    return fMin + f * (fMax - fMin);
}

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
	double dft_energy = 0;
	double tmp_dft_bands[900][56];
	double kx, ky, kz;
	double tmp_dft_path[900][3];
	double tmp_dft_eigenvectors[900][56][N];
	double s, py, pz, px, dxy, dyz, dz2, dxz, dx2;
	string throwaway;
	PROCAR >> throwaway >> throwaway >> throwaway;
	PROCAR >> throwaway >> throwaway >> throwaway;
	PROCAR >> num_kpoints >> throwaway >> throwaway >> throwaway >> num_bands >> throwaway >> throwaway >> throwaway >> num_ions;
	PROCAR >> throwaway >> throwaway >> throwaway >> kx >> ky >> kz >> throwaway >> throwaway >> throwaway;
	int counter = 4;
	int element = 0;
	int e_element = 0;
	int line = (num_ions + 4)*num_bands + 1;
	while(PROCAR){	
		if((counter-4)%line==0){
			//k-points
			if(counter-4 != 0){
				tmp_dft_path[element][0] = kx;
				tmp_dft_path[element][1] = ky;
				tmp_dft_path[element][2] = kz;
				element++;
				e_element = 0;
				PROCAR >> throwaway >> throwaway >> throwaway >> kx >> ky >> kz >> throwaway >> throwaway >> throwaway;
			}
		}
		if(((counter-4)%line)%(num_ions+3)==1){
			PROCAR >> throwaway >> throwaway >> throwaway >> throwaway >> dft_energy >> throwaway >> throwaway >> throwaway;
			tmp_dft_bands[element][e_element] = dft_energy;

		}
		if((((counter-4)%line)%(num_ions+3)==0) && (counter-4)%line!=0){
			PROCAR >> throwaway >> s >> py >> pz >> px >> dxy >> dyz >> dz2 >> dxz >> dx2 >> throwaway;
			tmp_dft_eigenvectors[element][e_element][0] = s;
			tmp_dft_eigenvectors[element][e_element][1] = py;
			tmp_dft_eigenvectors[element][e_element][2] = pz;
			tmp_dft_eigenvectors[element][e_element][3] = px;
			tmp_dft_eigenvectors[element][e_element][4] = dxy;
			tmp_dft_eigenvectors[element][e_element][5] = dyz;
			tmp_dft_eigenvectors[element][e_element][6] = dz2;
			tmp_dft_eigenvectors[element][e_element][7] = dxz;
			e_element++;
		}
		if(!((counter-4)%line==0) && !(((counter-4)%line)%(num_ions+3)==1) && !((((counter-4)%line)%(num_ions+3)==0) && (counter-4)%line!=0)){
			PROCAR >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway >> throwaway;
		}
		counter+=1;
	}

	double dft_path[P][3] = {0};
	double dft_bands[P][56] = {0};
	double dft_band_eigenvectors[P][56][N] = {0};

	counter = 0;
	for(int k = 0; k < P; k++){
		if(k == 100 || k == 199 || k == 299 || k == 399 || k == 499 || k == 599 || k == 699 || k == 799){
			counter++;
		}
		dft_path[k][0] = tmp_dft_path[counter][0];
		dft_path[k][1] = tmp_dft_path[counter][1];
		dft_path[k][2] = tmp_dft_path[counter][2];
		for(int i = 0; i < 56; i++){
			dft_bands[k][i] = tmp_dft_bands[counter][i];
			for(int j = 0; j < 8; j++){
				dft_band_eigenvectors[k][i][j] = tmp_dft_eigenvectors[counter][i][j];
			}
		}
		counter++;
	}
	
	double ticks[P];
	ticks[0] = 0;
	double tick = 0;

	for(int i = 1; i < P; i++){
		if(pow(pow(dft_path[i][0]-dft_path[i-1][0],2)+pow(dft_path[i][1]-dft_path[i-1][1],2)+pow(dft_path[i][2]-dft_path[i-1][2],2),0.5) < 1.0){
			tick += pow(pow(dft_path[i][0]-dft_path[i-1][0],2)+pow(dft_path[i][1]-dft_path[i-1][1],2)+pow(dft_path[i][2]-dft_path[i-1][2],2),0.5);
		}
		else{
			tick += 0.0000000000001;
		}

		ticks[i] = tick;
	}

	// for(int k = 0; k < 892; k++){
	// 	cout <<dft_path[k][0] << " " <<dft_path[k][1] << " " <<dft_path[k][2] << " " <<endl; 
	// 	for(int i = 0; i < 56; i++){
	// 		cout << dft_bands[k][i] << endl;
	// 		for(int j = 0; j < 8; j++){
	// 			cout << dft_band_eigenvectors[k][i][j] << " ";
	// 		}
	// 		cout << endl;
	// 	}
	// }

	///////////////////////////////////////////////////////////////////////////////////////
	//// Read Fourier transformed spherical harmonics at gauss quad points
	///////////////////////////////////////////////////////////////////////////////////////

	complex<double> sh[N][1331] = {0};
	double num;

	ifstream SH;
	SH.open("s_real.dat");
	counter = 0;
	SH >> num;
	sh[0][0].real(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[0][counter].real(num);
	}
	SH.close();

	SH.open("py_real.dat");
	counter = 0;
	SH >> num;
	sh[1][0].real(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[1][counter].real(num);
	}
	SH.close();

	SH.open("pz_real.dat");
	counter = 0;
	SH >> num;
	sh[2][0].real(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[2][counter].real(num);
	}
	SH.close();

	SH.open("px_real.dat");
	counter = 0;
	SH >> num;
	sh[3][0].real(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[3][counter].real(num);
	}
	SH.close();

	SH.open("dxy_real.dat");
	counter = 0;
	SH >> num;
	sh[4][0].real(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[4][counter].real(num);
	}
	SH.close();

	SH.open("dyz_real.dat");
	counter = 0;
	SH >> num;
	sh[5][0].real(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[5][counter].real(num);
	}
	SH.close();

	SH.open("dz2_real.dat");
	counter = 0;
	SH >> num;
	sh[6][0].real(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[6][counter].real(num);
	}
	SH.close();

	SH.open("dxz_real.dat");
	counter = 0;
	SH >> num;
	sh[7][0].real(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[7][counter].real(num);
	}
	SH.close();

	SH.open("s_imag.dat");
	counter = 0;
	SH >> num;
	sh[0][0].imag(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[0][counter].imag(num);
	}
	SH.close();

	SH.open("py_imag.dat");
	counter = 0;
	SH >> num;
	sh[1][0].imag(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[1][counter].imag(num);
	}
	SH.close();

	SH.open("pz_imag.dat");
	counter = 0;
	SH >> num;
	sh[2][0].imag(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[2][counter].imag(num);
	}
	SH.close();

	SH.open("px_imag.dat");
	counter = 0;
	SH >> num;
	sh[3][0].imag(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[3][counter].imag(num);
	}
	SH.close();

	SH.open("dxy_imag.dat");
	counter = 0;
	SH >> num;
	sh[4][0].imag(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[4][counter].imag(num);
	}
	SH.close();

	SH.open("dyz_imag.dat");
	counter = 0;
	SH >> num;
	sh[5][0].imag(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[5][counter].imag(num);
	}
	SH.close();

	SH.open("dz2_imag.dat");
	counter = 0;
	SH >> num;
	sh[6][0].imag(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[6][counter].imag(num);
	}
	SH.close();

	SH.open("dxz_imag.dat");
	counter = 0;
	SH >> num;
	sh[7][0].imag(num);
	while(SH){
		counter++;
		SH >> num; 
		sh[7][counter].imag(num);
	}
	SH.close();

	////////////////////////////////////////////////////////////////////////////////////////
	//**************************************************************************************
	////////////////////////////////////////////////////////////////////////////////////////

	///////////////////////////////////////////////////////////////////////////////////////
	//// Plot DFT Bandstructure
	///////////////////////////////////////////////////////////////////////////////////////
	// ofstream bandvecs_and_values_data_dft;
	// bandvecs_and_values_data_dft.open("BAND_STRUCTURE_DFT.dat");
	// counter = 0;

	// for(int j = 0; j < 56; j++){
	// 	for(int k = 0; k < P; k++){
	// 		bandvecs_and_values_data_dft << ticks[k] << "\t" << dft_bands[k][j] << endl;
	// 	}
	// 	bandvecs_and_values_data_dft << endl << endl;
	// }
	// bandvecs_and_values_data_dft.close();
	// string command= string("gnuplot -e \"") +
	// // string("set terminal postscript enhanced;") +
	// // string("set output \"./Band-Structure.png\";") +
	// string("unset key;") +
	// string("set xtics (\'L\' 0,\'G\' 0.866025,\'X\' 1.86603,\'W\' 3.36603,\'U,K\' 4.6279,\'L\' 5.24028,\'W\' 5.94738,\'X\' 7.44738,\'U,K\' 7.79737,\'G\' 8.85803);") +
	// string("plot for [i=0:55]\\\"BAND_STRUCTURE_DFT.dat\\\" index i u 1:2 w l;") + 
	// string("pause mouse\" 2>/dev/null;");
	// // string("eog ./Band-Structure.png");
	// const char* com = command.c_str();
	// system(com);

	////////////////////////////////////////////////////////////////////////////////////////
	//**************************************************************************************
	////////////////////////////////////////////////////////////////////////////////////////

	double neighbors[4][3] = {
	{0.25*a,0.25*a,0.25*a},
	{0.25*a,-0.25*a,-0.25*a},
	{-0.25*a,0.25*a,-0.25*a},
	{-0.25*a,-0.25*a,0.25*a}};

	// params es, ep, vss, vsp, vppS, vppP, vss2, vsp2, vppS2, vppP2
	// vxx = (4/3)(vppS+2vppP), vxy = (4/3)(vppS-vppP), vppS = (1/4)(vxx+2vxy), vppP = (1/4)(vxx-vxy)
	double es = 5.695;
	double ep = 7.2 + es;
	double params[6] = {es, ep,-8.13, 5.88, 1.71 , 7.51};
	//  5.69502933
	// 14.89480613
	// -8.12985530
	//  5.87977066
	//  1.71005241
	//  7.51014784

	// double params[6] = {es, ep,-8.23, 5.9134, 2.0471 , 4.2553};
	// double params[6] = {5.695,	13.0977,	-8.22962,	6.04601,	1.62654,	7.56112};
	// double params[6] = {1,-3.40277,-1.77416,-0.84730,1.53114,-0.716431,-0.0891011,-0.105639,0.388664,-0.126518};
	// double params[6] ={1.39458,	2.59735,	-4.0652,	-3.65098,	0.202101,	-4.21768,	-4.30094,	-2.95345,	-0.385795,	3.19677};
	// double params[6] = {0.567518,	3.28803,	12.9791,	-0.821661,	2.32315,	-11.4802,	-7.48998,	0.730091,	-0.500956,	3.66033	};
	// double params[6] = {1.21926,	2.64689,	10.9146,	7.39687,	5.84823,	-13.0186};
	
	// double fit_points[5][3] = {{0,0,0},
	// 							{0.5,0.5,0.5},
	// 							{0.0,0.0,1.0},
	// 							{0.257576,0.252525,0.989899},
	// 							{1.0,0.5,0.0}};
	// int points[5] = {99,0,198,395,297};

	int points[P];
	for(int p = 0; p < P;p++){	
		points[p] = p;
	}
	
	// double res =  residual(sh, params, neighbors, dft_path, points, dft_bands, dft_band_eigenvectors);
	// double min_res = res;
	// double new_params[6];	
	// int index = 0;
	// for(int n =0; n <1000; n++){
	// 	index = (rand()%5)+1;
	// 	for(int j = 0; j<6;j++){
	// 		new_params[j] = params[j];
	// 	}
	// 	new_params[index] += fRand(-0.1,0.1);
		
	// 	res = residual(sh, new_params, neighbors, dft_path, points, dft_bands, dft_band_eigenvectors);
	// 	if(res < min_res){
	// 		cout << n << endl;
	// 		min_res = res;
	// 		cout << res << endl;
	// 		cout << "{";
	// 		for(int l = 0; l<6;l++){
	// 			cout << params[l] << ",\t";
	// 			params[l] = new_params[l];
	// 		}
	// 		cout << "}"<<endl;
	// 	}
	// }

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
	// n = 6;
	// start = new double[n];	
	// step = new double[n];
	// xmin = new double[n];
 // 	reqmin = 1.0E-08;
	// konvge = 10;
	// kcount = 10000;
	// cout << "\n";
	// cout << "  Starting point X:\n";
	// cout << "\n";
	// for (int  i = 0; i < n; i++ )
	// {	
	// 	step[i] = 2;
	// 	start[i] = params[i];
	// 	cout << "  " << setw(14) << start[i] << "\n";
	// }
	// ynewlo = residual(sh, params, neighbors,  dft_path, points, dft_bands, dft_band_eigenvectors);
	// cout << "\n";
	// cout << "  F(X) = " << fixed << setprecision(8) << ynewlo << "\n";



	// nelmin (sh, dft_path, points, n, start, neighbors, dft_bands, dft_band_eigenvectors, xmin, &ynewlo, reqmin, step, konvge, kcount, &icount, &numres, &ifault );
	// cout << "\n";
	// cout << "  Return code IFAULT = " << ifault << "\n";
	// cout << "\n";
	// cout << "  Estimate of minimizing value X*:\n";
	// cout << "\n";
	// for (int i = 0; i < n; i++ )
	// {
	// cout << "  " << setw(14) << xmin[i] << "\n";
	// params[i] = xmin[i];
	// }

	// cout << "\n";
	// cout << "  F(X*) = " << fixed << setprecision(8) <<ynewlo << "\n";

	// cout << "\n";
	// cout << "  Number of iterations = " << icount << "\n";
	// cout << "  Number of restarts =   " << numres << "\n";

	// delete [] start;
	// delete [] step;
	// delete [] xmin;

	// ///////////////////////////////////////////////////////////////////////////////////////
	// //// Plot TB Bandstructure
	// ///////////////////////////////////////////////////////////////////////////////////////

	eigen BANDS; 

	double bands[P][N] = {0};

	for(int i = 0; i < P; i++){
		BANDS = band_eigenvalues_and_eigenvectors(params, neighbors,  dft_path[i][0], dft_path[i][1], dft_path[i][2]);
		for(int j = 0; j < N; j++){
			bands[i][j] = BANDS.evalues[j];
			cout << bands[i][j] << "\t";
		}
		cout << endl;
	}

	ofstream bandvecs_and_values_data_tb;
	bandvecs_and_values_data_tb.open("BAND_STRUCTURE_TB.dat");

	for(int j = 0; j < N; j++){
		for(int k = 0; k < P; k++){
			bandvecs_and_values_data_tb << ticks[k] << "\t" << bands[k][j] << endl;
		}
		bandvecs_and_values_data_tb << endl << endl;
	}

	// for(int k = 0; k < P; k++){
	// 	cout << k << " " << ticks[k] << " "  << dft_path[k][0] <<" "  << dft_path[k][1] <<" "  << dft_path[k][2] << endl; 
	// }
	bandvecs_and_values_data_tb.close();
	string command2= string("gnuplot -e \"") +
	// string("set terminal postscript enhanced;") +
	// string("set output \"./Band-Structure.png\";") +
	string("unset key;") +
	string("set yrange [-10:27];") +
	string("set xtics (\'L\' 0,\'G\' 0.866025,\'X\' 1.86603,\'W\' 3.36603,\'U,K\' 4.6279,\'L\' 5.24028,\'W\' 5.94738,\'X\' 7.44738,\'U,K\' 7.79737,\'G\' 8.85803);") +
	// string("plot \\\"BAND_STRUCTURE_DFT.dat\\\" u 1:2 w l lc 7 lw 1,") +
	string("plot \\\"BAND_STRUCTURE_TB.dat\\\" u 1:2 w l lc 8 lw 2;") +  
	string("pause mouse\"");
	// string("eog ./Band-Structure.png");
	const char* com2 = command2.c_str();
	system(com2);

	////////////////////////////////////////////////////////////////////////////////////////
	//**************************************************************************************
	////////////////////////////////////////////////////////////////////////////////////////

	return 0;
}
