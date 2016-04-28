/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstdio>
#include <cstdlib>
#include <cmath>
#include <ctime>
#include <cfloat>

#include "general_constants.h"
#include "matrix_ops.h"
#include "random.h"

void seedrand()
{			
	srand((unsigned int)time(NULL));	
}
void seedrand(unsigned int seed)
{		
	srand(seed);	
}
double urand()
{		
	double u = (double)rand() / ((double)RAND_MAX + 1.0);
	if (u == 0.0){
		u = DBL_EPSILON;
	}
	return u;
}
double urand(double rmin, double rmax)
{	
	return rmin + urand()*(rmax-rmin);
}
size_t irand(size_t imin, size_t imax)
{
	return imin + rand() % (imax - imin + 1);
}
int irand(int imin, int imax)
{	
	return imin + rand() % (imax-imin+1);	
}
double nrand()
{	
	//Box-Muller transform
	//http://en.wikipedia.org/wiki/Box%E2%80%93Muller_transform
	double u,v,s;
	
	do{
		u=urand(-1.0,1.0);
		v=urand(-1.0,1.0);
		s=u*u+v*v;
	}while(s==0 || s>=1.0);

	return u*sqrt(-2.0*log(s)/s);

}
std::vector<double> nrand(size_t n)
{
	std::vector<double> x(n);
	for (size_t i=0; i<n; i++){
		x[i] = nrand();
	}
	return x;
}
void nrand(size_t n, double* x)
{	
	for(size_t i=0;i<n;i++){
		x[i]=nrand();
	}
}
void nrand(size_t n, double* x, double mean, double std)
{	
	for(size_t i=0;i<n;i++){
		x[i]=nrand()*std + mean;
	}
}
std::vector<double> mvnrand_lowercholesky(const MatrixDouble& L)
{
	//L - lower diagonal of the Cholesky decomposition
	size_t n = L.num_rows();
	VectorDouble u(n);
	nrand(n,u,0.0,1.0);
	VectorDouble v = L*u;
	std::vector<double> x(n);
	for(size_t i=0;i<n;i++){
		x[i]=v[i];
	}
	return x;
}
std::vector<double> mvnrand_covariance(const MatrixDouble& C)
{
	MatrixDouble L = lower_cholesky(C);
	return mvnrand_lowercholesky(L);
}
double gaussian_pdf(double mean, double std, double x)
{			
	double p = exp(-0.5 * pow((x-mean)/std,2.0) ) / (sqrt(TWOPI)*std);	
	return p;
}
double mvgaussian_pdf(const VectorDouble& m0, const MatrixDouble& C, const VectorDouble& m)
{				
	size_t k = m0.dim();		
	MatrixDouble I = identitymatrix(k);

	TNT::Linear_Algebra::LU<double> lu(C);	
	MatrixDouble invC = lu.solve(I);	
	double detC = lu.det();

	VectorDouble dm = m-m0;
	double a = -0.5*mtAm(dm,invC);
	double pdf = exp(a) / sqrt(pow(TWOPI,k)*detC);
	return pdf;
}






