/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _RANDOM_H_
#define _RANDOM_H_
#include <vector>
#include "matrix_ops.h"
using namespace std;

void seedrand();
void seedrand(unsigned int seed);
double urand();
double urand(double rmin, double rmax);
size_t irand(size_t imin, size_t imax);
int    irand(int imin, int imax);
double nrand();
std::vector<double> nrand(size_t n);
void nrand(size_t n, double* x);
void nrand(size_t n, double* x, double mean, double std);

std::vector<double> mvnrand_lowercholesky(const MatrixDouble& L);
std::vector<double> mvnrand_covariance(const MatrixDouble& C);

double gaussian_pdf(double mean, double std, double x);
double mvgaussian_pdf(const VectorDouble& m0, const MatrixDouble& C, const VectorDouble& m);

#endif
