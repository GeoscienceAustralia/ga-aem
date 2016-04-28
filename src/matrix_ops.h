/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _matrix_ops_H
#define _matrix_ops_H
#include <vector>

#define TNT_SUBSCRIPT_TYPE size_t
#include "tnt_subscript.h"

#include "tnt_vector.h"
#include "tnt_matrix.h"

//Stop Visual Studio and GNU compilers from throwwing warnings in tnt_linalg.h (tnt_linalg.h should use TNT::Subscript instead of int as an index)
#ifdef _WIN32	
	#pragma warning( push, 3 )
	#pragma warning( disable : 4267 )
		#include "tnt_linalg.h"
	#pragma warning( default : 4267 )	
	#pragma warning( pop )
#else	
	#pragma GCC diagnostic push
	#pragma GCC diagnostic ignored "-Wsign-compare"
		#include "tnt_linalg.h"
	#pragma GCC diagnostic pop
#endif

typedef TNT::Vector<double> VectorDouble;
typedef TNT::Matrix<double> MatrixDouble;
///////////////////////////////////////////////////////////////////////////

void writetofile(const MatrixDouble& A, std::string path);
void writetofile(const VectorDouble& x, std::string path);
void writetofile(const std::vector<double>& x, std::string path);

std::vector<double> operator*(const MatrixDouble& A, const std::vector<double>& x);
double dot(const VectorDouble& a, const VectorDouble& b);
double mtDm(const VectorDouble& m, const MatrixDouble& D);
double mtDm(const std::vector<double>& m, const MatrixDouble& D);
double mtAm(const VectorDouble& m, const MatrixDouble& A);
double mtAm(const std::vector<double>& m, const MatrixDouble& A);

MatrixDouble submatrix(const MatrixDouble& X, const TNT::Subscript i0, const TNT::Subscript i1, const TNT::Subscript j0, const TNT::Subscript j1);
MatrixDouble pseudoinverse(const MatrixDouble& X);
MatrixDouble pseudoinverse_od(const MatrixDouble& X);
MatrixDouble lower_cholesky(const MatrixDouble& A);
MatrixDouble identitymatrix(TNT::Subscript n);
MatrixDouble inverse(const MatrixDouble& A);
void print(const MatrixDouble& A, const std::string& name="");

#endif


