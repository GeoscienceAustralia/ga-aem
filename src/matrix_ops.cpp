/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cmath>
#include "general_utils.h"
#include "file_utils.h"

#include "matrix_ops.h"
using namespace TNT;

void writetofile(const MatrixDouble& A, std::string path)
{	
	FILE* fp = fileopen(path.c_str(),"w");
	for(Subscript i=0; i<A.num_rows(); i++){
		for(Subscript j=0; j<A.num_cols(); j++){
			fprintf(fp,"%lu\t%lu\t%lg\n",i+1,j+1,A[i][j]);
		}
	}
	fclose(fp);
	
}
void writetofile(const VectorDouble& x, std::string path) 
{	
	FILE* fp = fileopen(path.c_str(),"w");
	for(Subscript i=0; i<x.dim(); i++){		
		fprintf(fp,"%lu\t%lg\n",i+1,x[i]);	
	}
	fclose(fp);
}
void writetofile(const std::vector<double>& x, std::string path) 
{	
	FILE* fp = fileopen(path.c_str(),"w");
	for(Subscript i=0; i<x.size(); i++){		
		fprintf(fp,"%lu\t%lg\n",i+1,x[i]);	
	}
	fclose(fp);
}

std::vector<double> operator*(const MatrixDouble& A, const std::vector<double>& x)
{	
	Subscript M=A.num_rows();
	Subscript N=A.num_cols();
	std::vector<double> b(M);	
	for(Subscript i=0;i<M;i++){
		b[i]=0.0;
		for(Subscript j=0;j<N;j++){
			b[i] += A[i][j]*x[j];
		}
	}
	return b;
}
double dot(const VectorDouble& a, const VectorDouble& b)
{
	double sum=0.0;
	for(Subscript k=0;k<a.dim();k++){
		sum += a[k]*b[k];
	}
	return sum;
}
double mtDm(const VectorDouble& m, const MatrixDouble& D)
{
	//D must be diagonal
	Subscript n=m.dim();
	double sum=0.0;
	for(Subscript i=0;i<n;i++)sum+=(m[i]*m[i]*D[i][i]);
	return sum;
	
}
double mtDm(const std::vector<double>& m, const MatrixDouble& D)
{
	//D must be diagonal
	Subscript n = (Subscript)m.size();
	double sum=0.0;
	for(Subscript i=0;i<n;i++){
		sum+=(m[i]*m[i]*D[i][i]);
	}
	return sum;	
}
double mtAm(const VectorDouble& m, const MatrixDouble& A)
{					
	return dot(m,A*m);	
}
double mtAm(const std::vector<double>& m, const MatrixDouble& A)
{	
	std::vector<double> a=A*m;
	double sum=0.0;
	Subscript n = (Subscript)m.size();
	for(Subscript i=0; i<n; i++){
		sum += m[i]*a[i];
	}
	return sum;	
}

MatrixDouble submatrix(const MatrixDouble& X, const Subscript i0, const Subscript i1, const Subscript j0, const Subscript j1)
{
	Subscript m = i1-i0+1;
	Subscript n = j1-j0+1;
	MatrixDouble M(m,n);
	for(Subscript i=0; i<m; i++){
		for(Subscript j=0; j<n; j++){
			M[i][j] = X[i0+i][j0+j];
		}
	}
	return M;
}

MatrixDouble pseudoinverse(const MatrixDouble& X)
{
	if(X.num_rows() >= X.num_cols()){
		return pseudoinverse_od(X);
	}
	else{
		return transpose(pseudoinverse_od(transpose(X)));
	}
}
MatrixDouble pseudoinverse_od(const MatrixDouble& X)
{
	//Only for nrows >= ncols
	assert(X.num_rows() >= X.num_cols());
	
	MatrixDouble U;//m x n	
	MatrixDouble V;//n x n 		
	VectorDouble s;

	TNT::Linear_Algebra::SVD<double> svd(X);
	svd.getU(U);	
	svd.getV(V);	
	svd.getSingularValues(s);
	
	//pinvA = V     (1/S)    U'
	//m     = (nxn) (nxn)  (nxm)
	Subscript rank = (Subscript)svd.rank();
	Subscript nrows = V.num_rows();
	Subscript ncols = U.num_rows();
	MatrixDouble P(nrows,ncols);
	for(Subscript i=0;i<nrows;i++){
		for(Subscript j=0;j<ncols;j++){
			double sum=0.0;
			for(Subscript k=0;k<rank;k++){				
				sum += V[i][k] * U[j][k] / s[k];
			}
			P[i][j] = sum;
		}
	}			
	return P;
}

MatrixDouble lower_cholesky(const MatrixDouble& A)
{	
	TNT::Linear_Algebra::Cholesky<double> chol(A);	
	MatrixDouble L = chol.getL();
	return L;
}
MatrixDouble identitymatrix(Subscript n)
{	
	MatrixDouble I(n,n,0.0);
	for(Subscript k=0;k<n;k++)I[k][k]=1.0;
	return I;	
}
MatrixDouble inverse(const MatrixDouble& A)
{	
	Subscript m=A.num_rows(); 		
	MatrixDouble I = identitymatrix(m);
	TNT::Linear_Algebra::LU<double> lu(A);
	MatrixDouble B = lu.solve(I);
	return B;
}

void print(const MatrixDouble& A, const std::string& name)
{
	printf(name.c_str());
	printf("\n");
	for(Subscript i=0;i<A.num_rows();i++){
		for(Subscript j=0;j<A.num_cols();j++){
			printf("%10lf  ",A[i][j]);
		}
		printf("\n");
	}
}



