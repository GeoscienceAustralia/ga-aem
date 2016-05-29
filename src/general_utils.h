/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _general_utils_H
#define _general_utils_H

#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <cstring>
#include <ctime>
#include <climits>
#include <vector>
#include <algorithm>
#include "general_constants.h"
#include "general_types.h"

std::string commandlinestring(int argc, char** argv);
std::string versionstring(const std::string& version, const std::string& compiletime, const std::string& compiledate);

bool mpi_initialised();
std::string mpi_processername();
int my_size();
int my_rank();
void my_barrier();

std::string strprint_va(const char* fmt, va_list argptr);
std::string strprint(const char* fmt, ...);

void message(const char* fmt, ...);
void message(FILE* fp,const char* fmt, ... );
void rootmessage(const char* fmt, ...);
void rootmessage(FILE* fp, const char* fmt, ...);
void warningmessage(const char* fmt, ...);
void errormessage(const char* fmt, ... );
void warningmessage(FILE* fp, const char* fmt, ...);
void errormessage(FILE* fp, const char* fmt, ...);

void prompttocontinue();
void prompttoexit();

void rb_sleep(double secs);
void debug(const char* msg);

template <typename T>
int roundnearest(const T& x, int nearest)
{
	return nearest * (int) std::round(x / nearest);
}

template <typename T>
T roundnearest(const T& x, const double& nearest)
{
	return (T)(nearest * std::round(x / nearest) );
}

template <typename T>
T roundupnearest(const T& x, const double& nearest)
{
	return (T)(nearest * std::ceil(x / nearest) );
}

template <typename T>
T rounddownnearest(const T& x, const double& nearest)
{
	T f = std::floor(x / nearest);
	T v = nearest*f;
	return v;
}

template <typename T> 
T pow10(const T& x)
{
	return std::pow((T)10.0, x);
};

template <typename T>
T distance(const T& x1, const T& y1, const T& x2, const T& y2)
{
	return std::sqrt(std::pow(x2 - x1, 2.0) + std::pow(y2 - y1, 2.0));
}

template <typename T>
T distance(const T& x, const T& y)
{
	return std::sqrt(x*x + y*y);
}

//==========
bool wildcmp(const char* wildpattern, const char* stringpattern);
double correlation_coefficient(std::vector<double>x, std::vector<double>y);

template <typename T>
T swap_endian(const T u)
{
	//from stackoverflow.com/questions/105252/how-do-i-convert-between-big-endian-and-little-endian-values-in-c	

	static_assert (CHAR_BIT == 8, "CHAR_BIT != 8");
	union
	{
		T u;
		unsigned char u8[sizeof(T)];
	} source, dest;
	source.u = u;
	for (size_t k = 0; k < sizeof(T); k++){
		dest.u8[k] = source.u8[sizeof(T) - k - 1];
	}
	return dest.u;
}

template <typename T>
void swap_endian(T* array, size_t num)
{
	for (size_t i = 0; i < num; i++){
		array[i] = swap_endian(array[i]);
	}
}

template <typename T>
void swap_endian(std::vector<T>& array)
{
	for (size_t i = 0; i < array.size(); i++){
		array[i] = swap_endian(array[i]);
	}
}

int log10stretch(double val, double lowclip, double highclip);
double inverselog10stretch(int bin, double lowclip, double highclip);
int linearstretch(double val, double lowclip, double highclip);
double inverselinearstretch(int bin, double lowclip, double highclip);
bool regression(double* x, double*y, size_t n, double* gradient, double* intercept);
bool regression(const std::vector<double>& x, const std::vector<double>& y, double& gradient, double& intercept);
bool bestfitlineendpoints(const std::vector<double>& x, const std::vector<double>& y, double& x1, double& y1, double& x2, double& y2);

std::string stringvalue(const double value, const char* fmt=(char*)NULL);
std::string stringvalue(const size_t value, const char* fmt = (char*)NULL);
std::string stringvalue(const int value,    const char* fmt=(char*)NULL);
std::string stringvalue(const bool value);

int strcasecmp(const std::string& A, const std::string& B);
int strncasecmp(const std::string& A, const std::string& B, size_t n);

std::vector<std::string> fieldparsestring_old(const char* str, const char delim);
std::vector<std::string> fieldparsestring(const char* str, const char* delims);
std::vector<double> getdoublevector(const char* str, const char* delims);

std::string timestamp();
std::string timestring(const std::string format, std::time_t t = 0);

int isinsidepolygon(int npol, double *xp, double *yp, double x, double y);

bool eq(double& a, double& b);
bool lt(double& a, double& b);
bool gt(double& a, double& b);
bool le(double& a, double& b);
bool ge(double& a, double& b);
void planeequation(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2, const double& x3, const double& y3, const double& z3, double& A, double& B, double& C, double& D);

int findindex(const size_t n, const double* x, const double& xtarget);
int findindex(const std::vector<double>& x, const double& xtarget);
double linearinterp(const double& x1, const double& y1, const double& x2, const double& y2, const double& x);
double linearinterp(const std::vector<double>& x, const std::vector<double>& y, const double& xtarget);
double linearinterp(const size_t n, const double* x, const double* y, const double& xtarget);
void   linearinterp(const size_t n, const double* x, const double* y, size_t ni, const double* xi, double* yi);
std::vector<double> linearinterp(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xi);

std::vector<sRange> parserangelist(std::string& str);
std::vector<std::string> parsestrings(const std::string& str, const std::string& delims);

bool isreportable(int rec);
bool isinrange(const sRange& r, const int& i);

double overlap(const double al, const double ah, const double bl, const double bh);
double fractionaloverlap(const double al, const double ah, const double bl, const double bh);

std::vector<double> overlaps(const double& a1, const double& a2, const std::vector<double>& b);
std::vector<double> fractionaloverlaps(const double& a1, const double& a2, const std::vector<double>& b);

std::vector<std::vector<double>> overlaps(const std::vector<double>& a, const std::vector<double>& b);
std::vector<std::vector<double>> fractionaloverlaps(const std::vector<double>& a, const std::vector<double>& b);

int findandreplace(std::string& source, const std::string& find, const std::string& replace);
double gettime();

char* temppath(const char* s, int set);
void settemppath(const char* s);
std::string gettemppath();

template<typename T>
inline const T sign(const T &a, const T &b) {return b >= 0 ? (a >= 0 ? a : -a) : (a >= 0 ? -a : a);}

void sort(char** strings, const size_t n);
void sort(int* x, const size_t n);
void sort(float* x, const size_t n);
void sort(double* x, const size_t n);
void indexsort(double* x, const size_t n);

template<typename T>
double covariance(const std::vector<T>& x, const std::vector<T>& y)
{
	size_t n = x.size();
	cStats<T> sx(x);
	cStats<T> sy(y);
	std::vector<T> v(n);
	for (size_t i = 0; i < x.size(); i++){
		v[i] = (x[i] - sx.mean) * (y[i] - sy.mean);
	}
	cStats<T> sv(v);	
	return sv.mean;
}

template<typename T>
double correlation(const std::vector<T>& x, const std::vector<T>& y)
{
	size_t n = x.size();
	cStats<T> sx(x);
	cStats<T> sy(y);
	std::vector<double> v(n);
	for (size_t i = 0; i < x.size(); i++){
		v[i] = (x[i] - sx.mean) / sx.std * (y[i] - sy.mean) / sy.std;
	}	
	cStats<T> sv(y);
	return sv.mean;
}

double reportusage();
double percentmemoryused();
void guage(int ntot, int n, int pdiv1, int pdiv2);
unsigned int factorial(unsigned int n);

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems);
std::vector<std::string> split(const std::string &s, char delim);


void allocate1darray(int*& a,const size_t n);
void allocate1darray(double*& a, const size_t n);
void deallocate1darray(double*& a);
void deallocate1darray(int*& a);
void allocate2darray(double**& a, const size_t nrows, const size_t ncols);
void deallocate2darray(double**& a, const size_t nrows);

double median(double* v,size_t n);
std::string trim(const std::string& s);
std::string stripquotes(const std::string& s);
std::vector<std::string> tokenize(const std::string& str);

#define SORT_UP 0
#define SORT_DOWN 1
template<typename T> void quicksortindex(T* a, int* index, const int& leftarg, const int& rightarg, int sortupordown)
{	
	if (leftarg < rightarg) {
		T pivotvalue = a[leftarg];
		int left  = leftarg - 1;
		int right = rightarg + 1;
		for(;;) {
			if(sortupordown==SORT_UP){
				while (a[--right] > pivotvalue);
				while (a[++left] < pivotvalue);
			}
			else{
				while (a[--right] < pivotvalue);
				while (a[++left] > pivotvalue);
			}
			if (left >= right) break;			

			T temp   = a[right];
			a[right] = a[left];
			a[left]  = temp;

			if(index){
				int tempind  = index[right];
				index[right] = index[left];
				index[left]  = tempind;
			}
		}
		int pivot = right;
		quicksortindex(a, index, leftarg, pivot, sortupordown);
		quicksortindex(a, index, pivot + 1, rightarg, sortupordown);
	}
}

template<typename T>
bool bwrite(FILE* fp, const T& v){
	size_t status = fwrite(&v, sizeof(T), 1, fp);
	if (status != 1){
		printf("Error in bwrite() writing to binary file\n");
		throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
	}
	return true;
}

template<typename T>
bool bwrite(FILE* fp, const std::vector<T>& v){
	size_t status = fwrite(&(v[0]), sizeof(T), v.size(), fp);
	if (status != v.size()){
		printf("Error in bwrite() writing to binary file\n");
		throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
	}
	return true;
}

std::vector<float> dvec2fvec(std::vector<double>& vd);

std::vector<double> linspace(const double x1, const double x2, const size_t n);
std::vector<double> log10space(const double x1, const double x2, const size_t n);

bool filegetline(FILE* fp, std::string& str);

int LevenshteinDistance(char* s, int len_s, char* t, int len_t);

#endif

