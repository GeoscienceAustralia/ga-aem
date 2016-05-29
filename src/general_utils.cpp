/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstdio>
#include <cfloat>
#include <cstring>
#include <ctime>
#include <sys/timeb.h>
#include <iterator>
#include <sstream>

#if defined _WIN32
	#define NOMINMAX 
	#include <windows.h>
	#include <conio.h>
#else
	#include <unistd.h>	
	#include <sys/types.h>	
	#include <sys/resource.h>	
	#include <errno.h>
#endif


#if defined _MPI_ENABLED
#include <mpi.h>
#endif

#if defined _OPENMP
#include <omp.h>
#endif

#if defined MATLAB_MEX_FILE
#include "mex.h"
#endif

#include "general_utils.h"
#include "file_utils.h"
using namespace std;

std::string commandlinestring(int argc, char** argv){
	std::string str = "Executing:";
	for (int i = 0; i < argc; i++){
		str += strprint(" %s", argv[i]);
	}
	return str;
};

std::string versionstring(const std::string& version, const std::string& compiletime, const std::string& compiledate){
	std::string s = strprint("Version: %s Compiled at %s on %s\n", version.c_str(), compiletime.c_str(), compiledate.c_str());
	return s;
};


bool mpi_initialised(){
#if defined _MPI_ENABLED
	int initialised;
	MPI_Initialized(&initialised);
	if (initialised)return true;
	else return false;
#else
	return false;
#endif
};

std::string mpi_processername(){
#if defined _MPI_ENABLED
	int len;
	char pname[MPI_MAX_PROCESSOR_NAME+1];
	MPI_Get_processor_name(pname, &len);
	pname[len] = '\0';
	return std::string(pname);
#else
	return std::string(" ");
#endif
};

int my_size(){
#if defined _MPI_ENABLED
	if (mpi_initialised()){
		int size;
		MPI_Comm_size(MPI_COMM_WORLD, &size);
		return size;
	}
	else return 1;
#else
	return 1;
#endif
};

int my_rank(){

	int threadnum = 0;
#if defined _OPENMP
	threadnum = omp_get_thread_num();
#endif

#if defined _MPI_ENABLED
	if (mpi_initialised()){
		int rank;
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);
		return rank;
	}
	else return threadnum;
#else
	return threadnum;
#endif
}

void my_barrier(){
#if defined _MPI_ENABLED
	if (mpi_initialised()){
		MPI_Barrier(MPI_COMM_WORLD);		
	}	
#endif
}

void rb_sleep(double secs)
{

#if defined _WIN32
	int ms = (int)(secs * 1000);
	Sleep((DWORD)ms);
#else        
	int s = (int)(secs);
	sleep(s);
#endif
}

void debug(const char* msg)
{
	message("Debug: %s\n", msg);
	fflush(stdout);
}

std::string strprint_va(const char* fmt, va_list vargs)
{
	size_t bufsize = 8192;
	std::string str;
	str.resize(bufsize);

	int status = vsnprintf(&(str[0]), bufsize, fmt, vargs);
	if (status < 0){
		str = std::string("");
		return str;
	}

	size_t len = (size_t)status;
	if (len < bufsize){
		str.resize(len);
	}
	else{
		bufsize = len + 1;
		str.resize(bufsize);
		status = vsnprintf(&(str[0]), bufsize, fmt, vargs);
		str.resize((size_t)status);
	}
	return str;
}

std::string strprint(const char* fmt, ...)
{
	va_list vargs;
	va_start(vargs, fmt);
	std::string s = strprint_va(fmt, vargs);
	va_end(vargs);
	return s;
}

void messageoutput(const std::string& msg)
{
	#if defined MATLAB_MEX_FILE
		mexPrintf("%s",msg.c_str());		
	#else			
		printf("%s", msg.c_str());
		fflush(stdout);
	#endif		
}

void messageoutput(const std::string& msg, FILE* fp)
{
	if (fp != NULL){
		fprintf(fp, "%s", msg.c_str());
		fflush(fp);
	}
}

void message(const char* fmt, ...)
{
	va_list vargs;
	va_start(vargs, fmt);
	std::string msg = strprint_va(fmt, vargs);
	va_end(vargs);
	messageoutput(msg);
}

void message(FILE* fp, const char* fmt, ...)
{
	va_list vargs;
	va_start(vargs, fmt);
	std::string msg = strprint_va(fmt, vargs);
	va_end(vargs);	
	messageoutput(msg);
	messageoutput(msg,fp);	
}

void rootmessage(const char* fmt, ...){
	va_list vargs;
	va_start(vargs, fmt);
	std::string msg = strprint_va(fmt, vargs);
	va_end(vargs);
	if (my_rank() == 0) messageoutput(msg);
}

void rootmessage(FILE* fp, const char* fmt, ...)
{
	va_list vargs;
	va_start(vargs, fmt);
	std::string msg = strprint_va(fmt, vargs);
	va_end(vargs);

	if (my_rank() == 0) messageoutput(msg);
	messageoutput(msg, fp);
}

void warningmessage(const char* fmt, ...)
{
	va_list vargs;
	va_start(vargs, fmt);
	std::string msg = "**Warning: " + strprint_va(fmt, vargs);
	va_end(vargs);

#if defined MATLAB_MEX_FILE
	mexWarnMsgTxt(msg.c_str());	
#else
	printf(msg.c_str());
	fflush(stdout);
#endif	
}

void warningmessage(FILE* fp, const char* fmt, ...)
{
	va_list vargs;
	va_start(vargs, fmt);
	std::string msg = "**Warning: " + strprint_va(fmt, vargs);
	va_end(vargs);

	fprintf(fp,msg.c_str());
	fflush(fp);

#if defined MATLAB_MEX_FILE
	mexWarnMsgTxt(msg.c_str());
#else
	printf(msg.c_str());
	fflush(stdout);
#endif	
}

void errormessage(const char* fmt, ...)
{
	va_list vargs;
	va_start(vargs, fmt);
	std::string msg = "**Error: " + strprint_va(fmt, vargs);
	va_end(vargs);
#if defined MATLAB_MEX_FILE
	mexErrMsgTxt(msg.c_str());				
#else
	printf(msg.c_str());
	fflush(stdout);
	throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
#endif		
}

void errormessage(FILE* fp, const char* fmt, ...){
	va_list vargs;
	va_start(vargs, fmt);
	std::string msg = "**Error: " + strprint_va(fmt, vargs);
	va_end(vargs);
	
	fprintf(fp,msg.c_str());
	fflush(fp);

#if defined MATLAB_MEX_FILE
	mexErrMsgTxt(msg.c_str());				
#else	
	printf(msg.c_str());
	fflush(stdout);
	throw(strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__));
#endif	

}

void prompttocontinue()
{
#if defined MATLAB_MEX_FILE				
	return;
#endif

	printf("Press any key to continue...\n");
#if defined _WIN32
	_getch();
#else
	char c;
	scanf("%c",&c);
#endif

}

void prompttoexit()
{
#if defined MATLAB_MEX_FILE		
	mexErrMsgTxt("Error");
	return;
#endif

#if defined _WIN32
	printf("Press any key to exit...\n");
	_getch();
	//#else
	//char c;
	//scanf("%c",&c);
#endif
	exit(0);
}

bool wildcmp(const char* wildpattern, const char* stringpattern)
{
	char* wild = new char[strlen(wildpattern) + 1];
	strcpy(wild, wildpattern);
	char* owild = wild;//oringinal pointer for delete

	char* str = new char[strlen(stringpattern) + 1];
	strcpy(str, stringpattern);
	char* ostring = str;//oringinal pointer for delete


	while ((*str) && (*wild != '*')) {
		if ((*wild != *str) && (*wild != '?')) {
			delete owild;
			delete ostring;
			return false;
		}
		wild++;
		str++;
	}

	char *cp = (char*)NULL;
	char *mp = (char*)NULL;
	while (*str) {
		if (*wild == '*') {
			if (!*++wild) {
				delete owild;
				delete ostring;
				return true;
			}
			mp = wild;
			cp = str + 1;
		}
		else if ((*wild == *str) || (*wild == '?')) {
			wild++;
			str++;
		}
		else {
			wild = mp;
			str = cp++;
		}
	}

	while (*wild == '*') {
		wild++;
	}
	bool rval = (bool)!*wild;
	delete owild;
	delete ostring;
	return rval;
}

double correlation_coefficient(std::vector<double>x, std::vector<double>y)
{
	double N = (double)x.size();
	double sumx = 0.0;
	double sumy = 0.0;
	double sumxx = 0.0;
	double sumyy = 0.0;
	double sumxy = 0.0;
	for (size_t i = 0; i < N; i++){
		sumx += x[i];
		sumy += y[i];
		sumxx += x[i] * x[i];
		sumyy += y[i] * y[i];
		sumxy += x[i] * y[i];
	}
	return (N*sumxy - sumx*sumy) / sqrt((N*sumxx - sumx*sumx)*(N*sumyy - sumy*sumy));

}

int log10stretch(double val, double lowclip, double highclip)
{
	if (val <= lowclip)return  0;
	if (val >= highclip)return 255;

	int bin;
	double logl = log10(lowclip);
	double logh = log10(highclip);
	if (val <= 0.0) bin = 0;
	else bin = (int)(256.0 * (log10(val) - logl) / (logh - logl));
	if (bin > 255)bin = 255;
	if (bin < 0)bin = 0;
	return bin;
}

double inverselog10stretch(int bin, double lowclip, double highclip)
{
	double logl = log10(lowclip);
	double logh = log10(highclip);
	double lval = ((double)bin * (logh - logl) / 255.0) + logl;
	return pow(10.0, lval);
}

int linearstretch(double val, double lowclip, double highclip)
{
	int bin;
	if (val <= lowclip) bin = 0;
	else if (val >= highclip) bin = 255;
	else bin = (int)(256.0 * (val - lowclip) / (highclip - lowclip));
	if (bin > 255)bin = 255;
	if (bin < 0)bin = 0;
	return bin;
}

double inverselinearstretch(int bin, double lowclip, double highclip){
	return ((double)bin * (highclip - lowclip) / 255.0) + lowclip;
}

bool regression(double* x, double*y, size_t n, double* gradient, double* intercept)
{
	//regression line y=mx + b // Thomas/Finney pp887
	//user must ensure line is not vertical	
	double sx = 0.0;
	double sy = 0.0;
	double sxy = 0.0;
	double sxx = 0.0;
	for (size_t i = 0; i < n; i++){
		sx += x[i];
		sy += y[i];
		sxy += x[i] * y[i];
		sxx += x[i] * x[i];
	}

	if ((sx*sx - n*sxx) == 0.0)return false;

	double m = (sx*sy - n*sxy) / (sx*sx - n*sxx);
	double b = (sy - m*sx) / (double)n;

	*gradient = m;
	*intercept = b;
	return true;
}

bool regression(const std::vector<double>& x, const std::vector<double>& y, double& gradient, double& intercept)
{
	//regression line y=mx + b // Thomas/Finney pp887
	//user must ensure line is not vertical
	double n = (double)x.size();
	if (n <= 0)return false;

	double sx = 0.0;
	double sy = 0.0;
	double sxy = 0.0;
	double sxx = 0.0;
	for (size_t i = 0; i < x.size(); i++){
		sx += x[i];
		sy += y[i];
		sxy += x[i] * y[i];
		sxx += x[i] * x[i];
	}

	if ((sx*sx - n*sxx) == 0.0)return false;

	gradient = (sx*sy - n*sxy) / (sx*sx - n*sxx);
	intercept = (sy - gradient*sx) / (double)n;

	return true;
}

bool bestfitlineendpoints(const std::vector<double>& x, const std::vector<double>& y, double& x1, double& y1, double& x2, double& y2)
{
	size_t n = x.size();
	double m = 0, c = 0;
	if (fabs(x[n - 1] - x[0]) > fabs(y[n - 1] - y[0])){
		//predominantly EW line
		regression(x, y, m, c);
		x1 = (x[0] + m*(y[0] - c)) / (m*m + 1);
		x2 = (x[n - 1] + m*(y[n - 1] - c)) / (m*m + 1);
		y1 = x1*m + c;
		y2 = x2*m + c;
	}
	else{
		//predominantly NS line
		regression(y, x, m, c);
		y1 = (y[0] + m*(x[0] - c)) / (m*m + 1);
		y2 = (y[n - 1] + m*(x[n - 1] - c)) / (m*m + 1);
		x1 = y1*m + c;
		x2 = y2*m + c;
	}
	return true;
}

std::string stringvalue(const double value, const char* fmt)
{
	if (value == -DBL_MAX)return std::string("Undefined");
	if (fmt == NULL) return strprint("%lf", value);
	return strprint(fmt, value);
}

std::string stringvalue(const size_t value, const char* fmt)
{	
	if (fmt == NULL) return strprint("%lu", value);
	return strprint(fmt, value);
}

std::string stringvalue(const int value, const char* fmt)
{	
	if (fmt == NULL) return strprint("%d", value);
	return strprint(fmt, value);
}

std::string stringvalue(const bool value)
{
	if (value == true)return std::string("True");
	return std::string("False");
}

int strcasecmp(const std::string& A, const std::string& B)
{
#if defined _MSC_VER //Microsoft Visual Studio compiler does not seem to define strcasecmp
	return _stricmp(A.c_str(), B.c_str());
#else
	return strcasecmp(A.c_str(), B.c_str());
#endif
}

int strncasecmp(const std::string& A, const std::string& B, const size_t n)
{
#if defined _MSC_VER //Microsoft Visual Studio compiler does not seem to define strncasecmp
	return _strnicmp(A.c_str(), B.c_str(), n);
#else
	return strncasecmp(A.c_str(), B.c_str(), n);
#endif
}

std::vector<std::string> fieldparsestring_old(const char* str, const char delim)
{
	int numquotes = 0;
	bool priorwhitespace = true;

	size_t len = 0;
	size_t k = 0;
	while (str[k] != '\0')
	{
		len++;
		if (str[k] == '"')numquotes++;
		k++;
	}

	if (numquotes % 2 != 0){
		warningmessage("fieldparsestd::string(): Unmatched quotes\n");
		std::vector<std::string> v;
		return v;
	}

	if (len == 0){
		warningmessage("fieldparsestd::string(): Zero length std::string\n");
		std::vector<std::string> v;
		return v;
	}

	std::vector<int> delimpos;
	delimpos.push_back(-1);
	for (size_t i = 0; i < len; i++){
		char c = str[i];

		if (c == ' ' || c == '\t')
		{
			if (priorwhitespace == false){
				delimpos.push_back((int)i);
			}
			priorwhitespace = true;
		}
		else if (c == delim)
		{
			delimpos.push_back((int)i);
			priorwhitespace = true;
		}
		else if (c == '\r' || c == '\n' || c == '\0')
		{
			delimpos.push_back((int)i);
			break;
		}
		else if (i == len - 1)
		{
			delimpos.push_back((int)i + 1);
			break;
		}
		else{
			priorwhitespace = false;
		}
	}

	size_t nfields = delimpos.size() - 1;
	std::vector<std::string> v(nfields);
	for (size_t i = 0; i < nfields; i++){
		size_t a = (size_t)(delimpos[i] + 1);
		size_t b = (size_t)(delimpos[i + 1] - 1);
		v[i] = std::string(&(str[a]), b - a + 1);
	}
	return v;
}

std::vector<std::string> fieldparsestring(const char* s, const char* delims)
{
	std::vector<std::string> fields;
	fields.reserve(400);
	std::string work = std::string(s);
	char* p = strtok(&(work[0]), delims);
	while (p != NULL)
	{
		fields.push_back(std::string(p));
		p = strtok(NULL, delims);
	}
	return fields;
}

std::vector<double> getdoublevector(const char* str, const char* delims)
{
	double v;
	std::vector<double> vec;
	std::vector<std::string> fields = fieldparsestring(str, delims);
	for (size_t i = 0; i < fields.size(); i++){
		sscanf(fields[i].c_str(), "%lf", &v);
		vec.push_back(v);
	}
	return vec;
}

std::string timestamp()
{	
	time_t ltime;
	time(&ltime);
	std::string str = std::string(ctime(&ltime));
	str.erase(str.length() - 1, 1);
	return str;
}

std::string timestring(const std::string format,  std::time_t t){
	
	if (t == 0) t = std::time(NULL);
	//"%Y%m%d";
	char str[100];	
	std::strftime(str, sizeof(str), format.c_str(), std::localtime(&t));
	return std::string(str);

}

int isinsidepolygon(int npol, double *xp, double *yp, double x, double y)
{
	//The following code is by Randolph Franklin, it returns 1 for interior points and 0 for exterior points. 
	int i, j, c = 0;
	for (i = 0, j = npol - 1; i < npol; j = i++) {
		if ((((yp[i] <= y) && (y < yp[j])) ||
			((yp[j] <= y) && (y < yp[i]))) &&
			(x < (xp[j] - xp[i]) * (y - yp[i]) / (yp[j] - yp[i]) + xp[i]))
			c = !c;
	}
	return c;
}

bool eq(double& a, double& b)
{
	if (fabs(a - b) < DBL_EPSILON)return true;
	else return false;
}

bool gt(double& a, double& b)
{
	if ((a - b) > DBL_EPSILON)return true;
	else return false;
}

bool lt(double& a, double& b)
{
	if ((b - a) > DBL_EPSILON)return true;
	else return false;
}

bool le(double& a, double& b)
{
	if (eq(a, b) || lt(a, b))return true;
	else return false;
}

bool ge(double& a, double& b)
{
	if (eq(a, b) || gt(a, b))return true;
	else return false;
}

void planeequation(const double& x1, const double& y1, const double& z1, const double& x2, const double& y2, const double& z2, const double& x3, const double& y3, const double& z3, double& A, double& B, double& C, double& D)
{
	//Ax+By+Cz+D=0
	double ax, ay, az, bx, by, bz;
	ax = x2 - x1;
	ay = y2 - y1;
	az = z2 - z1;

	bx = x3 - x1;
	by = y3 - y1;
	bz = z3 - z1;

	A = ay*bz - by*az;
	B = az*bx - bz*ax;
	C = ax*by - bx*ay;
	D = -(A*x1 + B*y1 + C*z1);
}

bool interplineline(const std::vector<double>& xin, const std::vector<double>& yin, std::vector<double>& xout, std::vector<double>& yout, double& dl)
{
	double d, x1, y1, x2, y2;
	double dx, dy, gradient, intercept;
	size_t n = xin.size();

	gradient = 0.0;
	intercept = 0.0;
	if (fabs(xin[n - 1] - xin[0]) > fabs(yin[n - 1] - yin[0])){
		regression(xin, yin, gradient, intercept);
		x1 = xin[0];
		y1 = x1*gradient + intercept;
		x2 = xin[n - 1];
		y2 = x2*gradient + intercept;

		d = sqrt(pow(2.0, (x2 - x1)) + pow(2.0, (y2 - y1)));

		size_t nl = (size_t)floor(0.5 + d / dl);
		xout.resize(nl);
		yout.resize(nl);

		dx = (x2 - x1) / (double)(nl - 1);
		for (size_t i = 0; i < nl; i++){
			xout[i] = x1 + (double)i*dx;
			yout[i] = xout[i] + gradient*intercept;
		}
	}
	else{
		regression(yin, xin, gradient, intercept);
		y1 = yin[0];
		x1 = y1*gradient + intercept;
		y2 = yin[n - 1];
		x2 = y2*gradient + intercept;

		d = sqrt(pow(2.0, (x2 - x1)) + pow(2.0, (y2 - y1)));

		size_t nl = (size_t)floor(0.5 + d / dl);
		xout.resize(nl);
		yout.resize(nl);

		dy = (y2 - y1) / (double)(nl - 1);
		for (size_t i = 0; i < nl; i++){
			yout[i] = y1 + (double)i*dy;
			xout[i] = yout[i] + gradient*intercept;
		}
	}

	return true;
}

int findindex(const size_t n, const double* x, const double& xtarget)
{
	int N = (int)n;

	if (xtarget < x[0])return -1;
	if (xtarget > x[n - 1])return N;

	int lo = 0;
	int hi = N - 1;

	while (lo <= hi){

		if (hi == lo + 1)return lo;

		int mid = (lo + hi) / 2;
		if (xtarget <= x[mid]){
			hi = mid;
		}
		else {
			lo = mid;
		}
	}
	return lo;
}

int findindex(const std::vector<double>& x, const double& xtarget)
{
	return findindex(x.size(), &(x[0]), xtarget);
}

double linearinterp(const size_t n, const double* x, const double* y, const double& xtarget)
{
	int k = findindex(n, x, xtarget);
	if (k < 0)k = 0;
	else if (k >= (int)(n)-1)k = (int)(n)-2;
	return linearinterp(x[k], y[k], x[k + 1], y[k + 1], xtarget);
}

double linearinterp(const double& x1, const double& y1, const double& x2, const double& y2, const double& xtarget)
{
	return ((y2 - y1) / (x2 - x1)) * (xtarget - x1) + y1;
}

double linearinterp(const std::vector<double>& x, const std::vector<double>& y, const double& xtarget)
{
	return linearinterp(x.size(), &(x[0]), &(y[0]), xtarget);
}

void   linearinterp(const size_t n, const double* x, const double* y, size_t ni, const double* xi, double* yi)
{
	for (size_t i = 0; i < ni; i++){
		yi[i] = linearinterp(n, x, y, xi[i]);
	}
}

std::vector<double> linearinterp(const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xi)
{
	std::vector<double> yi(xi.size());
	for (size_t i = 0; i < xi.size(); i++){
		yi[i] = linearinterp(x, y, xi[i]);
	}
	return yi;
}

std::vector<sRange> parserangelist(std::string& str)
{
	std::vector<std::string> items = parsestrings(str, ",");

	std::vector<sRange> list;
	sRange r;

	for (size_t i = 0; i < items.size(); i++){
		int f, t;
		int n = sscanf(items[i].c_str(), "%d to %d", &f, &t);
		//printf("%s\n",items[i].c_str());
		if (n == 1){
			r.from = f;
			r.to = f;
			list.push_back(r);
		}
		else if (n == 2){
			r.from = f;
			r.to = t;
			list.push_back(r);
		}
	}
	return list;
}

std::vector<std::string> parsestrings(const std::string& str, const std::string& delims)
{
	size_t len = strlen(str.c_str());
	char* workstr = new char[len + 1];
	strcpy(workstr, str.c_str());

	std::vector<std::string> list;

	char* token = strtok(workstr, delims.c_str());
	while (token != NULL){
		list.push_back(std::string(token));
		token = strtok(NULL, delims.c_str());
	}
	delete workstr;
	return list;
}

size_t bytesallocated(const std::vector<int>& v)
{
	return v.capacity()*sizeof(int);
}

size_t bytesallocated(const std::vector<double>& v)
{
	return v.capacity()*sizeof(double);
}

bool isreportable(int rec)
{
	int k, interval;

	if (rec < 1000)k = 2;
	else        k = (int)log10((double)rec);

	interval = (int)pow(10.0, (double)k);

	if (rec%interval == 0) return true;
	else return false;
}

bool isinrange(const sRange& r, const int& i)
{
	if (i<r.from || i>r.to)return false;
	else return true;
}

double overlap(const double al, const double ah, const double bl, const double bh)
{
	//this returns the amount of overlap between the ranges al-ah with bl-bh
	if (al >= bh || ah <= bl)return 0;//no overlap 
	else if (al <= bl){
		if (ah >= bh)return bh - bl;//a overlaps all of b
		else return ah - bl;//a overlaps low end of b 
	}
	else{
		if (ah <= bh)return ah - al;//a is all inside b
		else return bh - al;//a overlaps high end of b
	}
	warningmessage("overlap(): Found undefined condition al=%lf ah=%lf bl=%lf bh=%lf\n", al, ah, bl, bh);
	return 0;
}

double fractionaloverlap(const double al, const double ah, const double bl, const double bh)
{
	double o = overlap(al, ah, bl, bh);
	return o / (ah - al);
}

std::vector<double> overlaps(const double& a1, const double& a2, const std::vector<double>& b)
{
	std::vector<double> o(b.size() - 1);	
	o.resize(b.size() - 1);
	for (size_t bi = 0; bi < b.size() - 1; bi++){
		o[bi] = overlap(a1,a2, b[bi], b[bi + 1]);
	}	
	return o;
}

std::vector<double> fractionaloverlaps(const double& a1, const double& a2, const std::vector<double>& b)
{
	std::vector<double> o(b.size() - 1);
	o.resize(b.size() - 1);
	for (size_t bi = 0; bi < b.size() - 1; bi++){
		o[bi] = fractionaloverlap(a1, a2, b[bi], b[bi + 1]);
	}
	return o;
}

std::vector<std::vector<double>> overlaps(const std::vector<double>& a, const std::vector<double>& b)
{
	std::vector<std::vector<double>> o(a.size() - 1);
	for(size_t ai = 0; ai < a.size()-1; ai++){
		o[ai].resize(b.size() - 1);
		for (size_t bi = 0; bi < b.size()-1; bi++){
			o[ai][bi] = overlap(a[ai], a[ai + 1], b[bi], b[bi + 1]);
		}	
	}
	return o;
}
std::vector<std::vector<double>> fractionaloverlaps(const std::vector<double>& a, const std::vector<double>& b)
{
	std::vector<std::vector<double>> o(a.size() - 1);
	for (size_t ai = 0; ai < a.size() - 1; ai++){
		o[ai].resize(b.size() - 1);
		for (size_t bi = 0; bi < b.size() - 1; bi++){
			o[ai][bi] = overlap(a[ai], a[ai + 1], b[bi], b[bi + 1]);
		}
	}
	return o;
}

double gettime(){
	struct timeb t;
	ftime(&t);
	return (double)t.time + 0.001*(double)t.millitm;
}

char* temppath(const char* s, int set)
{
	static char* p = (char*)NULL;
	if (set == 1){
		if (p)delete[]p;
		p = new char[strlen(s) + 1];
		strcpy(p, s);
		return (char*)NULL;
	}
	else return p;
}

void settemppath(const char* s)
{
	temppath(s, 1);
}

std::string gettemppath()
{
	char* s = (char*)NULL;
	char* p = temppath(s, 0);
	return std::string(p);
}

int floatcompare(const void *pa, const void *pb)
{
	float& a = *(float*)pa; float& b = *(float*)pb;
	if (a < b)return -1;
	else if (a == b)return 0;
	else return 1;
}

void sort(float* x, const size_t n)
{
	qsort(x, n, sizeof(float), floatcompare);
}

int doublecompare(const void *pa, const void *pb)
{
	double& a = *(double*)pa; double& b = *(double*)pb;
	if (a < b)return -1;
	else if (a == b)return 0;
	else return 1;
}

void sort(double* x, const size_t n)
{
	qsort(x, n, sizeof(double), doublecompare);
}

int stringcompare(const void *pa, const void *pb)
{
	char** a = (char**)pa;
	char** b = (char**)pb;
	return strcmp(*a, *b);
}

void sort(char** strings, const size_t n)
{
	qsort(strings, n, sizeof(char*), stringcompare);
}

int intcompare(const void *pa, const void *pb)
{
	int& a = *(int*)pa; int& b = *(int*)pb;
	return a - b;
}

void sort(int* x, const size_t n)
{
	qsort(x, n, sizeof(int), intcompare);
}

#if defined _WIN32
double reportusage()
{
	MEMORYSTATUS memstatus;
	GlobalMemoryStatus(&memstatus);
	return memstatus.dwMemoryLoad;
}
#else
double reportusage()
{	
	int pid = getpid();	
	double vsize,pcpu,pmem;	
	std::string tmpfile = strprint("ps.%d.tmp",pid);
	std::string cmd     = strprint("ps --pid %d --format pcpu,vsize,pmem > %s\n",pid,tmpfile.c_str());
	system(cmd.c_str());
	FILE* fp=fileopen(tmpfile,"r");
	char buf[201];
	fgets(buf,200,fp);
	fgets(buf,200,fp);
	sscanf(buf,"%lf %lf %lf",&pcpu,&vsize,&pmem);
	fclose(fp);
	deletefile(tmpfile);
	message("Percent CPU used: %.2lf\n",pcpu);
	message("Percent memory used: %.2lf\n",pmem);
	message("Virtual memory used (Mb): %.2lf\n",vsize/1000.0);	
	return pmem;
}
#endif

#if defined _WIN32
double percentmemoryused()
{
	MEMORYSTATUS memstatus;
	GlobalMemoryStatus(&memstatus);
	return memstatus.dwMemoryLoad;
}
#else
double percentmemoryused()
{	
	int pid = getpid();	
	double vsize,pcpu,pmem;	
	std::string tmpfile = strprint("ps.%d.tmp",pid);
	std::string cmd = strprint("ps --pid %d --format pcpu,vsize,pmem > %s\n",pid,tmpfile.c_str());
	system(cmd.c_str());
	FILE* fp=fileopen(tmpfile.c_str(),"r");
	char buf[201];
	fgets(buf,200,fp);
	fgets(buf,200,fp);
	sscanf(buf,"%lf %lf %lf",&pcpu,&vsize,&pmem);
	fclose(fp);
	deletefile(tmpfile.c_str());	
	return pmem;
}
#endif

void guage(int ntot, int n, int pdiv1, int pdiv2)
{
	double d1 = ceil(((double)pdiv1 / 100.0)*(double)ntot);
	double d2 = ceil(((double)pdiv2 / 100.0)*d1);

	if (n == 0){
		printf("<");
	}
	else if (n >= ntot - 1){
		printf(">");
	}
	else if (n % (int)d1 == 0){
		printf("+");
	}
	else if (n % (int)d2 == 0){
		printf(".");
	}
}

void allocate1darray(int*& a, const size_t n)
{
	a = new int[n];
	if (a == (int*)NULL){
		printf("allocate1darray(int*& a, int n) returned NULL\n");
		prompttocontinue();
	}
}

void allocate1darray(double*& a, const size_t n)
{
	a = new double[n];
	if (a == (double*)NULL){
		printf("allocate1darray(double*& a, int n) returned NULL\n");
		prompttocontinue();
	}
}

void allocate2darray(double**& a, const size_t nrows, const size_t ncols)
{
	a = new double*[nrows];
	for (size_t i = 0; i < nrows; i++)a[i] = new double[ncols];
}

void deallocate2darray(double**& a, const size_t nrows)
{
	if (a){
		for (size_t i = 0; i < nrows; i++){
			if (a[i])delete[]a[i];
		}
		delete[]a;
		a = (double**)NULL;
	}
}

void deallocate1darray(double*& a)
{
	if (a){
		delete[]a;
		a = (double*)NULL;
	}
}

void deallocate1darray(int*& a)
{
	if (a){
		delete[]a;
		a = (int*)NULL;
	}
}

double median(double* v, const size_t n)
{
	double* d = new double[n];
	for (size_t i = 0; i < n; i++)d[i] = v[i];

	sort(d, n);

	double m = d[n / 2];

	delete[]d;
	return m;
}

std::string trim(const std::string& s)
{
	if (s.length() == 0) return s;
	size_t index1 = s.find_first_not_of(" \t\r\n");
	size_t index2 = s.find_last_not_of(" \t\r\n");
	if (index1 == std::string::npos || index2 == std::string::npos){
		return std::string("");
	}
	else{
		return s.substr(index1, index2 - index1 + 1);
	}
}

std::string stripquotes(const std::string& s)
{
	if (s.length() == 0) return s;
	if (s[0] == '"' && s[s.size() - 1] == '"'){
		return s.substr(1,s.size()-2);
	}	
	return s;
}

std::vector<std::string> tokenize(const std::string& str)
{
	std::stringstream strstr(str);
	istream_iterator<std::string> it(strstr);
	istream_iterator<std::string> end;
	std::vector<std::string> results(it, end);
	return results;
}

std::vector<float> dvec2fvec(std::vector<double>& vd)
{
	std::vector<float> vf(vd.size());
	for (size_t i = 0; i < vd.size(); i++)vf[i] = (float)vd[i];
	return vf;
}

std::vector<double> linspace(const double x1, const double x2, const size_t n)
{
	double dx = (x2 - x1) / (double)(n - 1);
	std::vector<double> xi(n);
	xi[0] = x1;
	for (size_t i = 1; i < n; i++){
		xi[i] += xi[i - 1] + dx;
	}
	return xi;
};

std::vector<double> log10space(const double x1, const double x2, const size_t n)
{
	vector<double> v = linspace(std::log10(x1), std::log10(x2), n);
	for (size_t i = 0; i < n; i++){
		v[i] = std::pow(10.0,v[i]);
	}
	return v;
};

bool filegetline(FILE* fp, std::string& str)
{
	size_t buflen = 8192;//buffer length
	str.clear();
	std::vector<char> buf(buflen + 1);
	while (fgets(&(buf[0]), (int)buflen, fp) != NULL){
		str += std::string(&(buf[0]));
		size_t len = strlen(&(buf[0]));
		if (len < buflen - 1){
			if (str[str.length() - 1] == 10){
				//strip linefeed character
				str.resize(str.length() - 1);
			}
			return true;
		}
	}
	return false;
}

unsigned int factorial(unsigned int n){
	if (n <= 1) return 1;
	unsigned int fact = 1;
	for (unsigned int i = 2; i <= n; i++){
		fact *= i;
	}
	return fact;
}

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

int LevenshteinDistance(char* s, int len_s, char* t, int len_t)
{
	//Adapted from https://en.wikipedia.org/wiki/Levenshtein_distance#Example
	int cost;
	/* base case: empty strings */
	if (len_s == 0) return len_t;
	if (len_t == 0) return len_s;
	/* test if last characters of the strings match */
	if (s[len_s - 1] == t[len_t - 1])cost = 0;
	else cost = 1;

	/* return minimum of delete char from s, delete char from t, and delete char from both */
	//return std::min(LevenshteinDistance(s, len_s - 1, t, len_t) + 1, LevenshteinDistance(s, len_s, t, len_t - 1) + 1, LevenshteinDistance(s, len_s - 1, t, len_t - 1) + cost);

	int a = LevenshteinDistance(s, len_s - 1, t, len_t) + 1;
	int b = LevenshteinDistance(s, len_s, t, len_t - 1) + 1;
	int c = LevenshteinDistance(s, len_s - 1, t, len_t - 1) + cost;
	if (a <= b && a <= c)return a;
	else if (b <= c)return b;
	return c;
}