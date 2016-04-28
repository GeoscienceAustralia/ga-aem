/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _general_types_H
#define _general_types_H

#include <climits>
#include <cstring>
#include <complex>
#include <vector>

typedef std::vector<double>  dvector;
typedef std::vector<dvector> dmatrix;

typedef std::complex<double> cdouble;
typedef std::vector<cdouble> cvector;
typedef std::vector<cvector> cmatrix;

template<typename T>
class cHistogramStats{

public:
	size_t nbins;//number of bins in histogram
	size_t nsamples;//number of samples in histogram
	T min;//minimum of sample
	T max;//maximum of sample
	T mean;
	T std;
	T var;	
	T mode;
	T p10;
	T p50;
	T p90;

	cHistogramStats(){
		nbins = 0;
	}

	template<typename U> 
	cHistogramStats(const std::vector<T>& bins, const std::vector<U>& counts){
		compute(bins, counts);
	};

	template<typename U>
	void compute(const std::vector<T>& bins, const std::vector<U>& counts){
		
		nbins = bins.size();

		min = bins[nbins-1];
		max = bins[0];

		T sum = 0.0;
		nsamples = 0;
		std::vector<size_t> cumcounts(nbins + 1, 0);
		for (size_t i = 0; i < nbins; i++){
			U c = counts[i];
			nsamples += c;
			cumcounts[i + 1] = cumcounts[i] + c;

			T s = bins[i] * (T)c;
			sum += s;

			if (c>0){
				if (bins[i] < min)min = bins[i];
				if (bins[i] > max)max = bins[i];
			}
		}
		mean = sum / (T)nsamples;
		T np10 = 0.1*(T)nsamples;
		T np50 = 0.5*(T)nsamples;
		T np90 = 0.9*(T)nsamples;

		size_t modebin = 0;
		T sumdsqr = 0.0;
		for (size_t i = 0; i < nbins; i++){
			T d = bins[i] - mean;
			sumdsqr += d*d*(T)counts[i];

			if (cumcounts[i] <= np10 && cumcounts[i + 1] >= np10)p10 = bins[i];
			if (cumcounts[i] <= np50 && cumcounts[i + 1] >= np50)p50 = bins[i];
			if (cumcounts[i] <= np90 && cumcounts[i + 1] >= np90)p90 = bins[i];

			if (counts[i] > counts[modebin]){
				mode = bins[i];
				modebin = i;
			}
		}
		var = sumdsqr / (T)nsamples;
		std = sqrt(var);		
	}
};

template<typename T, typename U>
class cHistogram{

public:
	size_t nbins;
	std::vector<T> edge;
	std::vector<T> centre;
	std::vector<U> count;

	cHistogram(){
		nbins = 0;
	}

	cHistogram(const std::vector<T>& v, T hmin, T hmax, size_t _nbins){
		compute(v, hmin, hmax, _nbins);
	}


	void compute(const std::vector<T>& v, T hmin, T hmax, size_t _nbins)
	{
		nbins = _nbins;
		T dx = (hmax - hmin) / (nbins);
		T e = hmin;
		edge.push_back(e);
		for (size_t i = 0; i < nbins; i++){
			centre.push_back(e + dx / 2.0);
			e += dx;
			edge.push_back(e);
		}
		count.resize(nbins, 0);

		for (size_t i = 0; i < v.size(); i++){
			if (v[i] < hmin)count[0]++;
			if (v[i] > hmax)count[nbins - 1]++;
			size_t b = (size_t)floor((v[i] - hmin) / dx);
			count[b]++;
		}		
	}
};


template<typename T>
class cStats{	

public:
	size_t nulls;
	size_t nonnulls;
	T min;
	T max;
	T mean;
	T var;
	T std;		

	cStats(){
		nulls = 0;
		nonnulls = 0;
	}

	cStats(const std::vector<T>& v){
		compute(v);
	}

	cStats(const std::vector<T>& v, T nullvalue){
		compute_with_nulls(v,nullvalue);
	}

	void compute(const std::vector<T>& v)
	{
		nulls    = 0;
		nonnulls = 0;
		min = v[0];
		max = v[0];

		double sx  = 0.0;
		double sx2 = 0.0;
		for (size_t i = 0; i < v.size(); i++){
			nonnulls++;

			if (v[i] < min) min = v[i];
			else if (v[i] > max) max = v[i];

			sx  += v[i];
			sx2 += (v[i] * v[i]);
		};		
		mean = sx / nonnulls;
		var  = (sx2 - (sx*sx) / nonnulls) / (nonnulls - 1.0);
		std  = sqrt(var);
	}

	void compute_with_nulls(const std::vector<T>& v, const T nullvalue)
	{
		nulls    = 0;
		nonnulls = 0;

		double sx  = 0.0;
		double sx2 = 0.0;
		for (size_t i = 0; i < v.size(); i++){
			if (v[i] == nullvalue){
				nulls++;
				continue;
			}
			nonnulls++;

			if(nonnulls==1){
				min = v[i];
				max = v[i];
			}
			else if (v[i] < min) min = v[i];
			else if (v[i] > max) max = v[i];

			sx  += v[i];
			sx2 += (v[i] * v[i]);
		};		
		mean = sx / nonnulls;
		var  = (sx2 - (sx*sx) / nonnulls) / (nonnulls - 1.0);
		std  = sqrt(var);
	}
};

struct sRange{
	int from;
	int to;
};

struct sBoundingBox{
	double xlow;
	double xhigh;
	double ylow;
	double yhigh;
};

class cPoint{

public:
	double x;
	double y;
	cPoint(){};
	cPoint(const double& px, const double& py){
		x = px;
		y = py;
	}

	bool operator==(const cPoint& p){
		if (x == p.x && y == p.y)return true;
		return false;
	}
	
};

class cEarth1D{ 
	public:
	std::vector<double> conductivity;	
	std::vector<double> thickness;
	size_t nlayers(){return conductivity.size();}
	void print(){
		for(size_t i=0; i<nlayers()-1; i++){
			printf("%d\t%8.6lf\t%6.2lf\n",(int)i,conductivity[i],thickness[i]);
		}
		printf("%d\t%8.6lf\n\n",(int)nlayers(),conductivity[nlayers()-1]);
	}
};

struct sAirborneSampleId{
  size_t uniqueid;
  size_t surveynumber;
  size_t daynumber;
  size_t flightnumber;
  size_t linenumber;
  double fidnumber;  
};

struct sAirborneSampleLocation{  
	double x;
	double y;  	
	double z;	
	double groundelevation;	
};


template <typename T>
class c3DArray{

private:
	std::vector<std::vector<std::vector<T>>> data;

public:
	c3DArray<T>(int ni = 0, int nj = 0, int nk = 0){ resize(ni, nj, nk); }
	void resize(int ni, int nj, int nk){
		data.resize(ni);
		for (int i = 0; i < ni; ++i){
			data[i].resize(nj);
			for (int j = 0; j < nj; ++j){
				data[i][j].resize(nk);								
			}
		}
	}

	void initialise(const T& v){		
		for (int i = 0; i < ni(); ++i){			
			for (int j = 0; j < nj(); ++j){				
				for (int k = 0; k < nk(); ++k){
					data[i][j][k] = v;
				}
			}
		}
	}

	std::vector<std::vector<T>>& operator[](int index) {
		return data[index];
	}

	c3DArray<T>& operator=(const T& v) {
		initialise(v);
		return *this;
	}

	int ni() const { return (int) data.size(); }
	int nj() const { return (int) data[0].size(); }
	int nk() const { return (int) data[0][0].size(); }
};

#endif
