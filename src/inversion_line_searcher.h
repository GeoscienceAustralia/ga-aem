/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _inversion_line_searcher_H
#define _inversion_line_searcher_H

#include <algorithm>
#include <vector>
#include <file_utils.h>

class cInversionLineSearcher{

	double current;
	double target;
	std::vector<std::pair<double, double>> trials;
	double ytol;

public:

	cInversionLineSearcher(double _current, double _target){
		current = _current;
		target = _target;
		ytol = target*0.05;
		addtrial(0.0, current);
	};

	bool next(double& nextx){
		nextx = 0.0;
		if (trials.size() >= 10)return false;

		double xnearest, ynearest;
		nearestindex(xnearest, ynearest);
		if (std::fabs(ynearest - target) <= ytol){
			return false;
		}

		if (trials.size() == 1){
			//if (cMpiEnv::world_rank() == 0)printf("1st it\n");
			nextx = 0.01;
			return true;
		}

		if (trials.size() == 2){
			double y0 = trials[0].second;
			double y1 = trials[1].second;
			if (y1 < y0){
				//if (cMpiEnv::world_rank() == 0)printf("2nd it linear\n");
				nextx = linear_estimate(0, 1);
				return true;
			}
			//if (cMpiEnv::world_rank() == 0)printf("2nd it fail\n");
			return false;
		}

		size_t ibrak;
		if (istargetbraketed(ibrak)){
			//if (cMpiEnv::world_rank() == 0)printf("tarbrak\n");
			nextx = linear_estimate(ibrak, ibrak + 1);
			return true;
		}

		size_t imin;
		double xmin, ymin;
		if (isminbraketed(imin, xmin, ymin)){
			//if (cMpiEnv::world_rank() == 0)printf("minbrak qest %lf\n", trials[imin].first);
			nextx = quadratic_estimate(imin - 1, imin, imin + 1);
			return true;

			double dx0 = trials[imin].first - trials[imin - 1].first;
			double dx1 = trials[imin + 1].first - trials[imin].first;
			if (dx0 > dx1){
				//if (cMpiEnv::world_rank() == 0)printf("minbrak left %lf\n", trials[imin].first);
				nextx = trials[imin].first - dx0 / 2.0;
			}
			else{
				//if (cMpiEnv::world_rank() == 0)printf("minbrak right %lf\n", trials[imin].first);
				nextx = trials[imin].first + dx1 / 2.0;
			}
			return true;
		}

		if (ymin>target){
			if (imin == 0){
				//if (cMpiEnv::world_rank() == 0)printf("imin=low\n");
				nextx = trials[1].first / 2.0;
				return true;
			}

			if (imin == trials.size() - 1){
				//if (cMpiEnv::world_rank() == 0)printf("imin=high\n");
				nextx = trials[trials.size() - 1].first * 2.0;
				return true;
			}
		}
		printf("Error unexpected circumstance\n");
		return false;
	};

	void addtrial(double x, double y){
		std::pair<double, double> p(x, y);
		trials.push_back(p);
		std::sort(trials.begin(), trials.end());
	};

	double linear_estimate(const size_t i0, const size_t i1){
		double x0 = trials[i0].first;
		double y0 = trials[i0].second;
		double x1 = trials[i1].first;
		double y1 = trials[i1].second;
		double m = (y1 - y0) / (x1 - x0);
		double c = y0 - m*x0;
		double x3 = (target - c) / m;
		return x3;
	}

	double quadratic_estimate(const size_t i1, const size_t i2, const size_t i3){
		const double& x1 = trials[i1].first;
		const double& x2 = trials[i2].first;
		const double& x3 = trials[i3].first;
		const double& y1 = trials[i1].second;
		const double& y2 = trials[i2].second;
		const double& y3 = trials[i3].second;

		double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
		double a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
		double b = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
		//double c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
		double xmin = -b / (2.0 * a);
		//double ymin = c - b*b / (4 * a);
		return xmin;
	}

	bool istargetbraketed(size_t& ibrak){
		for (size_t k = 0; k < trials.size() - 1; k++){
			if (trials[k].second > target && trials[k + 1].second < target){
				ibrak = k;
				return true;
			}
			else if (trials[k].second < target && trials[k + 1].second > target){
				ibrak = k;
				return true;
			}
		}
		return false;
	}

	bool isminbraketed(size_t& imin, double& xmin, double& ymin){
		imin = minindex(xmin, ymin);
		if (imin > 0 && imin < trials.size() - 1){
			return true;
		}
		return false;
	}

	size_t nearestindex(double& xnearest, double& ynearest){
		size_t index = 0;
		double mind = std::fabs(trials[0].second - target);
		xnearest = trials[0].first;
		ynearest = trials[0].second;
		for (size_t k = 1; k < trials.size(); k++){
			double d = std::fabs(trials[k].second - target);
			if (d < mind){
				mind = d;
				index = k;
				xnearest = trials[k].first;
				ynearest = trials[k].second;
			}
		}
		return index;
	};

	size_t minindex(double& xmin, double& ymin){
		size_t index = 0;
		xmin = trials[0].first;
		ymin = trials[0].second;
		for (size_t k = 1; k < trials.size(); k++){
			if (trials[k].second < ymin){
				xmin = trials[k].first;
				ymin = trials[k].second;
				index = k;
			}
			else{
				break;
			}
		}
		return index;
	};

	void printtrials(){
		for (size_t k = 0; k < trials.size(); k++){
			printf("sf=%lf phid=%lf\n", trials[k].first, trials[k].second);
		}
	}

	void writetextfile(const std::string& name){
		FILE* fp = fileopen(name.c_str(), "w");
		for (size_t k = 0; k < trials.size(); k++){
			fprintf(fp, "%lf %lf\n", trials[k].first, trials[k].second);
		}
		fclose(fp);
	}

};


#endif

