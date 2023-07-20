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

using DoublePair = std::pair<double, double>;

class cInversionLineSearcher{

	size_t maxtrials = 20;	
	double target;
	std::vector<std::pair<double, double>> trials;
	double xtol = 1e-8;
	double ytol =  0.0;	

public:

	cInversionLineSearcher(const double& _target){
		target  = _target;
		ytol    = target*0.01;		
	};

	~cInversionLineSearcher() {
		//std::cerr << *this;
		//std::cout << *this;
	};

	void set_ytol(const double& _ytol) {
		ytol = _ytol;
	}

	void set_maxtrials(const size_t& _maxtrials) {
		maxtrials = _maxtrials;
	}

	DoublePair get_pair(const size_t& i) {
		return trials[i];
	};

	bool xtol_ok(const double& x) const {
		for (size_t i = 0; i < trials.size(); i++) {
			if (std::fabs(x - trials[i].first) < xtol) {
				return false;
			}
		}
		return true;
	}

	bool add_pair(const double& x, const double& y) {
		if (xtol_ok(x)) {
			trials.push_back(DoublePair(x, y));
			std::sort(trials.begin(), trials.end());
			return true;
		}
		return false;
	};
	
	bool next_x(double& nextx) const {
		
		if (trials.size() >= maxtrials)return false;

		nextx = 0.0;
		DoublePair p = nearest();
		if (std::fabs(p.second - target) <= ytol) {
			return false;
		}

		if (trials.size() == 1) {
			nextx = 0.01;
			return xtol_ok(nextx);
		}
		else if (trials.size() == 2) {
			double y0 = trials[0].second;
			double y1 = trials[1].second;
			if (y1 < y0) {
				nextx = linear_estimate(0, 1);
				return xtol_ok(nextx);				
			}
			return false;
		}

		size_t ibrak;
		if (is_target_braketed(ibrak)) {
			nextx = linear_estimate(ibrak, ibrak + 1);
			return xtol_ok(nextx);
		}

		size_t imin;
		double xmin, ymin;
		if (is_min_braketed(imin, xmin, ymin)) {
			nextx = quadratic_estimate(imin - 1, imin, imin + 1);
			return xtol_ok(nextx);
		}

		if (ymin >= target) {
			if (imin == 0) {	
				const double& x0 = trials[0].first;
				const double& x1 = trials[1].first;
				const double dx = x1-x0;
				nextx = x0 - dx / 2.0;
				return xtol_ok(nextx);
			}
			else if (imin == trials.size() - 1) {
				const double& x0 = trials[trials.size() - 2].first;
				const double& x1 = trials[trials.size() - 1].first;								
				const double dx = x1 - x0;
				nextx = x1 + dx / 2.0;
				return xtol_ok(nextx);
			}
		}
		else {			
			double xmax, ymax;
			size_t imax = max_index(xmax, ymax);
			if (imax == 0) {
				const double dx = trials[1].first - trials[0].first;
				nextx = trials[0].first - dx / 2.0;				
				return xtol_ok(nextx);
			}
			else if (imax == trials.size() - 1) {
				const double dx = trials[trials.size() - 1].first - trials[trials.size() - 2].first;
				nextx = trials[trials.size() - 1].first + dx / 2.0;				
				return xtol_ok(nextx);
			}
			else {	
				//Must be convex up 
				std::cerr << "Warning unexpected circumstance" << std::endl;;
				std::cout << "Warning unexpected circumstance" << std::endl;;
			    //char c; std::cin >> c;
				return false;
			}
		}
		return false;
	};
	
	double linear_estimate(const size_t i0, const size_t i1) const {
		double x0 = trials[i0].first;
		double y0 = trials[i0].second;
		double x1 = trials[i1].first;
		double y1 = trials[i1].second;
		if (x1==x0) {
			glog.errormsg(_SRC_ + "\nDivide by zero error 'a' in linear estimate");
		}
		double m = (y1 - y0) / (x1 - x0);
		double c = y0 - m*x0;
		double x3 = (target - c) / m;
		return x3;
	}

	double quadratic_estimate(const size_t i1, const size_t i2, const size_t i3) const {
		const double& x1 = trials[i1].first;
		const double& x2 = trials[i2].first;
		const double& x3 = trials[i3].first;
		const double& y1 = trials[i1].second;
		const double& y2 = trials[i2].second;
		const double& y3 = trials[i3].second;

		double denom = (x1 - x2) * (x1 - x3) * (x2 - x3);
		if (denom == 0) {
			glog.errormsg(_SRC_ + "\nDivide by zero error 'denom' in quadratic estimate");
		}
		double a = (x3 * (y2 - y1) + x2 * (y1 - y3) + x1 * (y3 - y2)) / denom;
		double b = (x3*x3 * (y1 - y2) + x2*x2 * (y3 - y1) + x1*x1 * (y2 - y3)) / denom;
		//double c = (x2 * x3 * (x2 - x3) * y1 + x3 * x1 * (x3 - x1) * y2 + x1 * x2 * (x1 - x2) * y3) / denom;
		double xmin = -b / (2.0 * a);
		if (a == 0) {
			glog.errormsg(_SRC_ + "\nDivide by zero error 'a' in quadratic estimate");
		}
		//double ymin = c - b*b / (4 * a);
		return xmin;
	}

	bool is_target_braketed(size_t& ibrak)const {
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

	bool is_min_braketed(size_t& imin, double& xmin, double& ymin)const {
		imin = min_index(xmin, ymin);
		if (imin > 0 && imin < trials.size() - 1){
			return true;
		}
		return false;
	}

	size_t nearest_index() const {
		size_t index = 0;
		double mind = std::fabs(trials[0].second - target);
		double xnearest = trials[0].first;
		double ynearest = trials[0].second;
		for (size_t k = 1; k < trials.size(); k++) {
			double d = std::fabs(trials[k].second - target);
			if (d < mind) {
				mind = d;
				index = k;
				xnearest = trials[k].first;
				ynearest = trials[k].second;
			}
		}
		return index;
	};

	DoublePair nearest() const {
		return trials[nearest_index()];		
	};

	size_t min_index(double& xmin, double& ymin) const {
		size_t index = 0;
		xmin = trials[0].first;
		ymin = trials[0].second;
		for (size_t k = 1; k < trials.size(); k++){
			if (trials[k].second < ymin){
				xmin = trials[k].first;
				ymin = trials[k].second;
				index = k;
			}			
		}
		return index;
	};

	size_t max_index(double& xmax, double& ymax) const {
		size_t index = 0;
		xmax = trials[0].first;
		ymax = trials[0].second;
		for (size_t k = 1; k < trials.size(); k++) {
			if (trials[k].second > ymax) {
				xmax = trials[k].first;
				ymax = trials[k].second;
				index = k;
			}
		}
		return index;
	};

	friend std::ostream& operator<<(std::ostream& os, const cInversionLineSearcher& s) {		
		os << "*Target=" << s.target << std::endl;
		for (size_t k = 0; k < s.trials.size(); k++) {
			os << "Frac=" << std::setw(10) << s.trials[k].first << "\tPhiD=" << std::setw(10) << s.trials[k].second << std::endl;
		}
		return os;
	}
};


#endif

