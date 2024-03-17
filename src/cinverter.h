/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _cinverter_H
#define _cinverter_H

#include <stdio.h>
#include <sstream>
#include <vector>
#include <cstring>
#include <algorithm>

#include <Eigen/Dense>
#include "general_constants.h"
#include "general_types.h"
#include "inputmanager.h"
#include "outputmanager.h"
#include "blocklanguage.h"
#include "eigen_utils.h"
#include "gaaem_version.h"
#include "samplebunch.h"
#include "inversion_line_searcher.h"

#if defined _OPENMP
#include <omp.h>
//This external thread lock must be set when fftw is being initialised
extern omp_lock_t fftw_thread_lock;
#endif

enum class eBracketResult { BRACKETED, MINBRACKETED, ALLABOVE, ALLBELOW };
enum class eNormType { L1, L2 };

class cInvertibleFieldDefinition {

public:
	int  offset = -1;
	bool solve = false;
	cFieldDefinition input;
	cFieldDefinition ref;
	cFieldDefinition std;
	cFieldDefinition min;
	cFieldDefinition max;
	cFieldDefinition tfr;

	cInvertibleFieldDefinition() {};

	cInvertibleFieldDefinition(const cBlock& parent, const std::string& key) {
		initialise(parent, key);
	};

	bool initialise(const cBlock& parent, const std::string& key) {
		std::string id = parent.findkey(key);
		if (id.compare(undefinedvalue<std::string>()) != 0) {
			return entryinit(parent, key);
		}
		else {
			cBlock b = parent.findblock(key);
			if (b.empty() == true) {
				std::string msg;
				msg += strprint("Could not find control file block: %s\n", key.c_str());
				glog.errormsg(msg);
			}

			if (blockinit(b) == false) {
				std::string msg;
				msg += strprint("Could not parse control file block: %s\n", key.c_str());
				std::cerr << msg;
				glog.errormsg(msg);
			}
			return true;
		}
	};

	bool entryinit(const cBlock& b, const std::string& key) {
		offset = -1;
		input.initialise(b, key);
		return true;
	};

	bool blockinit(const cBlock& b) {
		b.get("solve", solve, false);
		input.initialise(b, "input");
		ref.initialise(b, "ref");
		std.initialise(b, "std");
		min.initialise(b, "min");
		max.initialise(b, "max");
		tfr.initialise(b, "tfr");
		return true;
	};

	bool bound() const {

		if (solve && min.isinitialised() && max.isinitialised()) {
			return true;
		}
		else return false;
	}
};

class cTrial {

public:
	size_t insert_order = 0;
	double lambda = 0.0;
	double stepfactor = 0.0;
	double phid = 0.0;

	static bool lambda_compare(const cTrial& a, const cTrial& b)
	{
		if (a.lambda < b.lambda) return true;
		if (a.lambda > b.lambda) return false;
		if (a.stepfactor < b.stepfactor) return true;
		return false;
	}

	static bool phid_compare(const cTrial& a, const cTrial& b)
	{
		if (a.phid < b.phid) return true;
		if (a.phid > b.phid) return false;
		if (a.stepfactor < b.stepfactor) return true;
		return false;
	}
};

class cTrialCache {

private:
	std::vector<cTrial> trials;

public:

	double target = 0;

	const size_t size() const {
		return trials.size();
	}

	const cTrial& trial(const size_t& index) const {
		return trials[index];
	}

	void insert(cTrial& t) {
		t.insert_order = size();
		trials.push_back(t);
	}

	double sfsearch(const double& x) {
		for (size_t k = 0; k < trials.size(); ++k) {
			if (x == trials[k].stepfactor)return trials[k].phid;
		}
		return -1.0;
	}

	size_t minphidindex() {
		size_t ind = 0;
		for (size_t k = 1; k < trials.size(); ++k) {
			if (trials[k].phid < trials[ind].phid) {
				ind = k;
			}
		}
		return ind;
	}

	size_t maxphidindex() {
		size_t ind = 0;
		for (size_t k = 1; k < trials.size(); ++k) {
			if (trials[k].phid > trials[ind].phid) {
				ind = k;
			}
		}
		return ind;
	}

	size_t maxlambdaindex() {
		size_t ind = 0;
		for (size_t k = 1; k < trials.size(); ++k) {
			if (trials[k].lambda > trials[ind].lambda) {
				ind = k;
			}
		}
		return ind;
	}

	cTrial findlambda(const double lambda) {
		for (size_t i = 0; i < trials.size(); i++) {
			if (trials[i].lambda == lambda) {
				return trials[i];
			}
		}
		return trials[0];
	}

	double minphid() { return trials[minphidindex()].phid; }

	double maxphid() { return trials[maxphidindex()].phid; }

	cTrial minphidtrial() { return trials[minphidindex()]; }

	cTrial maxphidtrial() { return trials[maxphidindex()]; }

	cTrial maxlambdatrial() { return trials[maxlambdaindex()]; }

	bool is_target_braketed()
	{
		bool below = false;
		bool above = false;
		for (size_t i = 0; i < trials.size(); i++) {
			if (trials[i].phid <= target) {
				below = true;
			}
			else if (trials[i].phid >= target) {
				above = true;
			}
			if (above && below) return true;
		}
		return false;
	}

	bool is_min_braketed()
	{
		size_t index = minphidindex();
		if (index == 0 || index == trials.size() - 1) {
			return false;
		}

		double fa = trials[index - 1].phid;
		double fb = trials[index].phid;
		double fc = trials[index + 1].phid;
		if ((fb < fa) && (fb < fc)) {
			return true;
		}
		return false;
	}

	void sort_lambda()
	{
		std::sort(trials.begin(), trials.end(), cTrial::lambda_compare);
	}

	void sort_phid()
	{
		std::sort(trials.begin(), trials.end(), cTrial::phid_compare);
	}

	std::string info_string(const double& current_lambda, const double& current_phid)
	{
		sort_lambda();
		std::ostringstream s;
		s << "Trials: Current Lambda = " << current_lambda << " Current Phid = " << current_phid << " Target = " << target << std::endl;;
		s << "   N           Lambda       Stepfactor             Phid\n";
		for (size_t i = 0; i < trials.size(); i++) {
			s << ixd(4) << trials[i].insert_order << " " <<
				std::setw(16) << trials[i].lambda << " " <<
				std::setw(16) << trials[i].stepfactor << " " <<
				std::setw(16) << trials[i].phid << std::endl;
		}
		return s.str();
	}

	void print(const double& current_lambda, const double& current_phid)
	{
		std::cout << info_string(current_lambda, current_phid);
		std::cerr << info_string(current_lambda, current_phid);
	}
};

//State at end of an Iteration
class cIterationState {

public:
	size_t iteration = 0;
	double lambda = 0.0;
	double targetphid = 0.0;
	double phid = 0.0;
	double phim = 0.0;

	Vector pred;
	Vector param;

	std::string info_string() const {
		std::ostringstream ss;
		ss << "Iteration " << iteration << std::endl;
		ss << "Lambda " << lambda << std::endl;
		ss << "PhiD " << phid << std::endl;
		ss << "TargetPhiD " << targetphid << std::endl;
		ss << "PhiM " << phim << std::endl;
		//ss << "PhiC " << phic << std::endl;
		//ss << "PhiT " << phit << std::endl;
		//ss << "PhiG " << phig << std::endl;
		//ss << "PhiVC " << phivc << std::endl;
		//ss << "PhiQC " << phiqc << std::endl;
		//ss << "PhiLC " << philc << std::endl;
		//ss << "PhiLG " << philg << std::endl;
		return ss.str();
	};
};

class cInverter {

private:

	cTrial stepfactor_trial(const double& lambda, const Vector& dm, const double& stepfactor)
	{
		cTrial t;
		Vector p = CIS.param + stepfactor * dm;
		Vector g(nData);
		forwardmodel(p, g);

		t.phid = phiData(g);
		t.lambda = lambda;
		t.stepfactor = stepfactor;
		return t;
	}

	cTrial stepfactor_search(const double& lambda, const Vector& dm, const double& target)
	{
		cInversionLineSearcher s(target);
		s.set_maxtrials(20);
		s.set_ytol(target * 0.01);
		s.add_pair(0.0, CIS.phid);

		cTrial t;

		//bookmark
		double maxsf = 1.0;

		t = stepfactor_trial(lambda, dm, 0.10);
		s.add_pair(t.stepfactor, t.phid);

		if (t.phid > CIS.phid) {
			t = stepfactor_trial(lambda, dm, 0.01);
			s.add_pair(t.stepfactor, t.phid);
			if (t.phid > CIS.phid) {
				t = stepfactor_trial(lambda, dm, 0.001);
				s.add_pair(t.stepfactor, t.phid);
			}
		}

		if (t.phid <= CIS.phid && t.phid > target) {
			t = stepfactor_trial(lambda, dm, 0.25);
			s.add_pair(t.stepfactor, t.phid);
		}

		if (t.phid <= CIS.phid && t.phid > target) {
			t = stepfactor_trial(lambda, dm, 0.5);
			s.add_pair(t.stepfactor, t.phid);
		}

		if (t.phid <= CIS.phid && t.phid > target) {
			t = stepfactor_trial(lambda, dm, maxsf);
			s.add_pair(t.stepfactor, t.phid);
		}

		double sf;
		while (s.next_x(sf)) {
			if (sf >= maxsf || sf <= 0.0) {
				if (Verbose) {
					std::ostringstream msg;
					msg << "\nAuto step factor " << sf << " outside range 0.0 to " << maxsf << std::endl;
					std::cout << msg.str();
					std::cerr << msg.str();
				}
				break;//only search in the range 0 - 1
			}
			t = stepfactor_trial(lambda, dm, sf);
			s.add_pair(t.stepfactor, t.phid);
		}


		DoublePair p = s.nearest();
		t.lambda = lambda;
		t.stepfactor = p.first;
		t.phid = p.second;

		if (Verbose) {
			std::ostringstream msg;
			msg << "Trial Lambda=" << lambda << std::endl;
			msg << s << std::endl;
			msg << "Found ";
			msg << "  Lambda: " << t.lambda;
			msg << "  Step Factor: " << t.stepfactor;
			msg << "  Phid: " << t.phid << std::endl;
			std::cout << msg.str();
			std::cerr << msg.str();
		};
		return t;
	}

	double lambda_trial_function(cTrialCache& T, const double& lambda, bool dostepfactorsearch)
	{
		Vector dm = parameter_change(lambda, CIS.param, CIS.pred);
		cTrial t;
		if (dostepfactorsearch) {
			t = stepfactor_search(lambda, dm, T.target);
		}
		else {
			t = stepfactor_trial(lambda, dm, 1.0);
		}
		T.insert(t);
		return t.phid;
	}

	double goldensearch(double a, double b, double c, double xtol, const double lambda, const Vector& m, const Vector& dm, Vector& g, cTrialCache& cache)
	{
		//adapted from http://en.wikipedia.org/wiki/Golden_section_search	
		const double resphi = 2 - ((1 + sqrt(5.0)) / 2.0);
		double x;
		if (c - b > b - a) {
			x = b + resphi * (c - b);
		}
		else {
			x = b - resphi * (b - a);
		}

		//if(fabs(c - a) < tau * (fabs(b) + fabs(x))){
		//  return (c + a) / 2; 
		//}
		if (fabs(c - a) < xtol) {
			return (c + a) / 2;
		}

		double fx = cache.sfsearch(x);
		if (fx < 0) {
			cTrial t;
			Vector p = m + x * dm;
			forwardmodel(p, g);
			fx = phiData(g);
			t.stepfactor = x;
			t.phid = fx;
			t.lambda = lambda;
			cache.insert(t);
		}

		double fb = cache.sfsearch(b);
		if (fb < 0) {
			cTrial t;
			Vector p = m + b * dm;
			forwardmodel(p, g);
			fb = phiData(g);
			t.stepfactor = b;
			t.phid = fb;
			t.lambda = lambda;
			cache.insert(t);
		}

		assert(fx != fb);
		if (fx < fb) {
			if (c - b > b - a) {
				return goldensearch(b, x, c, xtol, lambda, m, dm, g, cache);
			}
			else {
				return goldensearch(a, x, b, xtol, lambda, m, dm, g, cache);
			}
		}
		else {
			if (c - b > b - a) {
				return goldensearch(a, b, x, xtol, lambda, m, dm, g, cache);
			}
			else {
				return goldensearch(x, b, c, xtol, lambda, m, dm, g, cache);
			}
		}
	}

	eBracketResult brackettarget(cTrialCache& T, const double target, const double currentlambda)
	{
		bool dostepfactorsearch = true;
		cInversionLineSearcher s(target);
		s.set_ytol(target * 0.01);
		s.set_maxtrials(10);
		std::vector<double> xlist;
		if (CIS.iteration == 0) {
			xlist = increment<double>(32, 8.0, -0.25);
			s.set_maxtrials(xlist.size() * 2);
		}
		else {
			s.set_maxtrials(10);
			xlist = {
				std::log10(currentlambda * 2.0),
				std::log10(currentlambda),
				std::log10(currentlambda * 0.5),
				std::log10(currentlambda * 0.25)
			};
		}

		double x;
		for (size_t i = 0; i < xlist.size(); i++) {
			double x = xlist[i];
			double phid = lambda_trial_function(T, pow10(x), dostepfactorsearch);
			if (Verbose) T.print(CIS.lambda, CIS.phid);
			s.add_pair(x, phid);
			if (T.is_target_braketed()) {
				return eBracketResult::BRACKETED;//target bracketed		
			}
		}

		while (s.next_x(x)) {
			double phid = lambda_trial_function(T, pow10(x), dostepfactorsearch);
			if (Verbose) T.print(CIS.lambda, CIS.phid);
			s.add_pair(x, phid);
			if (T.is_target_braketed()) {
				return eBracketResult::BRACKETED;//target bracketed		
			}
		}

		if (T.maxphid() < target) {
			return eBracketResult::ALLBELOW;//all below target	
		}
		else if (T.is_min_braketed()) {
			return eBracketResult::MINBRACKETED;//min bracketed											
		}
		else return eBracketResult::ALLABOVE;//all above target
	}

	double brentsmethod(cTrialCache& T, const double target, double& newphid)
	{
		//Adapted from http://en.wikipedia.org/wiki/Brent's_method
		double xerrorTol = 0.01;//in log10 
		double yerrorTol = target * 0.01;//1% accuracy is good enough

		T.sort_lambda();
		size_t index = T.size() - 1;
		for (size_t i = T.size() - 1; i >= 1; i--) {
			double f1 = T.trial(i).phid - target;
			double f2 = T.trial(i - 1).phid - target;
			if (f1 * f2 <= 0.0) {
				index = i;
				break;
			}
		}

		double a = std::log10(T.trial(index - 1).lambda);
		double b = std::log10(T.trial(index).lambda);
		double fa = T.trial(index - 1).phid - target;
		double fb = T.trial(index).phid - target;
		if (fa * fb >= 0.0) {
			std::string msg = strprint("Warning: Target must be bracketed\n");
			glog.logmsg(msg);
			std::cerr << msg << std::endl;
		}

		double c = 0;
		double d = std::numeric_limits<double>::max();
		double fc = 0;
		double s = 0;
		double fs = 0;

		// if f(a) f(b) >= 0 then error-exit
		if (fa * fb >= 0)
		{
			if (fa < fb) {
				newphid = fa + target;
				return pow10(a);
			}
			else {
				newphid = fb + target;
				return pow10(b);
			}
		}

		// if |f(a)| < |f(b)| then swap (a,b) end if
		if (fabs(fa) < fabs(fb)) {
			double tmp;
			tmp = a;   a = b;  b = tmp;
			tmp = fa; fa = fb; fb = tmp;
		}

		c = a;
		fc = fa;
		bool mflag = true;
		int i = 0;

		while (!(fb == 0) && (fabs(a - b) > xerrorTol))
		{
			if ((fa != fc) && (fb != fc))
				// Inverse quadratic interpolation
				s = a * fb * fc / (fa - fb) / (fa - fc) + b * fa * fc / (fb - fa) / (fb - fc) + c * fa * fb / (fc - fa) / (fc - fb);
			else
				// Secant Rule
				s = b - fb * (b - a) / (fb - fa);

			double tmp2 = (3 * a + b) / 4;
			if ((!(((s > tmp2) && (s < b)) || ((s < tmp2) && (s > b)))) || (mflag && (fabs(s - b) >= (fabs(b - c) / 2))) || (!mflag && (fabs(s - b) >= (fabs(c - d) / 2))))
			{
				s = (a + b) / 2;
				mflag = true;
			}
			else
			{
				if ((mflag && (fabs(b - c) < xerrorTol)) || (!mflag && (fabs(c - d) < xerrorTol)))
				{
					s = (a + b) / 2;
					mflag = true;
				}
				else
					mflag = false;
			}
			fs = lambda_trial_function(T, pow10(s), true) - target;
			if (fabs(fs) < yerrorTol) {
				newphid = fs + target;
				return pow10(s);
			}

			d = c;
			c = b;
			fc = fb;
			if (fa * fs < 0) { b = s; fb = fs; }
			else { a = s; fa = fs; }

			// if |f(a)| < |f(b)| then swap (a,b) end if
			if (fabs(fa) < fabs(fb))
			{
				double tmp = a; a = b; b = tmp; tmp = fa; fa = fb; fb = tmp;
			}
			i++;
			if (i > 20) {
				newphid = fb + target;

				std::string msg = strprint("Warning: Too many bisections\n");
				glog.logmsg(msg);
				std::cerr << msg << std::endl;
				return pow10(b);
			}
		}
		if (fb < fa) {
			newphid = fb + target;
			return pow10(b);
		}
		else {
			newphid = fa + target;
			return pow10(a);
		}
	}

protected:
	size_t nForwards = 0;
	size_t nJacobians = 0;

	std::string CommandLine;
	int Size;
	int Rank;
	bool UsingOpenMP;
	bool Verbose = false;
	std::string OutputMessage;
	cBlock Control;

	size_t nData;
	size_t nParam;

	cIterationState CIS;
	Vector Obs;
	Vector Err;
	Vector RefParam;
	Vector RefParamStd;
	Vector ParameterSensitivity;
	Vector ParameterUncertainty;
	Vector Param_Min;
	Vector Param_Max;

	Matrix J;
	Matrix Wd;
	Matrix Wm;

	eNormType  NormType;

	double MinimumPhiD;//overall	
	double MinimumImprovement;//			
	size_t MaxIterations;
	std::string TerminationReason;
	cSampleBunch Bunch;
	std::unique_ptr<cInputManager> IM;
	std::unique_ptr<cOutputManager> OM;
	std::vector<size_t> ActiveData;//indices of the data that are not culled	

public:

	cInverter(const std::string& controlfile, const int& size, const int& rank, const bool& usingopenmp, const std::string commandline)
	{
		Size = size;
		Rank = rank;
		UsingOpenMP = usingopenmp;
		CommandLine = commandline;
	};

	virtual ~cInverter()
	{
		glog.close();
		//std::cout << "Destroying cInverter\n";
	};

	void initialise(const std::string& controlfile)
	{
		_GSTITEM_
			loadcontrolfile(controlfile);
		set_field_definitions();
		setup_data();
		setup_parameters();
		execute();
	}

	virtual void loadcontrolfile(const std::string& filename) = 0;
	virtual void set_field_definitions() = 0;
	virtual void setup_data() = 0;
	virtual void setup_parameters() = 0;
	virtual int  execute() = 0;
	virtual std::string dumppath() const = 0;

	virtual void forwardmodel(const Vector& parameters, Vector& predicted) = 0;
	virtual void forwardmodel_and_jacobian(const Vector& parameters, Vector& predicted, Matrix& jacobian) = 0;
	virtual Vector parameter_change(const double& lambda, const Vector& m_old, const Vector& g_old) = 0;

	inline static bool isnull(const double& val) {
		//if (val == 29322580645161.28516) return true;
		//else if (val == 293196483870967808.00000) return true;		
		if (val > 1e14) return true;
		else if (val < -1e14) return true;
		else if (isdefined(val)) return false;
		else return true;
	}

	void set_fftw_lock() {
#if defined _OPENMP
		//If OpenMP is being used set the thread lock while FFTW initialises
		if (UsingOpenMP) {
			omp_set_lock(&fftw_thread_lock);
		}
#endif	
	}

	void unset_fftw_lock() {
#if defined _OPENMP			
		if (UsingOpenMP) {
			omp_unset_lock(&fftw_thread_lock);
		}
#endif	
	}

	double l1_norm(const Vector& g)
	{
		double l1 = 0.0;
		for (size_t i = 0; i < nData; i++) {
			l1 += std::abs(Obs[i] - g[i]) / Err[i];
		}
		return l1 / nData;
	}

	double l2_norm(const Vector& g)
	{
		Vector v = Obs - g;
		Vector a = Wd * v;
		double l2 = mtDm(v, Wd);
		return l2;
	}

	double phiData(const Vector& g)
	{
		double phid;
		if (NormType == eNormType::L1) {
			phid = l1_norm(g);
		}
		else {
			phid = l2_norm(g);
		}
		//This reports invalid models 
		if (phid < 0.0) {
			phid = 1e9;
			std::string msg = strprint("Caught invalid PhiD\n");
			glog.logmsg(msg);
			std::cerr << msg << std::endl;
		}
		return phid;
	}

	virtual double phiModel(const Vector& m) = 0;

	virtual Vector solve_linear_system(const double& lambda, const Vector& param, const Vector& pred) = 0;

	void write(const Vector& v, std::string path) const
	{
		FILE* fp = fileopen(path, "w");
		for (auto i = 0; i < v.rows(); i++) {
			fprintf(fp, "%le\n", v[i]);
		}
		fclose(fp);
	}

	void write(const std::vector<double>& v, std::string path) const
	{
		FILE* fp = fileopen(path, "w");
		for (size_t i = 0; i < v.size(); i++) {
			fprintf(fp, "%le\n", v[i]);
		}
		fclose(fp);
	}

	std::string rec_it_str() const
	{
		std::ostringstream os;
		os << "Record " << Bunch.master_record() << " It " << CIS.iteration;
		return os.str();
	};

	cTrial lambda_search_targetphid(const double& currentlambda, const double& targetphid)
	{
		cTrialCache T;
		T.target = targetphid;
		eBracketResult b = brackettarget(T, targetphid, currentlambda);
		cTrial t;
		if (b == eBracketResult::BRACKETED) {
			//bracketed target - find with Brents Method
			double newphid = undefinedvalue<double>();
			double lambda = brentsmethod(T, targetphid, newphid);
			t = T.findlambda(lambda);
		}
		else if (b == eBracketResult::MINBRACKETED) {
			//bracketed minimum but above target - take smallest phid
			t = T.minphidtrial();
		}
		else if (b == eBracketResult::ALLBELOW) {
			//all below target	- take the largest lambda			
			t = T.maxlambdatrial();
		}
		else if (b == eBracketResult::ALLABOVE) {
			//all above target - take smallest phid
			t = T.minphidtrial();
		}
		else {
			glog.errormsg(_SRC_, "lambda_search_target(): Unknown value %d returned from target brackettarget()\n", b);
		}
		if (Verbose) T.print(CIS.lambda, CIS.phid);
		return t;
	}

	Vector compute_parameter_sensitivity()
	{
		Vector s = Vector::Zero(nParam);
		for (size_t pi = 0; pi < nParam; pi++) {
			for (size_t di = 0; di < nData; di++) {
				s[pi] += (std::fabs(J(di, pi)) * std::sqrt((double)nData * Wd(di, di)));
			}
		}
		return s;
	}

	Vector compute_parameter_uncertainty()
	{
		Matrix JtWdJ = J.transpose() * Wd * J;
		Matrix iCm = Matrix::Zero(nParam, nParam);
		for (size_t i = 0; i < nParam; i++) {
			iCm(i, i) = 1.0 / (RefParamStd[i] * RefParamStd[i]);
		}
		Matrix X = (double)nData * JtWdJ + iCm;
		Matrix pinvX = pseudoInverse(X);
		Vector s(nParam);
		for (size_t i = 0; i < nParam; i++) {
			s[i] = std::sqrt(pinvX(i, i));
		}
		return s;
	}

	void initialise_Wd() {
		Wd = Matrix::Zero(nData, nData);
		double s = 1.0 / (double)nData;
		for (size_t i = 0; i < nData; i++) {
			Wd(i, i) = s / (Err[i] * Err[i]);
		}
	}
};

#endif
