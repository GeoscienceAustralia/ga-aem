/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cmath>
#include <ctime>
#include <cfloat>
#include <algorithm>

#if defined(_MPI_ENABLED)
	#include "mpi.h"
#endif
#include "general_utils.h"
#include "random_utils.h"
//#include "eigen_utils.h"
#include "rjmcmc1d.h"

rjMcMC1DSampler::rjMcMC1DSampler()
{

}

void rjMcMC1DSampler::reset()
{
	description = "";
	starttime = "";
	endtime = "";
	samplingtime = 0;
	np_valuechange = 0; np_birth = 0; np_death = 0; np_move = 0; np_nuisancechange = 0;
	na_valuechange = 0; na_birth = 0; na_death = 0; na_move = 0; na_nuisancechange = 0;
	pmap.resettozero();
	nmap.resettozero();	
}
bool rjMcMC1DSampler::isreportable(size_t si)
{
	if (si == 0 || si == nsamples - 1)return true;
	size_t k = (size_t)pow(10, floor(log10((double)si)));
	if (k > 1000)k = 1000;
	if (si%k == 0)return true;
	return false;
}
void rjMcMC1DSampler::printstats(size_t ci, size_t si, size_t np, double nmf)
{
	printf("ci=%4zu si=%7zu np=%2zu nmf=%10g vc=%5.2lf m=%5.2lf b=%5.2lf d=%5.2lf n=%5.2lf\n", ci, si, np, nmf, ar_valuechange(), ar_move(), ar_birth(), ar_death(), ar_nuisancechange());
}
size_t rjMcMC1DSampler::nnuisances()
{
	return nuisance_init.size();
}
bool rjMcMC1DSampler::isinbounds(double& bmin, double& bmax, double& b)
{
	if (b<bmin || b>bmax)return false;
	else return true;
}
bool rjMcMC1DSampler::includeinmap(size_t si)
{
	if (si < nburnin)return false;
	size_t k = si - nburnin;
	if ((k%thinrate) == 0){
		return true;
	}
	return false;
}
double rjMcMC1DSampler::l2misfit(const std::vector<double>& g)
{
	double sum = 0.0;
	for (size_t di = 0; di < ndata; di++){
		double nr = (obs[di] - g[di]) / err[di];
		sum += nr*nr;
	}
	return sum;
}
double rjMcMC1DSampler::computemisfit(const rjMcMC1DModel& m)
{
	std::vector<double> pred = forwardmodel(m);
	return l2misfit(pred);
}

void rjMcMC1DSampler::set_misfit(rjMcMC1DModel& m)
{
	double mfit = computemisfit(m);
	//double mfit = (double)ndata;
	m.setmisfit(mfit);
}
void rjMcMC1DSampler::writeresults(FILE* fp)
{
	writeheader(fp);
	writemodel(fp, mHighestLikelihood);
	writemodel(fp, mLowestMisfit);
	pmap.writedata(fp);
	nmap.writedata(fp);
	for (size_t ci = 0; ci < nchains; ci++){
		mChainInfo[ci].writeconvergencedata(ci,fp);
		mChainInfo[ci].writemodelchain_binary(ci,fp);
	}
}
void rjMcMC1DSampler::writeheader(FILE* fp)
{

	std::string pparm = "LINEAR";
	if (param_position == LOG10) pparm = "LOG10";

	std::string vparm = "LINEAR";
	if (param_value == LOG10) vparm = "LOG10";

	fprintf(fp, "%s\n", description.c_str());
	fprintf(fp, "Start time %s\n", starttime.c_str());
	fprintf(fp, "End time %s\n", endtime.c_str());
	fprintf(fp, "%lf sec sampling time \n", samplingtime);

	fprintf(fp, "%zu nchains\n", nchains);
	fprintf(fp, "%zu nsamples\n", nsamples);
	fprintf(fp, "%zu nburnin\n", nburnin);
	fprintf(fp, "%zu nnuisance\n", nnuisances());

	fprintf(fp, "%zu nl_min\n", nl_min);
	fprintf(fp, "%zu nl_max\n", nl_max);

	fprintf(fp, "%s position parameterization\n", pparm.c_str());
	fprintf(fp, "%lf pmin\n", 0.0);
	fprintf(fp, "%lf pmax\n", pmax);
	fprintf(fp, "%zu npcells\n", pmap.npbins());

	fprintf(fp, "%s value parameterization\n", vparm.c_str());
	fprintf(fp, "%lf vmin\n", vmin);
	fprintf(fp, "%lf vmax\n", vmax);
	fprintf(fp, "%zu nvcells\n", pmap.nvbins());

	fprintf(fp, "%lf sd_value\n", sd_valuechange);
	fprintf(fp, "%lf sd_move\n", sd_move);
	fprintf(fp, "%lf sd_bd_value\n", sd_bd_valuechange);

	fprintf(fp, "%lf Acceptance_Rate_birth\n", ar_birth());
	fprintf(fp, "%lf Acceptance_Rate_death\n", ar_death());
	fprintf(fp, "%lf Acceptance_Rate_value\n", ar_valuechange());
	fprintf(fp, "%lf Acceptance_Rate_move\n", ar_move());
	fprintf(fp, "%lf Acceptance_Rate_nuisance\n", ar_nuisancechange());

}
void rjMcMC1DSampler::writemodel(FILE* fp, const rjMcMC1DModel& m)
{
	fprintf(fp, "%.6e %.6e %zu", m.misfit(), m.logppd(), m.nlayers());
	for (size_t li = 0; li < m.nlayers(); li++){
		fprintf(fp, " %.6e", m.layers[li].ptop);
	}
	for (size_t li = 0; li < m.nlayers(); li++){
		fprintf(fp, " %.6e", m.layers[li].value);
	}
	fprintf(fp, "\n");
}
bool rjMcMC1DSampler::saveconvergencerecord(size_t si)
{
	if (si == 0 || si == nsamples - 1)return true;
	size_t k = (size_t)pow(10, floor(log10((double)si)));
	if (k > thinrate) k = thinrate;
	if (si%k == 0)return true;
	return false;
}

bool rjMcMC1DSampler::propose_valuechange(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
{
	np_valuechange++;

	size_t index = irand((size_t)0,mcur.nlayers() - 1);

	double logstd = DEFAULTLOGSTDDECADES;
	double vold = mcur.layers[index].value;
	double vnew;
	double pqratio;
	if (param_value == LINEAR){
		double m = (std::pow(10.0, logstd) - std::pow(10.0, -logstd)) / 2.0;
		vnew = vold + m * vold * nrand<double>();
		double qpdfforward = gaussian_pdf(vold, m*vold, vnew);
		double qpdfreverse = gaussian_pdf(vnew, m*vnew, vold);
		pqratio = qpdfreverse / qpdfforward;
	}
	else{
		vnew = vold + logstd*nrand<double>();
		pqratio = 1.0;
	}

	bool isvalid = isinbounds(vmin, vmax, vnew);
	if (isvalid == false)return false;
	mpro.layers[index].value = vnew;

	set_misfit(mpro);
	double logpqratio = log(pqratio);
	double logliker = -(mpro.misfit() - mcur.misfit()) / 2.0;
	double logar = logpqratio + logliker;
	if (std::log(urand<double>()) < logar){
		na_valuechange++;
		return true;
	}
	return false;
}
bool rjMcMC1DSampler::propose_move(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
{
	np_move++;

	size_t n = mcur.nlayers();
	if (n <= 1)return false;

	size_t index = irand((size_t)1, n - 1);
	double pold = mcur.layers[index].ptop;

	//double std = sd_move;	
	double std = DEFAULTMOVESTDFRACTION*pold;
	double pnew = pold + std*nrand<double>();
	double qpdfforward = gaussian_pdf(pold, pold*DEFAULTMOVESTDFRACTION, pnew);
	double qpdfreverse = gaussian_pdf(pnew, pnew*DEFAULTMOVESTDFRACTION, pold);


	bool isvalid = mpro.move_interface(index, pnew);
	if (isvalid == false)return false;

	set_misfit(mpro);
	double pqratio = qpdfreverse / qpdfforward;
	double logpqratio = std::log(pqratio);

	double loglr = -(mpro.misfit() - mcur.misfit()) / 2.0;
	double logar = logpqratio + loglr;
	double logu = std::log(urand<double>());
	if (logu < logar){
		na_move++;
		return true;
	}
	return false;
}
bool rjMcMC1DSampler::propose_birth(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
{
	np_birth++;

	size_t n = mcur.nlayers();
	if (n >= nl_max)return false;

	double  pos = urand(0.0, pmax);
	size_t  index = mcur.which_layer(pos);
	double vold = mcur.layers[index].value;
	double vnew, pqratio;

	if (mBirthDeathFromPrior){
		vnew = urand(vmin, vmax);
		pqratio = 1.0;
	}
	else{
		double vcpdf;
		double logstd = DEFAULTLOGSTDDECADES;
		if (param_value == LINEAR){
			double m = (pow(10.0, logstd) - pow(10.0, -logstd)) / 2.0;
			vnew = vold + m*vold*nrand<double>();
			vcpdf = gaussian_pdf(vold, m*vold, vnew);
		}
		else{
			vnew  = vold + logstd*nrand<double>();
			vcpdf = gaussian_pdf(vold, logstd, vnew);
		}
		pqratio = 1.0 / ((vmax - vmin)*vcpdf);
	}
	//pqratio *= (double)(n) / double(n + 1);

	bool   isvalid = mpro.insert_interface(pos, vnew);
	if (isvalid == false)return false;
	set_misfit(mpro);

	double logpqratio = std::log(pqratio);
	double loglikeratio = -(mpro.misfit() - mcur.misfit()) / 2.0;
	double logar = logpqratio + loglikeratio;
	if (std::log(urand<double>()) < logar){
		na_birth++;
		return true;
	}
	return false;
}
bool rjMcMC1DSampler::propose_death(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
{
	np_death++;

	size_t n = mcur.nlayers();
	if (n <= nl_min)return false;

	size_t index = irand((size_t)1, n - 1);
	bool isvalid = mpro.delete_interface(index);
	if (isvalid == false)return false;
	set_misfit(mpro);

	double pqratio;
	if (mBirthDeathFromPrior){
		pqratio = 1.0;
	}
	else{
		double logstd = DEFAULTLOGSTDDECADES;
		double vnew = mcur.layers[index - 1].value;
		double vold = mcur.layers[index].value;
		double vcpdf;
		if (param_value == LINEAR){
			double m = (pow(10.0, logstd) - pow(10.0, -logstd)) / 2.0;
			vcpdf = gaussian_pdf(vnew, m*vnew, vold);
		}
		else{
			vcpdf = gaussian_pdf(vnew, logstd, vold);
		}
		pqratio = (vmax - vmin)*vcpdf;
	}
	//pqratio *= (double)(n) / double(n - 1);

	double logpqratio = log(pqratio);
	double loglikeratio = -(mpro.misfit() - mcur.misfit()) / 2.0;
	double logar = logpqratio + loglikeratio;
	if (log(urand<double>()) < logar){
		na_death++;
		return true;
	}
	return false;
}
bool rjMcMC1DSampler::propose_nuisancechange(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
{
	np_nuisancechange++;

	size_t ni = irand((size_t)0,mcur.nnuisances() - 1);
	double delta = nrand<double>() * mcur.nuisances[ni].sd_valuechange;
	double nv = mcur.nuisances[ni].value + delta;;
	bool isvalid = isinbounds(mcur.nuisances[ni].min, mcur.nuisances[ni].max, nv);
	if (isvalid == false)return false;

	mpro.nuisances[ni].value = nv;

	set_misfit(mpro);
	double logar = -(mpro.misfit() - mcur.misfit()) / 2.0;
	double logu = log(urand<double>());
	if (logu < logar){
		na_nuisancechange++;
		return true;
	}
	return false;
}
bool rjMcMC1DSampler::propose_independent(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
{
	mpro = choosefromprior();
	set_misfit(mpro);

	size_t npro = mpro.nlayers() - 1;
	size_t ncur = mcur.nlayers() - 1;
	double a = pow(1.0 / (vmax - vmin), npro - ncur);
	double b = pow(1.0 / pmax, npro - ncur);
	double c = (double)factorial((unsigned int)npro) / (double)factorial((unsigned int)ncur);
	double priorratio = a*b*c;
	priorratio = 1.0;

	double logpriorratio = std::log(priorratio);
	double logliker = -(mpro.misfit() - mcur.misfit()) / 2.0;
	double logar = logpriorratio + logliker;
	double logu = std::log(urand<double>());
	if (logu < logar){		
		return true;
	}
	return false;
}

rjMcMC1DModel rjMcMC1DSampler::choosefromprior()
{
	rjMcMC1DModel m;
	size_t nl = irand(nl_min, nl_max);
	m.initialise(pmax, vmin, vmax);
	for (size_t li = 0; li < nl; li++){
		bool status = false;
		while (status == false){
			double pos = urand(0.0, pmax);//position of interface		
			double value = urand(vmin, vmax);//value
			status = m.insert_interface(pos, value);
		}
	}
	m.nuisances = nuisance_init;
	return m;
}

void rjMcMC1DSampler::sample()
{	
	mBirthDeathFromPrior = false;
	
	starttime = timestamp();
	double t1 = gettime();
	mChainInfo.resize(nchains);
	for (size_t ci = 0; ci < nchains; ci++){				
		mChainInfo[ci].reset();
		//mChainInfo[ci].modelchain.resize(nsamples);
		rjMcMC1DModel mcur;		
		for (size_t si = 0; si < nsamples; si++){

			//Initialis chain
			if (si == 0){
				mcur = choosefromprior();
				set_misfit(mcur);
			}

			//Initialise "best" models
			if (ci == 0 && si==0){
				mHighestLikelihood = mcur;
				mLowestMisfit = mcur;
			}
			
			rjMcMC1DModel mpro = mcur;
			size_t nopt = 4;
			if (mcur.nnuisances() > 0)nopt = 5;
			size_t option = irand((size_t)0, nopt - 1);
			
			bool accept = false;
			if (option == 0){
				accept = propose_valuechange(mcur, mpro);
			}
			else if (option == 1){
				accept = propose_move(mcur, mpro);
			}
			else if (option == 2){
				accept = propose_birth(mcur, mpro);
			}
			else if (option == 3){
				accept = propose_death(mcur, mpro);
			}
			else if (option == 4){
				accept=propose_nuisancechange(mcur,mpro);				
			}
			else if (option == 5){
				accept = propose_independent(mcur, mpro);
			}
			else{
				exit(1);
				break;
			}

			if (accept) mcur = mpro;				
						
			if (includeinmap(si)){
				pmap.addmodel(mcur);
				nmap.addmodel(mcur);
				mChainInfo[ci].modelchain.push_back(mcur);
			}

			if (mcur.logppd() > mHighestLikelihood.logppd()){
				mHighestLikelihood = mcur;
			}

			if (mcur.misfit() < mLowestMisfit.misfit()){
				mLowestMisfit = mcur;
			}

			if (saveconvergencerecord(si)){
				mChainInfo[ci].sample.push_back((uint32_t)si);
				mChainInfo[ci].nlayers.push_back((uint32_t)mcur.nlayers());
				mChainInfo[ci].misfit.push_back((float)mcur.misfit());
				mChainInfo[ci].logppd.push_back((float)mcur.logppd());
				mChainInfo[ci].ar_valuechange.push_back((float)ar_valuechange());
				mChainInfo[ci].ar_move.push_back((float)ar_move());
				mChainInfo[ci].ar_birth.push_back((float)ar_birth());
				mChainInfo[ci].ar_death.push_back((float)ar_death());
				mChainInfo[ci].ar_nuisancechange.push_back((float)ar_nuisancechange());
			}

			#ifdef _WIN32
			if (reportstats){
				if (mRank == 0 && isreportable(si)){
					//printstats(ci, si, mcur.nlayers(), mcur.misfit() / (double)ndata);
				}
			}
			#endif
		}
	}
	double t2 = gettime();
	endtime = timestamp();
	samplingtime = t2 - t1;
}