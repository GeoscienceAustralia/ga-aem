/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _rjMcMC1D_H_
#define _rjMcMC1D_H_

#include <climits>
#include <cstdint>
#include "general_utils.h"
#include "random_utils.h"

#if defined(_MPI_ENABLED)
	#include "mpi.h"
#endif

enum proposal_type {VALUECHANGE, BIRTH, DEATH, MOVE, NUISANCE};
enum parameterization_type {NOTSET, LINEAR, LOG10};
enum nuisance_type {TX_HEIGHT,TX_ROLL,TX_PITCH,TX_YAW,TXRX_DX,TXRX_DY,TXRX_DZ,RX_ROLL,RX_PITCH,RX_YAW,TXRX_ANGLE,TXRX_DISTANCE,NTYPE_UNKNOWN};

class rjMcMC1DLayer{

public:
	double ptop;
	double value;
	int operator==(rjMcMC1DLayer a){
		if(ptop==a.ptop)return 1;
		else return 0;
	};
	int operator<(const rjMcMC1DLayer a) const {
		if(ptop<a.ptop)return 1;
		else return 0;
	};
	int operator>(const rjMcMC1DLayer a) const {
		if(ptop>a.ptop)return 1;
		else return 0;
	};
};

class rjMcMCNuisance{

public:	
	nuisance_type type;
	double value;
	double min;
	double max;
	double sd_valuechange;	
	
	static std::string ntype2str(nuisance_type t){
		if(t==TX_HEIGHT)return "tx_height";
		else if(t==TX_ROLL)return "tx_roll";
		else if(t==TX_PITCH)return "tx_pitch";
		else if(t==TX_YAW)return "tx_yaw";
		else if(t==TXRX_DX)return "txrx_dx";
		else if(t==TXRX_DY)return "txrx_dy";
		else if(t==TXRX_DZ)return "txrx_dz";
		else if(t==RX_ROLL)return "rx_roll";
		else if(t==RX_PITCH)return "rx_pitch";
		else if(t==RX_YAW)return "rx_yaw";
		else if(t==TXRX_DISTANCE)return "txrx_distance";
		else if(t==TXRX_ANGLE)return "txrx_angle";		
		else if(t==NTYPE_UNKNOWN)return "UnknownNuisanceType";
		else{
			return "UnknownNuisanceType";
		}			
	};


	static nuisance_type str2ntype(std::string s){
		for(size_t i=0;i<NTYPE_UNKNOWN;i++){
			nuisance_type nt = (nuisance_type)i;
			if(strcasecmp(s,ntype2str(nt))==0){
				return nt;
			}			
		}		
		return NTYPE_UNKNOWN;		
	}

	void settype(std::string s){
		type = str2ntype(s);		
	}

	std::string typestring() const{
		return ntype2str(type);
	}
};

class rjMcMC1DModel{
				
	double pmax;
	double vmin;
	double vmax;	
	double mMisfit; 	

public:		
	std::vector<rjMcMC1DLayer>  layers;
	std::vector<rjMcMCNuisance> nuisances;

	void initialise(double maxp, double minv, double maxv)
	{
		layers.clear();
		nuisances.clear();

		mMisfit = DBL_MAX;
		pmax = maxp;		

		vmin = minv;
		vmax = maxv;
	}
	size_t nlayers() const {return layers.size();}
	size_t nnuisances() const {return nuisances.size();}
	size_t nparams() const
	{
		size_t np = 2 * nlayers() + nnuisances();
		return np;
	}
	double misfit() const { return mMisfit; }
	void setmisfit(double mfit) { mMisfit=mfit; }
	double logppd() const { return -mMisfit / 2.0 - log((double)nparams()); }
	void sort_layers()
	{
		std::sort(layers.begin(), layers.end());
	}
	size_t which_layer(double pos) const
	{
		for (size_t li = 0; li<nlayers() - 1; li++){
			if (pos < layers[li + 1].ptop){
				return li;
			}
		}
		return nlayers() - 1;
	}
	bool move_interface(size_t index, double pnew)
	{
		mMisfit = DBL_MAX;
		if(index == 0)return false;
		if(index >= nlayers())return false;
		
		if(pnew<=0)return false;
		if(pnew>=pmax)return false;
				
		layers[index].ptop = pnew;
		sort_layers();
		return true;
	}
	bool insert_interface(double pos, double vbelow)
	{
		mMisfit = DBL_MAX;
		if (pos < 0.0 || pos > pmax)return false;
		if (vbelow < vmin || vbelow > vmax)return false;

		for (size_t li=0; li<nlayers(); li++){
			if (fabs(pos - layers[li].ptop) < DBL_EPSILON){
 				return false;
			}
		}

		rjMcMC1DLayer l;
		if (nlayers() == 0){
			l.ptop  = 0.0;
			l.value = vbelow;
			layers.push_back(l);
		}
		else{
			l.ptop  = pos;
			l.value = vbelow;
			layers.push_back(l);
			sort_layers();
		}
		return true;
	}	
	bool delete_interface(size_t index)
	{
		mMisfit = DBL_MAX;
		if (index <= 0)return false;
		if (index >= nlayers())return false;
		layers[index].value = DBL_MAX;
		layers[index].ptop  = DBL_MAX;
		sort_layers();
		layers.pop_back();
		return true;
	}
	double value(const size_t index){		
			return layers[index].value;
	};
	double thickness(const size_t index){
		if (index < nlayers() - 1){
			double p1 = layers[index].ptop;
			double p2 = layers[index + 1].ptop;			
			return p2 - p1;			
		}
		else if (index == nlayers() - 1){
			return DBL_MAX;
		}
		else{
			glog.errormsg(_SRC_,"Invalid thicknes index");
			exit(1);
		}
	};
	std::vector<double> getvalues() const
	{
		std::vector<double> v;
		v.resize(nlayers());		
		for (size_t i = 0; i<nlayers(); i++)v[i] = layers[i].value;		
		return v;
	}
	std::vector<double> getthicknesses() const
	{
		std::vector<double> t;
		t.resize(nlayers() - 1);		
		for (size_t i = 0; i<nlayers() - 1; i++){
			t[i] = layers[i + 1].ptop - layers[i].ptop;
		}		
		return t;
	}
	void printmodel() const
	{
		printf("nl=%zu\tnn=%zu\tlppd=%lf\tmisfit=%lf\n", nlayers(), nnuisances(), logppd(), misfit());
		for (size_t li = 0; li<nlayers(); li++){
			printf("%4zu\t%10lf\t%10lf\n", li, layers[li].ptop, layers[li].value);
		}
		printf("\n");
	}
	void printmodelex() const
	{
		std::vector<double> c = getvalues();
		std::vector<double> t = getthicknesses();
		printf("nl=%zu\tnn=%zu\tlppd=%lf\tmisfit=%lf\n", nlayers(), nnuisances(), logppd(), misfit());
		double top = 0;
		double bot = 0;
		for (size_t li = 0; li<nlayers(); li++){
			if (li<nlayers() - 1){
				bot = top + t[li];
				printf("%4zu top=%7.2lfm bot=%7.2lfm thk=%7.2lfm conductivity=%6.4lf\n", li, top, bot, t[li], c[li]);
			}
			else{
				printf("%4zu top=%7.2lfm bot= Inf thk= Inf conductivity=%6.4lf\n", li, top, c[li]);
			}
			top = bot;
		}
		for (size_t ni = 0; ni<nnuisances(); ni++){
			printf(" %s = %.3lf\n", nuisances[ni].typestring().c_str(), nuisances[ni].value);
		}
		printf("\n");
	}
};

class rjMcMC1DNuisanceMap{
	
private:
	size_t ns=0;
	std::vector<std::string> typestring;
	
public:
	size_t nsamples() const { return ns; }
	std::vector<std::vector<double>> nuisance;
	void resettozero(){
		for (size_t i = 0; i < nuisance.size(); i++){
			nuisance[i].clear();
		}
	}
	void addmodel(const rjMcMC1DModel& m)
	{
		if (nuisance.size() != m.nnuisances()){
			nuisance.resize(m.nnuisances());
			typestring.resize(m.nnuisances());
			for (size_t i = 0; i < nuisance.size(); i++){
				typestring[i] = m.nuisances[i].typestring();
			}
		}

		for (size_t i = 0; i<nuisance.size(); i++){
			nuisance[i].push_back(m.nuisances[i].value);
		};

		ns++;
	}
	void writedata(FILE* fp)
	{
		size_t nn = nuisance.size();
		for (size_t i = 0; i<nn; i++){
			cStats<double> s(nuisance[i]);
			cHistogram<double, size_t> hist(nuisance[i],s.min,s.max,17);

			fprintf(fp, "%s", typestring[i].c_str());
			fprintf(fp, " %lf", s.min);
			fprintf(fp, " %lf", s.max);
			fprintf(fp, " %lf", s.mean);
			fprintf(fp, " %lf", s.std);
			fprintf(fp, " %zu", nsamples());

			fprintf(fp, " %zu", hist.nbins);
			for (size_t j = 0; j<hist.nbins; j++){
				fprintf(fp, " %lf", hist.centre[j]);
			}
			for (size_t j = 0; j<hist.nbins; j++){
				fprintf(fp, " %zu", hist.count[j]);
			}
			fprintf(fp, "\n");
		}

		for (size_t i = 0; i<nn; i++){
			for (size_t j = 0; j<nn; j++){
				double cov = covariance(nuisance[i], nuisance[j]);
				fprintf(fp, " %15.6e", cov);
			}
			fprintf(fp, "\n");
		}

		for (size_t i = 0; i<nn; i++){
			for (size_t j = 0; j<nn; j++){
				double cor = correlation(nuisance[i], nuisance[j]);
				fprintf(fp, " %15.6e", cor);
			}
			fprintf(fp, "\n");
		}
				
		for (size_t i = 0; i<nn; i++){
			bwrite(fp, nuisance[i]);
		}

	}
};

class rjMcMC1DPPDMap{

private:	
	size_t nlmin;
	size_t nlmax;
	double pmax;
	double vmin;
	double vmax;
	size_t nsamples;//number of samples	
	size_t np;//number of positions
	size_t nv;//number of values
	double dp;
	double dv;	
	std::vector<double> pbin;//positions
	std::vector<double> vbin;//values
	std::vector<std::vector<uint32_t>> counts;//frequency
	std::vector<uint32_t> cpcounts;//changepoints
	std::vector<uint32_t> layercounts;//number of layers

public:

	size_t npbins(){ return np;}
	size_t nvbins(){ return nv;}
	double toppbin(size_t i){ return pbin[i] - dp/2.0;}
	const std::vector<uint32_t>& changepoint(){ return cpcounts; }

	inline size_t getvbin(double val){
		if (val < vmin)return 0;
		if (val >= vmax)return nv-1;				
		return (size_t)((val-vmin)/dv);
	}
	inline size_t getpbin(double pos){
		if (pos < 0.0)return 0;
		if (pos >= pmax)return np-1;
		return (size_t)(pos/dp);
	}
	void initialise(size_t _nlmin, size_t _nlmax, double _pmax, size_t _np, double _vmin, double _vmax, size_t _nv)
	{
		nlmin = _nlmin;
		nlmax = _nlmax;
		np = _np;
		nv = _nv;
		pmax = _pmax;
		vmin = _vmin;
		vmax = _vmax;

		dp = _pmax / (double)(np);
		dv = (vmax - vmin)/(double)(nv);
		
		pbin.resize(np);
		for (size_t i = 0; i<np; i++)pbin[i] = dp*((double)i+0.5);
		
		vbin.resize(nv);
		for (size_t i = 0; i<nv; i++)vbin[i] = vmin + dv*((double)i+0.5);
		resettozero();
	}
	void resettozero()
	{
		nsamples = 0;
		layercounts.resize(nlmax - nlmin + 1, 0);
		counts.resize(np);
		for (size_t i=0; i<np; i++){
			counts[i].assign(nv,0);			
		}
		cpcounts.assign(np,0);
	}
	
	void addmodel(const rjMcMC1DModel& m)
	{		
		nsamples++;
		
		layercounts[m.nlayers()-nlmin]++;

		for(size_t pi=0; pi<np; pi++){			
			size_t li = m.which_layer(pbin[pi]);
			size_t vi = getvbin(m.layers[li].value);
			counts[pi][vi]++;
		}

		for(size_t li=1; li<m.nlayers(); li++){			
			size_t pi = getpbin(m.layers[li].ptop);
			cpcounts[pi]++;
		}

	}
	
	std::vector<double> modelmap(const rjMcMC1DModel& m)
	{
		std::vector<double> model(np);		
		for (size_t pi=0;pi<np;pi++){
			size_t li = m.which_layer(pbin[pi]);			
			model[pi] = m.layers[li].value;
		}
		return model;
	}

	cHistogramStats<double> hstats(size_t pi){
		cHistogramStats<double> hs(vbin, counts[pi]);
		return hs;
	}

	std::vector<cHistogramStats<double>> hstats(){
		std::vector<cHistogramStats<double>> hs(np);
		for (size_t pi = 0; pi < np; pi++){
			hs[pi] = hstats(pi);
		}
		return hs;
	}

	void writedata(FILE* fp){				
		bwrite(fp, layercounts);
		bwrite(fp, cpcounts);
		for (size_t pi = 0; pi < np; pi++){
			bwrite(fp, counts[pi]);
		}
	}

};

class cChainInfo{

public:
	std::vector<rjMcMC1DModel> modelchain;
	std::vector<uint32_t> sample;
	std::vector<uint32_t> nlayers;
	std::vector<float> misfit;
	std::vector<float> logppd;
	std::vector<float> ar_valuechange;
	std::vector<float> ar_move;
	std::vector<float> ar_birth;
	std::vector<float> ar_death;
	std::vector<float> ar_nuisancechange;	

	void reset(){
		sample.clear();
		nlayers.clear();
		misfit.clear();
		logppd.clear();
		ar_valuechange.clear();
		ar_move.clear();
		ar_birth.clear();
		ar_death.clear();
		ar_nuisancechange.clear();
	}

	void writeconvergencedata(size_t ci, FILE* fp)
	{		
		bwrite(fp, (uint32_t)ci);
		bwrite(fp, (uint32_t)sample.size());
		bwrite(fp, sample);
		bwrite(fp, nlayers);
		bwrite(fp, misfit);
		bwrite(fp, logppd);
		bwrite(fp, ar_valuechange);
		bwrite(fp, ar_move);
		bwrite(fp, ar_birth);
		bwrite(fp, ar_death);
		bwrite(fp, ar_nuisancechange);
	}

	void writemodelchain_binary(size_t ci, FILE* fp)
	{
		bwrite(fp, (uint32_t)modelchain.size());
		for (size_t si = 0; si < modelchain.size(); si++){
			rjMcMC1DModel& m = modelchain[si];
			bwrite(fp, (uint32_t)ci);
			bwrite(fp, (uint32_t)si);
			bwrite(fp, (uint32_t)m.nlayers());
			bwrite(fp, (float)m.misfit());			
			for (size_t li = 0; li < m.nlayers(); li++){
				bwrite(fp, (float)m.layers[li].ptop);
			}			
			for (size_t li = 0; li<m.nlayers(); li++){
				bwrite(fp, (float)m.layers[li].value);				
			}			
		}
	}

	void writemodelchain_ascii(size_t ci, FILE* fp)
	{
		for (size_t si = 0; si < modelchain.size(); si++){
			rjMcMC1DModel& m = modelchain[si];

			fprintf(fp, "%zu %zu %zu %.6e\n", ci, si, m.nlayers(), m.misfit());
			for (size_t li = 0; li < m.nlayers(); li++){
				fprintf(fp, " %.6e", m.layers[li].ptop);
			}
			fprintf(fp, "\n");
			for (size_t li = 0; li<m.nlayers(); li++){
				fprintf(fp, " %.6e", m.layers[li].value);
			}
			fprintf(fp, "\n");//%end of line
		}
	}
};

class rjMcMC1DSampler{

public:
	
	double DEFAULTLOGSTDDECADES   = 0.05;
	double DEFAULTMOVESTDFRACTION = 0.05;
	size_t mRank;
	size_t mSize;	

	std::string description;
	std::string starttime;
	std::string endtime;
	double samplingtime;

	std::vector<rjMcMCNuisance> nuisance_init;
	size_t nl_min;
	size_t nl_max;	
	double  pmax;	
	double  vmin;
	double  vmax;	
	double sd_valuechange;//std deviation of gaussian proposal on value change	
	double sd_bd_valuechange;//std deviation of gaussian proposal on BIRTH % (this number is also present in the DEATH % acceptance term when taking in acount the reverse jump )
	double sd_move;//std deviation of gaussian proposal on move	
	
	size_t nchains;
	size_t nsamples;
	size_t nburnin;
	size_t thinrate;
	
	std::vector<cChainInfo>  mChainInfo;
	
	size_t np_valuechange;
	size_t np_birth;
	size_t np_death;
	size_t np_move;
	size_t np_nuisancechange;

	size_t na_valuechange;
	size_t na_birth;
	size_t na_death;
	size_t na_move;
	size_t na_nuisancechange;
	
	enum parameterization_type param_position;
	enum parameterization_type param_value;
	
	rjMcMC1DPPDMap pmap;
	rjMcMC1DNuisanceMap nmap;
	
	rjMcMC1DModel  mHighestLikelihood;
	rjMcMC1DModel  mLowestMisfit;
		
	bool reportstats;
	size_t ndata;		 	
	std::vector<double> obs;
	std::vector<double> err;
	bool mBirthDeathFromPrior;
	
	static double gaussian_pdf(double mean, double std, double x)
	{
		double p = exp(-0.5 * pow((x - mean) / std, 2.0)) / (sqrt(TWOPI)*std);
		return p;
	}
	
	void addmodel(const rjMcMC1DModel& m)
	{
		pmap.addmodel(m);
		nmap.addmodel(m);
	}

	double ar_valuechange(){
		if(np_valuechange==0)return -1;
		else return 100.0*(double)na_valuechange/(double)np_valuechange;
	}
	
	double ar_move(){
		if(np_move==0)return -1;
		return 100.0*(double)na_move/(double)np_move;
	}
	
	double ar_birth(){
		if(np_birth==0)return -1;
		return 100.0*(double)na_birth/(double)np_birth;
	}
	
	double ar_death(){
		if(np_death==0)return -1;
		return 100.0*(double)na_death/(double)np_death;
	}	
	
	double ar_nuisancechange(){
		if(np_nuisancechange==0)return -1;
		return 100.0*(double)na_nuisancechange/(double)np_nuisancechange;
	}
	
	virtual std::vector<double> forwardmodel(const rjMcMC1DModel& m) = 0;		
	
	void reset()
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
	
	bool isreportable(size_t si)
	{
		if (si == 0 || si == nsamples - 1)return true;
		size_t k = (size_t)pow(10, floor(log10((double)si)));
		if (k > 1000)k = 1000;
		if (si % k == 0)return true;
		return false;
	}
	
	void printstats(size_t ci, size_t si, size_t np, double nmf)
	{
		printf("ci=%4zu si=%7zu np=%2zu nmf=%10g vc=%5.2lf m=%5.2lf b=%5.2lf d=%5.2lf n=%5.2lf\n", ci, si, np, nmf, ar_valuechange(), ar_move(), ar_birth(), ar_death(), ar_nuisancechange());
	}
	
	size_t nnuisances()
	{
		return nuisance_init.size();
	}
	
	bool isinbounds(double& bmin, double& bmax, double& b)
	{
		if (b<bmin || b>bmax)return false;
		else return true;
	}
	
	bool includeinmap(size_t si)
	{
		if (si < nburnin)return false;
		size_t k = si - nburnin;
		if ((k % thinrate) == 0) {
			return true;
		}
		return false;
	}
	
	double l2misfit(const std::vector<double>& g)
	{
		double sum = 0.0;
		for (size_t di = 0; di < ndata; di++) {
			double nr = (obs[di] - g[di]) / err[di];
			sum += nr * nr;
		}
		return sum;
	}
	
	double computemisfit(const rjMcMC1DModel& m)
	{
		std::vector<double> pred = forwardmodel(m);
		return l2misfit(pred);
	}

	void set_misfit(rjMcMC1DModel& m)
	{
		double mfit = computemisfit(m);
		//double mfit = (double)ndata;
		m.setmisfit(mfit);
	}
	
	void writeresults(FILE* fp)
	{
		writeheader(fp);
		writemodel(fp, mHighestLikelihood);
		writemodel(fp, mLowestMisfit);
		pmap.writedata(fp);
		nmap.writedata(fp);
		for (size_t ci = 0; ci < nchains; ci++) {
			mChainInfo[ci].writeconvergencedata(ci, fp);
			mChainInfo[ci].writemodelchain_binary(ci, fp);
		}
	}
	
	void writeheader(FILE* fp)
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
	
	void writemodel(FILE* fp, const rjMcMC1DModel& m)
	{
		fprintf(fp, "%.6e %.6e %zu", m.misfit(), m.logppd(), m.nlayers());
		for (size_t li = 0; li < m.nlayers(); li++) {
			fprintf(fp, " %.6e", m.layers[li].ptop);
		}
		for (size_t li = 0; li < m.nlayers(); li++) {
			fprintf(fp, " %.6e", m.layers[li].value);
		}
		fprintf(fp, "\n");
	}
	
	bool saveconvergencerecord(size_t si)
	{
		if (si == 0 || si == nsamples - 1)return true;
		size_t k = (size_t)pow(10, floor(log10((double)si)));
		if (k > thinrate) k = thinrate;
		if (si % k == 0)return true;
		return false;
	}

	bool propose_valuechange(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
	{
		np_valuechange++;

		size_t index = irand((size_t)0, (size_t)(mcur.nlayers() - 1));

		double logstd = DEFAULTLOGSTDDECADES;
		double vold = mcur.layers[index].value;
		double vnew;
		double pqratio;
		if (param_value == LINEAR) {
			double m = (std::pow(10.0, logstd) - std::pow(10.0, -logstd)) / 2.0;
			vnew = vold + m * vold * nrand<double>();
			double qpdfforward = gaussian_pdf(vold, m * vold, vnew);
			double qpdfreverse = gaussian_pdf(vnew, m * vnew, vold);
			pqratio = qpdfreverse / qpdfforward;
		}
		else {
			vnew = vold + logstd * nrand<double>();
			pqratio = 1.0;
		}

		bool isvalid = isinbounds(vmin, vmax, vnew);
		if (isvalid == false)return false;
		mpro.layers[index].value = vnew;

		set_misfit(mpro);
		double logpqratio = std::log(pqratio);
		double logliker = -(mpro.misfit() - mcur.misfit()) / 2.0;
		double logar = logpqratio + logliker;
		if (std::log(urand<double>()) < logar) {
			na_valuechange++;
			return true;
		}
		return false;
	}

	bool propose_move(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
	{
		np_move++;

		size_t n = mcur.nlayers();
		if (n <= 1)return false;

		size_t index = irand((size_t)1, n - 1);
		double pold = mcur.layers[index].ptop;

		//double std = sd_move;	
		double std = DEFAULTMOVESTDFRACTION * pold;
		double pnew = pold + std * nrand<double>();
		double qpdfforward = gaussian_pdf(pold, pold * DEFAULTMOVESTDFRACTION, pnew);
		double qpdfreverse = gaussian_pdf(pnew, pnew * DEFAULTMOVESTDFRACTION, pold);


		bool isvalid = mpro.move_interface(index, pnew);
		if (isvalid == false)return false;

		set_misfit(mpro);
		double pqratio = qpdfreverse / qpdfforward;
		double logpqratio = std::log(pqratio);

		double loglr = -(mpro.misfit() - mcur.misfit()) / 2.0;
		double logar = logpqratio + loglr;
		double logu = std::log(urand<double>());
		if (logu < logar) {
			na_move++;
			return true;
		}
		return false;
	}
	
	bool propose_birth(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
	{
		np_birth++;

		size_t n = mcur.nlayers();
		if (n >= nl_max)return false;

		double  pos = urand(0.0, pmax);
		size_t  index = mcur.which_layer(pos);
		double vold = mcur.layers[index].value;
		double vnew, pqratio;

		if (mBirthDeathFromPrior) {
			vnew = urand(vmin, vmax);
			pqratio = 1.0;
		}
		else {
			double vcpdf;
			double logstd = DEFAULTLOGSTDDECADES;
			if (param_value == LINEAR) {
				double m = (pow(10.0, logstd) - pow(10.0, -logstd)) / 2.0;
				vnew = vold + m * vold * nrand<double>();
				vcpdf = gaussian_pdf(vold, m * vold, vnew);
			}
			else {
				vnew = vold + logstd * nrand<double>();
				vcpdf = gaussian_pdf(vold, logstd, vnew);
			}
			pqratio = 1.0 / ((vmax - vmin) * vcpdf);
		}
		//pqratio *= (double)(n) / double(n + 1);

		bool   isvalid = mpro.insert_interface(pos, vnew);
		if (isvalid == false)return false;
		set_misfit(mpro);

		double logpqratio = std::log(pqratio);
		double loglikeratio = -(mpro.misfit() - mcur.misfit()) / 2.0;
		double logar = logpqratio + loglikeratio;
		if (std::log(urand<double>()) < logar) {
			na_birth++;
			return true;
		}
		return false;
	}
	
	bool propose_death(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
	{
		np_death++;

		size_t n = mcur.nlayers();
		if (n <= nl_min)return false;

		size_t index = irand((size_t)1, n - 1);
		bool isvalid = mpro.delete_interface(index);
		if (isvalid == false)return false;
		set_misfit(mpro);

		double pqratio;
		if (mBirthDeathFromPrior) {
			pqratio = 1.0;
		}
		else {
			double logstd = DEFAULTLOGSTDDECADES;
			double vnew = mcur.layers[index - 1].value;
			double vold = mcur.layers[index].value;
			double vcpdf;
			if (param_value == LINEAR) {
				double m = (pow(10.0, logstd) - pow(10.0, -logstd)) / 2.0;
				vcpdf = gaussian_pdf(vnew, m * vnew, vold);
			}
			else {
				vcpdf = gaussian_pdf(vnew, logstd, vold);
			}
			pqratio = (vmax - vmin) * vcpdf;
		}
		//pqratio *= (double)(n) / double(n - 1);

		double logpqratio = log(pqratio);
		double loglikeratio = -(mpro.misfit() - mcur.misfit()) / 2.0;
		double logar = logpqratio + loglikeratio;
		if (log(urand<double>()) < logar) {
			na_death++;
			return true;
		}
		return false;
	}
	
	bool propose_nuisancechange(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
	{
		np_nuisancechange++;

		size_t ni = irand((size_t)0, mcur.nnuisances() - 1);
		double delta = nrand<double>() * mcur.nuisances[ni].sd_valuechange;
		double nv = mcur.nuisances[ni].value + delta;;
		bool isvalid = isinbounds(mcur.nuisances[ni].min, mcur.nuisances[ni].max, nv);
		if (isvalid == false)return false;

		mpro.nuisances[ni].value = nv;

		set_misfit(mpro);
		double logar = -(mpro.misfit() - mcur.misfit()) / 2.0;
		double logu = log(urand<double>());
		if (logu < logar) {
			na_nuisancechange++;
			return true;
		}
		return false;
	}
	
	bool propose_independent(rjMcMC1DModel& mcur, rjMcMC1DModel& mpro)
	{
		mpro = choosefromprior();
		set_misfit(mpro);

		size_t npro = mpro.nlayers() - 1;
		size_t ncur = mcur.nlayers() - 1;
		double a = pow(1.0 / (vmax - vmin), npro - ncur);
		double b = pow(1.0 / pmax, npro - ncur);
		double c = (double)factorial((unsigned int)npro) / (double)factorial((unsigned int)ncur);
		double priorratio = a * b * c;
		priorratio = 1.0;

		double logpriorratio = std::log(priorratio);
		double logliker = -(mpro.misfit() - mcur.misfit()) / 2.0;
		double logar = logpriorratio + logliker;
		double logu = std::log(urand<double>());
		if (logu < logar) {
			return true;
		}
		return false;
	}

	rjMcMC1DModel choosefromprior()
	{
		rjMcMC1DModel m;
		size_t nl = irand(nl_min, nl_max);
		m.initialise(pmax, vmin, vmax);
		for (size_t li = 0; li < nl; li++) {
			bool status = false;
			while (status == false) {
				double pos = urand(0.0, pmax);//position of interface		
				double value = urand(vmin, vmax);//value
				status = m.insert_interface(pos, value);
			}
		}
		m.nuisances = nuisance_init;
		return m;
	}

	void sample()
	{
		mBirthDeathFromPrior = false;

		starttime = timestamp();
		double t1 = gettime();
		mChainInfo.resize(nchains);
		for (size_t ci = 0; ci < nchains; ci++) {
			mChainInfo[ci].reset();
			//mChainInfo[ci].modelchain.resize(nsamples);
			rjMcMC1DModel mcur;
			for (size_t si = 0; si < nsamples; si++) {

				//Initialis chain
				if (si == 0) {
					mcur = choosefromprior();
					set_misfit(mcur);
				}

				//Initialise "best" models
				if (ci == 0 && si == 0) {
					mHighestLikelihood = mcur;
					mLowestMisfit = mcur;
				}

				rjMcMC1DModel mpro = mcur;
				size_t nopt = 4;
				if (mcur.nnuisances() > 0)nopt = 5;
				size_t option = irand((size_t)0, nopt - 1);

				bool accept = false;
				if (option == 0) {
					accept = propose_valuechange(mcur, mpro);
				}
				else if (option == 1) {
					accept = propose_move(mcur, mpro);
				}
				else if (option == 2) {
					accept = propose_birth(mcur, mpro);
				}
				else if (option == 3) {
					accept = propose_death(mcur, mpro);
				}
				else if (option == 4) {
					accept = propose_nuisancechange(mcur, mpro);
				}
				else if (option == 5) {
					accept = propose_independent(mcur, mpro);
				}
				else {
					exit(1);
					break;
				}

				if (accept) mcur = mpro;

				if (includeinmap(si)) {
					pmap.addmodel(mcur);
					nmap.addmodel(mcur);
					mChainInfo[ci].modelchain.push_back(mcur);
				}

				if (mcur.logppd() > mHighestLikelihood.logppd()) {
					mHighestLikelihood = mcur;
				}

				if (mcur.misfit() < mLowestMisfit.misfit()) {
					mLowestMisfit = mcur;
				}

				if (saveconvergencerecord(si)) {
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
				if (reportstats) {
					if (mRank == 0 && isreportable(si)) {
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
};


#endif
