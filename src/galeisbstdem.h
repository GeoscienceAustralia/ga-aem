/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _galeisbs_H
#define _galeisbs_H

#include <stdio.h>
#include <vector>
#include <cstring>
#include <algorithm>

#include <Eigen/Dense>
#include "inputmanager.h"
#include "outputmanager.h"
#include "blocklanguage.h"
#include "airborne_types.h"
#include "tdemsystem.h"
#include "eigen_utils.h"
#include "gaaem_version.h"

#if defined _OPENMP
	#include <omp.h>
	//This external thread lock must be set when fftw is being initialised
	extern omp_lock_t fftw_thread_lock;
#endif

enum class eNormType { L1, L2 };
enum class eSmoothnessMethod { DERIVATIVE_1ST, DERIVATIVE_2ND };
enum class eBracketResult { BRACKETED, MINBRACKETED, ALLABOVE, ALLBELOW };

class cTrial{

public:
	size_t order;
	double lambda;
	double phid;
	double phim;
	double stepfactor;

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

class cTrialCache{	

	public:

	double target;	
	std::vector<cTrial> trial;
	double sfsearch(const double& x){		
		for(size_t k=0; k<trial.size(); ++k){
			if(x == trial[k].stepfactor)return trial[k].phid;
		}		
		return -1.0;
	}
	
	size_t minphidindex(){
		size_t ind=0;		
		for(size_t k=1; k<trial.size(); ++k){						
			if( trial[k].phid < trial[ind].phid){
				ind=k;
			}		
		}		
		return ind;
	}

	size_t maxphidindex(){
		size_t ind=0;
		for(size_t k=1; k<trial.size(); ++k){			
			if( trial[k].phid > trial[ind].phid){
				ind=k;
			}			
		}	
		return ind;
	}

	size_t maxlambdaindex(){
		size_t ind=0;		
		for(size_t k=1; k<trial.size(); ++k){			
			if( trial[k].lambda > trial[ind].lambda){
				ind=k;
			}			
		}	
		return ind;
	}

	cTrial findlambda(const double lambda){
		for(size_t i=0; i<trial.size(); i++){
			if(trial[i].lambda == lambda){
			  return trial[i];			  
			}
		}
		return trial[0];
	}

	double minphid(){return trial[minphidindex()].phid;}
	double maxphid(){return trial[maxphidindex()].phid;}
	cTrial minphidtrial(){return trial[minphidindex()];}
	cTrial maxphidtrial(){return trial[maxphidindex()];}
	cTrial maxlambdatrial(){return trial[maxlambdaindex()];}

	void sort_lambda()
	{
		std::sort(trial.begin(), trial.end(), cTrial::lambda_compare);
	}

	void sort_phid()
	{
		std::sort(trial.begin(), trial.end(), cTrial::phid_compare);
	}

	void print(const double& lambda, const double& phid)
	{
		sort_lambda();
		printf("\n");
		printf("CurrentLambda = %lf CurrentPhid = %lf    Target = %lf\n", lambda, phid, target);
		printf("N    Stepfactor       Lambda          Phid\n");
		for (size_t i = 0; i < trial.size(); i++) {
			printf("%2zu %12g %12g %12g\n", trial[i].order, trial[i].stepfactor, trial[i].lambda, trial[i].phid);
		}
		printf("\n");
	}
};

class cOutputOptions {

private:
	std::string DumpBasePath;

public:	
	std::string LogFile;
	bool PositiveLayerTopDepths = false;
	bool NegativeLayerTopDepths = false;
	bool PositiveLayerBottomDepths = false;
	bool NegativeLayerBottomDepths = false;
	bool InterfaceElevations = false;
	bool ParameterSensitivity = false;
	bool ParameterUncertainty = false;
	bool ObservedData = false;
	bool NoiseEstimates = false;
	bool PredictedData = false;
	bool Dump = false;
	std::string DumpPath(const size_t datafilerecord, const size_t iteration) const
	{
		return DumpBasePath + pathseparatorstring() + 
			strprint("si%07d", (int)datafilerecord) + pathseparatorstring() +
			strprint("it%03d", (int)iteration) + pathseparatorstring();
	};
	cOutputOptions(){};
	cOutputOptions(const cBlock& b) {		
		LogFile = b.getstringvalue("LogFile");				
		fixseparator(LogFile);

		PositiveLayerTopDepths = b.getboolvalue("PositiveLayerTopDepths");
		NegativeLayerTopDepths = b.getboolvalue("NegativeLayerTopDepths");
		PositiveLayerBottomDepths = b.getboolvalue("PositiveLayerBottomDepths");
		NegativeLayerBottomDepths = b.getboolvalue("NegativeLayerBottomDepths");
		InterfaceElevations = b.getboolvalue("InterfaceElevations");
		ParameterSensitivity = b.getboolvalue("ParameterSensitivity");
		ParameterUncertainty = b.getboolvalue("ParameterUncertainty");
		ObservedData = b.getboolvalue("ObservedData");
		NoiseEstimates = b.getboolvalue("NoiseEstimates");
		PredictedData = b.getboolvalue("PredictedData");				

		Dump = b.getboolvalue("Dump");
		if (Dump) {
			DumpBasePath = b.getstringvalue("DumpPath");
			fixseparator(DumpBasePath);
			if (DumpBasePath[DumpBasePath.length() - 1] != pathseparator()) {
				DumpBasePath.append(pathseparatorstring());
			}
			makedirectorydeep(DumpBasePath.c_str());
		}
	}	
};

class cComponentInfo {

public:

	std::string Name;
	bool Use = false;
	double  oP;
	std::vector<double>  oS;
	std::vector<double>  oE;
	cFieldDefinition fd_oP;
	cFieldDefinition fd_oS;
	cFieldDefinition fd_oE;
	bool EstimateNoiseFromModel;
	std::vector<double> mn;
	std::vector<double> an;
	int dataindex = 0;

	cComponentInfo() {};

	cComponentInfo(const cBlock& b, std::string name, size_t nwindows, bool inverttotalfield)
	{
		Name = name;
		if (b.Entries.size() == 0) {
			Use = false;
			return;
		}
		Use = b.getboolvalue("Use");
		if (Use == false)return;

		EstimateNoiseFromModel = b.getboolvalue("EstimateNoiseFromModel");

		if (EstimateNoiseFromModel) {
			mn = b.getdoublevector("MultiplicativeNoise");
			an = b.getdoublevector("AdditiveNoise");
			if (an.size() == 1) {
				an = std::vector<double>(nwindows, an[0]);
			}
			if (an.size() != nwindows) {				
				glog.errormsg(_SRC_,"Must have exactly 1 or nwindows AdditiveNoise values\n");				
			};

			if (mn.size() == 1) {
				mn = std::vector<double>(nwindows, mn[0]);
			}			
			if (mn.size() != nwindows) {
				glog.errormsg(_SRC_,"Must have exactly 1 or nwindows MultiplicativeNoise values\n");				
			}
		}

		oS.resize(nwindows);
		oE.resize(nwindows);
		fd_oP.initialise(b, "Primary");
		fd_oS.initialise(b, "Secondary");
		fd_oE.initialise(b, "Noise");
	}

	size_t nw() const
	{
		return oS.size();
	}

	void readdata(const std::unique_ptr<cInputManager>& IM)
	{
		if (Use == false) return;
		IM->read(fd_oP, oP);
		IM->read(fd_oS, oS, nw());
		oE.resize(nw());
		if (EstimateNoiseFromModel) {
			for (size_t w = 0; w < nw(); w++) {
				const double v = 0.01*mn[w]*oS[w];
				oE[w] = std::hypot(an[w],v);
			}
		}
		else {
			IM->read(fd_oE, oE, nw());
		}
	}
};

class cTDEmSystemInfo {

public:

	cTDEmSystem T;
	std::string SystemFile;
	size_t nwindows;
	size_t ncomps;
	size_t nchans;
	cComponentInfo CompInfo[3];
	int xzIndex;//Start index of XZ in data array

	bool invertXPlusZ;
	bool invertPrimaryPlusSecondary;
	bool reconstructPrimary;		
	cTDEmData predicted;

	void initialise(const cBlock& b) {
		std::string dummy;
		if (b.getvalue("InvertTotalField", dummy)) {
			glog.errormsg(_SRC_,"InvertTotalField is no longer an option, use InvertPrimaryPlusSecondary instead\n");
		};

		std::string stmfile = b.getstringvalue("SystemFile");
		fixseparator(stmfile);

		glog.logmsg(0, "Reading system file %s\n", stmfile.c_str());
		T.readsystemdescriptorfile(stmfile);
		glog.log("==============System file %s\n", stmfile.c_str());
		glog.log(T.STM.get_as_string());
		glog.log("==========================================================================\n");
		nwindows = T.NumberOfWindows;

		invertXPlusZ = b.getboolvalue("InvertXPlusZ");
		invertPrimaryPlusSecondary = b.getboolvalue("InvertPrimaryPlusSecondary");
		reconstructPrimary = b.getboolvalue("ReconstructPrimaryFieldFromInputGeometry");

		CompInfo[0] = cComponentInfo(b.findblock("XComponent"), "X", nwindows, invertPrimaryPlusSecondary);
		CompInfo[1] = cComponentInfo(b.findblock("YComponent"), "Y", nwindows, invertPrimaryPlusSecondary);
		CompInfo[2] = cComponentInfo(b.findblock("ZComponent"), "Z", nwindows, invertPrimaryPlusSecondary);

		ncomps = 0;
		if (CompInfo[0].Use) ncomps++;
		if (CompInfo[1].Use) ncomps++;
		if (CompInfo[2].Use) ncomps++;

		if (invertXPlusZ) {
			CompInfo[0].Use = true;
			CompInfo[2].Use = true;
		}
		nchans = nwindows * ncomps;
	}
};

class cSBSInverter{

	private:
	std::vector<size_t> Q;//indices of the data that are not culled
	std::string CommandLine;
	int Size;
	int Rank;
	bool UsingOpenMP;

	public:		
	std::string OutputMessage;
	cBlock Control;
	
	std::vector<cTDEmSystemInfo> SV;
		
	size_t nsystems;	
	
	std::unique_ptr<cInputManager> IM;
	std::unique_ptr<cOutputManager> OM;
	size_t pointsoutput = 0;

	cOutputOptions OO;	
	
	//Column definitions
	cFieldDefinition sn, dn, fn, ln, fidn;
	cFieldDefinition xord, yord, elevation;	
	std::vector<cFieldDefinition> fd_GI;
	std::vector<cFieldDefinition> fd_GR;
	std::vector<cFieldDefinition> fd_GS;	
	std::vector<cFieldDefinition> fd_GTFR;	
	cFieldDefinition fd_ERc;
	cFieldDefinition fd_ERt;
	cFieldDefinition fd_ESc;
	cFieldDefinition fd_ESt;		
			
	sAirborneSampleId Id;
	sAirborneSampleLocation Location; 	
	cTDEmGeometry GI;//Input geometry	
	cTDEmGeometry GR;//Reference model geometry
	cTDEmGeometry GS;//Standard deviation geometry
	cTDEmGeometry GTFR;//Total field reconstruction geometry
	cTDEmGeometry GM;//Final inversion geometry
	
	cEarth1D ER;//Reference model earth
	cEarth1D ES;//Standard deviation earth
	cEarth1D EM;//Final inversion earth

	double min_conductivity = 1e-5;
	double max_conductivity = 10.0;

	bool solve_conductivity;
	bool solve_thickness;	

	bool solve_tx_height;
	bool solve_tx_roll;
	bool solve_tx_pitch;
	bool solve_tx_yaw;
	bool solve_txrx_dx;
	bool solve_txrx_dy;
	bool solve_txrx_dz;
	bool solve_rx_roll;
	bool solve_rx_pitch;
	bool solve_rx_yaw;	
	
	double AlphaC;
	double AlphaT;
	double AlphaG;
	double AlphaS;
	eSmoothnessMethod SmoothnessMethod;
	eNormType  NormType;

	size_t ndata_all;
	size_t ndata_inv;
	size_t nparam;
	size_t nlayers;	
	size_t ngeomparam;	
	
	VectorDouble Obs;
	VectorDouble Err;
	VectorDouble Pred;

	VectorDouble Param;
	VectorDouble RefParam;
	VectorDouble RefParamStd;

	VectorDouble ParameterSensitivity;
	VectorDouble ParameterUncertainty;
	
	MatrixDouble J;		
	MatrixDouble Wd;
	MatrixDouble Wc;
	MatrixDouble Wt;
	MatrixDouble Wg;
	MatrixDouble Wr;	
	MatrixDouble Ws;
	MatrixDouble Wm;
		
	double MinimumPhiD;//overall	
	double MinimumImprovement;//			
	size_t MaxIterations;
	std::string TerminationReason;

	double TargetPhiD;//for an iteration
	
	double LastPhiD;//for previous iteration
	double LastPhiM;
	double LastPhiC;
	double LastPhiT;
	double LastPhiG;
	double LastPhiS;
	double LastLambda;
	size_t LastIteration;
	
	
	size_t cIndex;//Conductivty parameters start index
    size_t tIndex;//Thickness parameters start index
	size_t gIndex;//Geometery parameters start index
	size_t tx_heightIndex;
	size_t txrx_dxIndex;
	size_t txrx_dyIndex;
	size_t txrx_dzIndex;
	size_t rx_pitchIndex;
	size_t rx_rollIndex;
	
	cSBSInverter(const std::string& controlfile, const int& size, const int& rank, const bool& usingopenmp, const std::string commandline)
	{		
		_GSTPUSH_
		try {					
			Size = size;
			Rank = rank;
			UsingOpenMP = usingopenmp;
			CommandLine = commandline;
			initialise(controlfile);			
		}
		catch (const std::string& msg) {
			std::cerr << msg;
			glog.logmsg(msg);			
		}
		catch (const std::runtime_error& e) {
			std::cerr << e.what();
			glog.logmsg(std::string(e.what()));
		}
		catch (const std::exception& e) {
			std::cerr << e.what();
			glog.logmsg(std::string(e.what()));
		}		
	};

	~cSBSInverter()
	{				
		glog.close();
	};

	void initialise(const std::string& controlfile)
	{
		_GSTITEM_
		loadcontrolfile(controlfile);
		set_field_definitions();
		setup_data();
		setup_parameters();		
		go();	
	}

	void loadcontrolfile(const std::string& filename)
	{
		glog.logmsg(0, "Loading control file %s\n", filename.c_str());
		Control = cBlock(filename);
		OO = cOutputOptions(Control.findblock("Output"));
		std::string suffix = stringvalue(Rank, ".%04d");				
		OO.LogFile = insert_after_filename(OO.LogFile,suffix);
		openlogfile(); //load this first to get outputlogfile opened

		//Load control file
		parseoptions();
		initialisesystems();

		cBlock ib = Control.findblock("Input");
		if (cInputManager::isnetcdf(ib)){
			#if !defined HAVE_NETCDF
			glog.errormsg(_SRC_, "Sorry NETCDF I/O is not available in this executable\n");
			#endif			
			IM = std::make_unique<cNetCDFInputManager>(ib);			
			std::string s  = IM->datafilename();
		}
		else {			
			IM = std::make_unique<cASCIIInputManager>(ib);
		}

		cBlock ob = Control.findblock("Output");
		if (cOutputManager::isnetcdf(ob)){
			#if !defined HAVE_NETCDF
			glog.errormsg(_SRC_, "Sorry NETCDF I/O is not available in this executable\n");
			#endif			
			OM = std::make_unique<cNetCDFOutputManager>(ob,Size,Rank);			
		}
		else {
			OM = std::make_unique<cASCIIOutputManager>(ob,Size,Rank);
		}				
		OM->opendatafile(IM->datafilename(), IM->subsamplerate());		
	}

	void openlogfile()
	{		
		glog.logmsg(0, "Opening log file %s\n", OO.LogFile.c_str());
		glog.open(OO.LogFile);
		glog.logmsg(0, "%s\n", CommandLine.c_str());
		glog.logmsg(0, "%s\n", versionstring(GAAEM_VERSION, __TIME__, __DATE__).c_str());				
		glog.logmsg(0, "Working directory %s\n", getcurrentdirectory().c_str());
		if (UsingOpenMP && Size > 1) {
			glog.logmsg(0, "Using OpenMP threading Processes=%d\tRank=%d\n", Size, Rank);
		}
		else if (Size>1) {
			glog.logmsg(0, "Using MPI Processes=%d\tRank=%d\n", Size, Rank);
		}
		else {
			glog.logmsg(0, "Standalone Processes=%d\tRank=%d\n", Size, Rank);
		}

		glog.logmsg(0, "Control file %s\n", Control.Filename.c_str());
		glog.log(Control.get_as_string());
		glog.flush();
	}

	int go()
	{
		_GSTITEM_
			int job = 0;
		while (IM->readnextrecord()) {
			job++;
			if (((job - 1) % Size) != Rank) continue;
			bool valid = IM->is_record_valid();
			if (valid == false) continue;

			if (OO.Dump) {
				FILE* fp = fileopen(dumppath() + "record.dat", "w");
				fprintf(fp, "Record\t%zu", IM->record());
				fclose(fp);
			}

			double t1 = gettime();
			if (invert()) {
				double t2 = gettime();
				double etime = t2 - t1;
				writeresult(job - 1);
				std::string msg = strprint("Rec %6zu  %3zu  %5zu  %10lf  Its=%3zu  PhiD=%6.2lf  time=%.1lfs  %s %s\n", 1 + IM->record(), Id.flightnumber, Id.linenumber, Id.fidnumber, LastIteration, LastPhiD, etime, TerminationReason.c_str(), OutputMessage.c_str());
				glog.logmsg(msg);
				if (OutputMessage.size() > 0) {
					std::cerr << msg;
				}
			}
			else {
				std::string msg = strprint("Rec %6zu  Skipping %s\n", 1 + IM->record(), OutputMessage.c_str());
				glog.logmsg(msg);
				std::cerr << msg;
			}
			}
		glog.close();
		return 0;
		}

	void parseoptions()
	{
		cBlock b = Control.findblock("Options");
		solve_conductivity = b.getboolvalue("SolveConductivity");
		solve_thickness = b.getboolvalue("SolveThickness");

		solve_tx_height = b.getboolvalue("SolveTX_Height");
		solve_tx_roll = b.getboolvalue("SolveTX_Roll");
		solve_tx_pitch = b.getboolvalue("SolveTX_Pitch");
		solve_tx_yaw = b.getboolvalue("SolveTX_Yaw");
		solve_txrx_dx = b.getboolvalue("SolveTXRX_DX");
		solve_txrx_dy = b.getboolvalue("SolveTXRX_DY");
		solve_txrx_dz = b.getboolvalue("SolveTXRX_DZ");
		solve_rx_roll = b.getboolvalue("SolveRX_Roll");
		solve_rx_pitch = b.getboolvalue("SolveRX_Pitch");
		solve_rx_yaw = b.getboolvalue("SolveRX_Yaw");

		AlphaC = b.getdoublevalue("AlphaConductivity");
		AlphaT = b.getdoublevalue("AlphaThickness");
		AlphaG = b.getdoublevalue("AlphaGeometry");
		AlphaS = b.getdoublevalue("AlphaSmoothness");


		NormType = eNormType::L2;//default
		std::string nt = b.getstringvalue("NormType");
		if (!isdefined(nt)) {
			NormType = eNormType::L2;
		}
		else if (strcasecmp(nt, "L1") == 0) {
			NormType = eNormType::L1;
		}
		else if (strcasecmp(nt, "L2") == 0) {
			NormType = eNormType::L2;
		}
		else {
			glog.errormsg("Unknown NormType %s\n", nt.c_str());
		}


		SmoothnessMethod = eSmoothnessMethod::DERIVATIVE_2ND;//default
		std::string sm = b.getstringvalue("SmoothnessMethod");
		if (!isdefined(sm)) {
			SmoothnessMethod = eSmoothnessMethod::DERIVATIVE_2ND;
		}
		else if (strcasecmp(sm, "Minimise1stDerivatives") == 0) {
			SmoothnessMethod = eSmoothnessMethod::DERIVATIVE_1ST;
		}
		else if (strcasecmp(sm, "Minimize1stDerivatives") == 0) {
			SmoothnessMethod = eSmoothnessMethod::DERIVATIVE_1ST;
		}
		else if (strcasecmp(sm, "Minimise2ndDerivatives") == 0) {
			SmoothnessMethod = eSmoothnessMethod::DERIVATIVE_2ND;
		}
		else if (strcasecmp(sm, "Minimize2ndDerivatives") == 0) {
			SmoothnessMethod = eSmoothnessMethod::DERIVATIVE_2ND;
		}
		else {
			glog.errormsg(_SRC_, "Unknown SmoothnessMethod %s\n", sm.c_str());
		}


		MaxIterations = b.getsizetvalue("MaximumIterations");
		MinimumPhiD = b.getdoublevalue("MinimumPhiD");
		MinimumImprovement = b.getdoublevalue("MinimumPercentageImprovement");
	}

	void set_field_definitions()
	{
		cBlock b = Control.findblock("Input.Columns");
		sn.initialise(b, "SurveyNumber");
		dn.initialise(b, "DateNumber");
		fn.initialise(b, "FlightNumber");
		ln.initialise(b, "LineNumber");
		fidn.initialise(b, "FidNumber");
		xord.initialise(b, "Easting");
		yord.initialise(b, "Northing");
		elevation.initialise(b, "GroundElevation");

		fd_GI = set_field_definitions_geometry(b);

		cBlock rm = b.findblock("ReferenceModel");
		fd_GR = set_field_definitions_geometry(rm);
		fd_ERc.initialise(rm, "Conductivity");
		fd_ERt.initialise(rm, "Thickness");

		cBlock sd = b.findblock("StdDevReferenceModel");
		fd_GS = set_field_definitions_geometry(sd);
		fd_ESc.initialise(sd, "Conductivity");
		fd_ESt.initialise(sd, "Thickness");

		cBlock tfr = b.findblock("TotalFieldReconstruction");
		fd_GTFR = set_field_definitions_geometry(tfr);		
	}

	std::vector<cFieldDefinition> set_field_definitions_geometry(const cBlock& b)
	{
		std::vector<cFieldDefinition> g(10);
		for (size_t i = 0; i < g.size(); i++) {
			g[i].initialise(b, cTDEmGeometry::fname(i));
		}
		return g;
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
		if(UsingOpenMP) {
			omp_unset_lock(&fftw_thread_lock);		
		}
		#endif	
	}

	void initialisesystems()
	{
		set_fftw_lock();				
		std::vector<cBlock> B = Control.findblocks("EMSystem");
		nsystems = B.size();
		SV.resize(nsystems);
		for (size_t si = 0; si < nsystems; si++) {
			SV[si].initialise(B[si]);
		}
		unset_fftw_lock();
	}

	void setup_data()
	{
		ndata_all = 0;
		for (size_t si = 0; si < nsystems; si++) {
			cTDEmSystemInfo& S = SV[si];
			if (S.invertXPlusZ) {
				S.xzIndex = (int) ndata_all;
				S.CompInfo[0].dataindex = -1;
				S.CompInfo[2].dataindex = -1;
				ndata_all += S.nwindows;

				if (S.CompInfo[1].Use) {
					S.CompInfo[1].dataindex = (int) ndata_all;
					ndata_all += S.nwindows;
				}
			}
			else {
				for (size_t i = 0; i < 3; i++) {
					if (S.CompInfo[i].Use) {
						S.CompInfo[i].dataindex = (int) ndata_all;
						ndata_all += S.nwindows;
					}
				}
			}
		}				
	}

	void setup_parameters()
	{
		double v;
		if (Control.getvalue("Earth.MinConductivity", v)) {
			min_conductivity = v;
		}

		if (Control.getvalue("Earth.MaxConductivity", v)) {
			max_conductivity = v;
		}

		nlayers = Control.getsizetvalue("Earth.NumberOfLayers");
		nparam = 0;
		ngeomparam = 0;
		if (solve_conductivity) {
			cIndex = nparam;
			nparam += nlayers;
		}

		if (solve_thickness) {
			tIndex = nparam;
			nparam += nlayers - 1;
		}


		//Geometry params
		gIndex = nparam;
		if (solve_tx_height) {
			tx_heightIndex = nparam;
			nparam++;
			ngeomparam++;
		}

		if (solve_txrx_dx) {
			txrx_dxIndex = nparam;
			nparam++;
			ngeomparam++;
		}

		if (solve_txrx_dy) {
			txrx_dyIndex = nparam;
			nparam++;
			ngeomparam++;
		}

		if (solve_txrx_dz) {
			txrx_dzIndex = nparam;
			nparam++;
			ngeomparam++;
		}

		if (solve_rx_pitch) {
			rx_pitchIndex = nparam;
			nparam++;
			ngeomparam++;
		}

		if (solve_rx_roll) {
			rx_rollIndex = nparam;
			nparam++;
			ngeomparam++;
		}

		///////////////
		Param.resize(nparam);
		RefParam.resize(nparam);
		RefParamStd.resize(nparam);
	}

	bool parserecord()
	{
		if (IM->parserecord() == false) return false;		
		bool readstatus = true;
		bool status;
		Id.uniqueid = IM->record();
		status = IM->read(sn, Id.surveynumber); if (status == false) readstatus = false;
		status = IM->read(dn, Id.daynumber); if (status == false) readstatus = false;
		status = IM->read(fn, Id.flightnumber); if (status == false) readstatus = false;
		status = IM->read(ln, Id.linenumber); if (status == false) readstatus = false;
		status = IM->read(fidn, Id.fidnumber); if (status == false) readstatus = false;
		status = IM->read(xord, Location.x); if (status == false) readstatus = false;
		status = IM->read(yord, Location.y); if (status == false) readstatus = false;
		status = IM->read(elevation, Location.groundelevation); if (status == false) readstatus = false;

		Location.z = ud_double();
		
		status = readgeometry(fd_GI,GI); if (status == false) readstatus = false;
		status = readgeometry(fd_GR,GR); 		
		GR.fillundefined(GI);
		status = readgeometry(fd_GTFR,GTFR); 
		GTFR.fillundefined(GI); 
		status = readgeometry(fd_GS,GS); 

		status = IM->read(fd_ERc, ER.conductivity, nlayers); if (status == false) readstatus = false;
		status = IM->read(fd_ERt, ER.thickness, nlayers - 1); if (status == false) readstatus = false;

		if (solve_conductivity) {
			status = IM->read(fd_ESc, ES.conductivity, nlayers); if (status == false) readstatus = false;
		}
		if (solve_thickness) {
			status = IM->read(fd_ESt, ES.thickness, nlayers - 1); if (status == false) readstatus = false;
		}

		for (size_t si = 0; si < nsystems; si++) {
			readsystemdata(si);
		}
		return readstatus;
	}

	void readsystemdata(size_t sysindex)
	{
		cTDEmSystemInfo& S = SV[sysindex];
		S.CompInfo[0].readdata(IM);
		S.CompInfo[1].readdata(IM);
		S.CompInfo[2].readdata(IM);
	}

	bool initialise_sample()
	{
		LastIteration = 0;
		LastLambda = 1e8;

		bool status = initialise_data();
		if (status == false) return false;

		initialise_parameters();		
		initialise_Wd();
		initialise_Wc();
		initialise_Wt();
		initialise_Wg();

		if (SmoothnessMethod == eSmoothnessMethod::DERIVATIVE_1ST) {
			initialise_L_Ws_1st_derivative();
		}
		else if (SmoothnessMethod == eSmoothnessMethod::DERIVATIVE_2ND) {
			initialise_L_Ws_2nd_derivative();
		}

		initialise_Wr_Wm();

		Param = RefParam;
		LastPhiM = phiModel(Param, LastPhiC, LastPhiT, LastPhiG, LastPhiS);
		TerminationReason = "Has not terminated";

		if (OO.Dump) {
			const std::string dp = dumppath();
			GR.write(dp+"geometry_start.dat");
			ER.write(dp+"earth_start.dat");

			GR.write(dp+"geometry_ref.dat");
			ER.write(dp+"earth_ref.dat");

			GS.write(dp+"geometry_std.dat");
			ES.write(dp+"earth_std.dat");

			FILE* fp = fileopen(dp+"Id.dat","w");
			fprintf(fp, "%zu\t%zu\t%zu\t%zu\t%zu\t%lf\t%lf\t%lf\t%lf\t%lf", Id.uniqueid, Id.surveynumber, Id.daynumber, Id.flightnumber, Id.linenumber, Id.fidnumber, Location.x, Location.y, Location.groundelevation, Location.z);
			fclose(fp);
		}
		return true;
	}

	bool initialise_data()
	{
		std::vector<double> obs(ndata_all);
		std::vector<double> err(ndata_all);
		std::vector<double> pred(ndata_all);
		for (size_t si = 0; si < nsystems; si++) {
			cTDEmSystemInfo& S = SV[si];
			cTDEmSystem& T = S.T;
			if (S.reconstructPrimary) {
				T.setgeometry(GTFR);
				T.LEM.calculation_type = CT_FORWARDMODEL;
				T.LEM.derivative_layer = INT_MAX;
				T.setprimaryfields();				
				S.CompInfo[0].oP = T.PrimaryX;
				S.CompInfo[1].oP = T.PrimaryY;
				S.CompInfo[2].oP = T.PrimaryZ;
			}

			if (S.invertXPlusZ){
				for (size_t wi = 0; wi < S.nwindows; wi++) {

					//X+Z Comp
					size_t di = wi + S.xzIndex;					
					double X = S.CompInfo[0].oS[wi];
					double Z = S.CompInfo[2].oS[wi];					
					if (S.invertPrimaryPlusSecondary) {
						X += S.CompInfo[0].oP;
						Z += S.CompInfo[2].oP;
					}
					obs[di] = std::hypot(X,Z);

					const double& Xerr = S.CompInfo[0].oE[wi];
					const double& Zerr = S.CompInfo[2].oE[wi];
					err[di] = std::hypot(X*Xerr,Z*Zerr)/obs[di];
					
					//Y Comp
					if (S.CompInfo[1].Use) {
						di = S.CompInfo[1].dataindex + wi;
						obs[di] = S.CompInfo[1].oS[wi];
						if (S.invertPrimaryPlusSecondary) {
							obs[di] += S.CompInfo[1].oP;
						}
						err[di] = S.CompInfo[1].oE[wi];
					}
				}
			}
			else {
				for (size_t ci = 0; ci < 3; ci++) {
					if (S.CompInfo[ci].Use == false) continue;
					for (size_t wi = 0; wi < S.nwindows; wi++) {
						size_t di = S.CompInfo[ci].dataindex + wi;
						obs[di] = S.CompInfo[ci].oS[wi];
						if (S.invertPrimaryPlusSecondary) {
							obs[di] += S.CompInfo[ci].oP;
						}
						err[di] = S.CompInfo[ci].oE[wi];
					}
				}
			}
		}
		
		//Work out indices to be culled				
		Q.clear();
		for (size_t i = 0; i < ndata_all; i++){
			if(!isnull(obs[i])  && !isnull(err[i])) Q.push_back(i);			
		}
		ndata_inv = Q.size();

		if (ndata_inv != ndata_all) {
			size_t ncull = ndata_all - ndata_inv;
			OutputMessage += strprint(", %d null data/noise were culled",(int)ncull);
		}
		Err = cull(err);
		Obs = cull(obs);

		//Check for zero err values		
		int nzeroerr = 0;
		for (size_t i = 0; i < Err.size(); i++) {
			if (Err[i] == 0.0) nzeroerr++;
		}
		if (nzeroerr > 0) {			
			OutputMessage += strprint(", Skipped %d noise values were 0.0", nzeroerr);			
			return false;
		}
		
		if (OO.Dump) {
			write(Obs, dumppath()+"observed.dat");
			write(Err, dumppath()+"observed_std.dat");
		}
		return true;
	}

	inline bool isnull(const double& val) {
		//if (val == 29322580645161.28516) return true;
		//else if (val == 293196483870967808.00000) return true;		
		if (val >  1e14) return true;
		else if (val < -1e14) return true;
		else if (isdefined(val)) return false;		
		else return true;				
	}
	
	void initialise_parameters()
	{
		if (solve_conductivity) {
			for (size_t i = 0; i < nlayers; i++) {
				RefParam[i + cIndex] = log10(ER.conductivity[i]);
				RefParamStd[i + cIndex] = ES.conductivity[i];
			}
		}

		if (solve_thickness) {
			for (size_t i = 0; i < nlayers - 1; i++) {
				RefParam[i + tIndex] = log10(ER.thickness[i]);
				RefParamStd[i + tIndex] = ES.thickness[i];
			}
		}

		if (solve_tx_height) {
			RefParam[tx_heightIndex] = GR.tx_height;
			RefParamStd[tx_heightIndex] = GS.tx_height;
		}

		if (solve_txrx_dx) {
			RefParam[txrx_dxIndex] = GR.txrx_dx;
			RefParamStd[txrx_dxIndex] = GS.txrx_dx;
		}
		if (solve_txrx_dy) {
			RefParam[txrx_dyIndex] = GR.txrx_dy;
			RefParamStd[txrx_dyIndex] = GS.txrx_dy;
		}
		if (solve_txrx_dz) {
			RefParam[txrx_dzIndex] = GR.txrx_dz;
			RefParamStd[txrx_dzIndex] = GS.txrx_dz;
		}

		if (solve_rx_pitch) {
			RefParam[rx_pitchIndex] = GR.rx_pitch;
			RefParamStd[rx_pitchIndex] = GS.rx_pitch;
		}

		if (solve_rx_roll) {
			RefParam[rx_rollIndex] = GR.rx_roll;
			RefParamStd[rx_rollIndex] = GS.rx_roll;
		}

	}

	void initialise_Wd()
	{		
		Wd = MatrixDouble::Zero(ndata_inv, ndata_inv);		
		double s = 1.0 / (double)ndata_inv;
		for (size_t i = 0; i < ndata_inv; i++) {
			Wd(i,i) = s / (Err[i] * Err[i]);
		}		
		if (OO.Dump) writetofile(Wd, dumppath() + "Wd.dat");
	}

	void initialise_Wc()
	{
		Wc = MatrixDouble::Zero(nparam, nparam);		
		if (solve_conductivity == false)return;

		std::vector<double> t(nlayers);
		if (nlayers == 1) {
			t[0] = 1;
		}
		else if (nlayers == 2) {
			t[0] = ER.thickness[0];
			t[1] = ER.thickness[0];
		}
		else {
			for (size_t i = 0; i < (nlayers - 1); i++) {
				t[i] = ER.thickness[i];
			}
			t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3])*t[nlayers - 2];
		}


		double tsum = 0.0;
		for (size_t i = 0; i < nlayers; i++)tsum += t[i];
		double tavg = tsum / (double)nlayers;

		double s = AlphaC / (double)(nlayers);
		for (size_t i = 0; i < nlayers; i++) {
			size_t p = i + cIndex;
			Wc(p,p) = s * (t[i] / tavg) / (RefParamStd[p] * RefParamStd[p]);
		}
	}

	void initialise_Wt()
	{
		Wt = MatrixDouble::Zero(nparam, nparam);
		if (solve_thickness == false)return;

		double s = AlphaT / (double)(nlayers - 1);
		for (size_t i = 0; i < nlayers - 1; i++) {
			size_t p = i + tIndex;
			Wt(p,p) = s / (RefParamStd[p] * RefParamStd[p]);
		}

	}

	void initialise_Wg()
	{
		Wg = MatrixDouble::Zero(nparam, nparam);
		if (ngeomparam <= 0)return;

		double s = AlphaG / (double)ngeomparam;
		for (size_t i = 0; i < ngeomparam; i++) {
			size_t p = i + gIndex;
			Wg(p,p) = s / (RefParamStd[p] * RefParamStd[p]);
		}
	}

	void initialise_L_Ws_1st_derivative()
	{
		Ws = MatrixDouble::Zero(nparam, nparam);
		if (AlphaS == 0 || nlayers < 3) return;
		if (solve_conductivity == false) return;

		std::vector<double> t(nlayers);
		for (size_t i = 0; i < (nlayers - 1); i++) {
			t[i] = ER.thickness[i];
		}
		t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3])*t[nlayers - 2];


		double tsum = 0.0;
		for (size_t i = 0; i < nlayers; i++)tsum += t[i];
		double tavg = tsum / (double)nlayers;

		MatrixDouble L = MatrixDouble::Zero(nlayers - 1, nparam);
		size_t neqn = 0;
		for (size_t li = 1; li < nlayers; li++) {
			size_t pindex = cIndex + li;
			double t1 = t[li - 1];
			double t2 = t[li];
			double d12 = (t1 + t2) / 2.0;
			double s = sqrt(t2 / tavg);//sqrt because it gets squared in L'L		
			L(neqn,pindex - 1) = -s / d12;
			L(neqn,pindex) = s / d12;
			neqn++;
		}
		Ws = L.transpose() * L;
		Ws *= (AlphaS / (double)(nlayers - 1));
	}

	void initialise_L_Ws_2nd_derivative()
	{
		Ws = MatrixDouble::Zero(nparam, nparam);
		if (AlphaS == 0 || nlayers < 3) return;
		if (solve_conductivity == false) return;

		std::vector<double> t(nlayers);
		if (nlayers == 1) {
			t[0] = 1.0;
		}
		else {
			for (size_t i = 0; i < (nlayers - 1); i++) {
				t[i] = ER.thickness[i];
			}
			t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3])*t[nlayers - 2];
		}

		double tsum = 0.0;
		for (size_t i = 0; i < nlayers; i++)tsum += t[i];
		double tavg = tsum / (double)nlayers;

		MatrixDouble L = MatrixDouble::Zero(nlayers - 2, nparam);
		size_t neqn = 0;
		for (size_t li = 1; li < nlayers - 1; li++) {
			size_t pindex = cIndex + li;
			double t1 = t[li - 1];
			double t2 = t[li];
			double t3 = t[li + 1];
			double d12 = (t1 + t2) / 2.0;
			double d23 = (t2 + t3) / 2.0;
			double s = sqrt(t2 / tavg);//sqrt because it gets squared in L'L		
			L(neqn,pindex - 1) = s / d12;
			L(neqn,pindex) = -s / d12 - s / d23;
			L(neqn,pindex + 1) = s / d23;
			neqn++;
		}
		Ws = L.transpose() * L;
		Ws *= (AlphaS / (double)(nlayers - 2));
	}

	void initialise_Wr_Wm()
	{
		Wr = MatrixDouble::Zero(nparam, nparam);
		if (AlphaC > 0.0) Wr += Wc;
		if (AlphaT > 0.0) Wr += Wt;
		if (AlphaG > 0.0) Wr += Wg;

		Wm = Wr + Ws;

		if (OO.Dump) {
			writetofile(Wc, dumppath() + "Wc.dat");
			writetofile(Wt, dumppath() + "Wt.dat");
			writetofile(Wg, dumppath() + "Wg.dat");
			writetofile(Wr, dumppath() + "Wr.dat");
			writetofile(Ws, dumppath() + "Ws.dat");
			writetofile(Wm, dumppath() + "Wm.dat");
		}

	}

	VectorDouble parameterchange(const double lambda)
	{
		VectorDouble x  = solve(lambda);
		VectorDouble dm = x - Param;

		if (solve_conductivity) {
			for (size_t li = 0; li < nlayers; li++) {
				size_t pindex = li + cIndex;
				if (Param[pindex] + dm[pindex] > log10(max_conductivity)) {
					//printf("upper limit li=%zu pindex=%zu dm=%lf\n",li,pindex,dm[pindex]);
					dm[pindex] = log10(max_conductivity) - Param[pindex];
				}
				else if (Param[pindex] + dm[pindex] < log10(min_conductivity)) {
					//printf("lower limit li=%zu pindex=%zu dm=%lf\n",li,pindex,dm[pindex]);
					dm[pindex] = log10(min_conductivity) - Param[pindex];
				}
			}
		}

		if (solve_thickness) {
			for (size_t li = 0; li < nlayers - 1; li++) {
				size_t pindex = li + tIndex;
				if (dm[pindex] > 0.5) {
					//printf("li=%zu pindex=%zu dm=%lf\n",li,pindex,dm[pindex]);
					dm[pindex] = 0.5;
				}
				else if (dm[pindex] < -0.5) {
					//printf("li=%zu pindex=%zu dm=%lf\n",li,pindex,dm[pindex]);
					dm[pindex] = -0.5;
				}
			}
		}

		if (solve_tx_height) {
			size_t pindex = tx_heightIndex;
			if (dm[pindex] > 0.5) {
				//printf("li=%zu pindex=%zu dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] = 0.5;
			}
			else if (dm[pindex] < -0.5) {
				//printf("li=%zu pindex=%zu dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] = -0.5;
			}

			if (Param[pindex] + dm[pindex] > 1000) {
				dm[pindex] = 1000 - Param[pindex];
			}
			else if (Param[pindex] + dm[pindex] < 10) {
				dm[pindex] = 10 - Param[pindex];
			}
		}
		return dm;
	}

	VectorDouble solve(const double lambda)
	{
		// Phi = (d-g(m)+Jm) Wd (d-g(m)+Jm) + lambda ( (m-m0)' Wr (m-m0) + m' Ws m) )
		//Ax = b
		//A = [J'WdJ + lambda (Wr + Ws)]
		//x = m(n+1)
		//b = J'Wd(d - g(m) + Jm) + lambda*Wr*m0
		//dm = m(n+1) - m = x - m

		const VectorDouble& m = Param;
		const VectorDouble& d = Obs;
		const VectorDouble& g = Pred;
		const VectorDouble& e = Err;
		const VectorDouble& m0 = RefParam;


		MatrixDouble V = Wd;
		if (NormType == eNormType::L2) {

		}
		else {
			for (size_t i = 0; i < ndata_inv; i++) {
				const double r = (d[i] - g[i]) / e[i];
				V(i,i) *= 1.0 / std::abs(r);
			}
		}

		MatrixDouble JtV = J.transpose() * V;
		MatrixDouble JtVJ = JtV * J;

		VectorDouble b = JtV * (d - g + J * m) + lambda * (Wr * m0);
		MatrixDouble A = JtVJ + lambda * Wm;
		VectorDouble x = pseudoInverse(A)*b;
		return x;
	}

	double l1_norm(const VectorDouble& g)
	{
		double l1 = 0.0;
		for (size_t i = 0; i < ndata_inv; i++) {
			l1 += std::abs(Obs[i] - g[i]) / Err[i];
		}
		return l1 / ndata_inv;
	}

	double l2_norm(const VectorDouble& g)
	{
		VectorDouble v = Obs - g;
		VectorDouble a = Wd * v;
		double l2 = mtDm(v, Wd);
		return l2;
	}

	double phiData(const VectorDouble& g)
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
			glog.warningmsg(_SRC_, "Caught invalid PhiD\n");
		}
		return phid;
	}

	double phiModel(const VectorDouble& p)
	{
		double phic, phit, phig, phis;
		return phiModel(p, phic, phit, phig, phis);
	}

	double phiModel(const VectorDouble& p, double& phic, double& phit, double& phig, double& phis)
	{
		phic = phiC(p);
		phit = phiT(p);
		phig = phiG(p);
		phis = phiS(p);

		double v = phic + phit + phig + phis;
		return v;
	}

	double phiC(const VectorDouble& p)
	{
		if (AlphaC == 0.0)return 0.0;
		if (solve_conductivity == false)return 0.0;
		VectorDouble v = p - RefParam;
		return mtDm(v, Wc);
	}

	double phiT(const VectorDouble& p)
	{
		if (AlphaT == 0.0)return 0.0;
		if (solve_thickness == false)return 0.0;
		VectorDouble v = p - RefParam;
		return mtDm(v, Wt);
	}

	double phiG(const VectorDouble& p)
	{
		if (AlphaG == 0.0)return 0.0;
		if (ngeomparam == 0)return 0.0;
		VectorDouble v = p - RefParam;
		return mtDm(v, Wg);
	}

	double phiS(const VectorDouble& p)
	{
		if (AlphaS == 0)return 0.0;
		else return mtAm(p, Ws);
	}

	cEarth1D get_earth(const VectorDouble& parameters)
	{
		cEarth1D e = ER;
		if (solve_conductivity) {
			for (size_t li = 0; li < nlayers; li++) {
				e.conductivity[li] = pow10(parameters[li + cIndex]);
			}
		}

		if (solve_thickness) {
			for (size_t li = 0; li < nlayers - 1; li++) {
				e.thickness[li] = pow10(parameters[li + tIndex]);
			}
		}
		return e;
	}

	cTDEmGeometry get_geometry(const VectorDouble& parameters)
	{
		cTDEmGeometry g = GI;
		if (solve_tx_height)g.tx_height = parameters[tx_heightIndex];
		if (solve_txrx_dx)g.txrx_dx = parameters[txrx_dxIndex];
		if (solve_txrx_dy)g.txrx_dy = parameters[txrx_dyIndex];
		if (solve_txrx_dz)g.txrx_dz = parameters[txrx_dzIndex];
		if (solve_rx_pitch)g.rx_pitch = parameters[rx_pitchIndex];
		if (solve_rx_roll)g.rx_roll = parameters[rx_rollIndex];
		return g;
	}

	void set_predicted()
	{
		for (size_t si = 0; si < nsystems; si++) {
			cTDEmSystemInfo& S = SV[si];
			cTDEmSystem& T = S.T;

			cTDEmData& d = S.predicted;
			d.xcomponent().Primary = T.PrimaryX;
			d.ycomponent().Primary = T.PrimaryY;
			d.zcomponent().Primary = T.PrimaryZ;
			d.xcomponent().Secondary = T.X;
			d.ycomponent().Secondary = T.Y;
			d.zcomponent().Secondary = T.Z;
		}
	}

	VectorDouble cull(const VectorDouble& vall) const {
		assert(Q.size() == ndata_inv);
		assert(vall.size() == ndata_all);
		VectorDouble vinv(ndata_inv);
		for (size_t i = 0; i < ndata_inv; i++) {
			vinv[i] = vall[Q[i]];
		}
		return vinv;
	}

	VectorDouble cull(const std::vector<double>& vall) const {
		assert(Q.size() == ndata_inv);
		assert(vall.size() == ndata_all);
		VectorDouble vinv(ndata_inv);
		for (size_t i = 0; i < ndata_inv; i++) {
			vinv[i] = vall[Q[i]];
		}
		return vinv;
	}

	MatrixDouble cull(const MatrixDouble& mall) const {
		assert(Q.size() == ndata_inv);
		assert(mall.rows() == ndata_all);
		assert(mall.cols() == nparam);		
		MatrixDouble minv(ndata_inv, nparam);
		for (size_t i = 0; i < ndata_inv; i++) {
			minv.row(i) = mall.row(Q[i]);
		}
		return minv;
	}


	void forwardmodel(const VectorDouble& parameters, VectorDouble& predicted, bool computederivatives)
	{
		VectorDouble pred(ndata_all);
		MatrixDouble M;
		if(computederivatives)M.resize(ndata_all, nparam);

		cEarth1D      e = get_earth(parameters);
		cTDEmGeometry g = get_geometry(parameters);
		for (size_t si = 0; si < nsystems; si++) {
			cTDEmSystemInfo& S = SV[si];
			cTDEmSystem& T = S.T;
			const size_t nw = T.NumberOfWindows;
			T.setconductivitythickness(e.conductivity, e.thickness);
			T.setgeometry(g);

			//Forwardmodel
			T.LEM.calculation_type = CT_FORWARDMODEL;
			T.LEM.derivative_layer = INT_MAX;
			T.setupcomputations();
			T.setprimaryfields();
			T.setsecondaryfields();

			std::vector<double> xfm = T.X;
			std::vector<double> yfm = T.Y;
			std::vector<double> zfm = T.Z;
			std::vector<double> xzfm;
			if (S.invertPrimaryPlusSecondary) {
				xfm += T.PrimaryX;
				yfm += T.PrimaryY;
				zfm += T.PrimaryZ;
			}

			if (S.invertXPlusZ) {
				xzfm.resize(T.NumberOfWindows);
				for (size_t wi = 0; wi < T.NumberOfWindows; wi++) {
					xzfm[wi] = std::hypot(xfm[wi], zfm[wi]);
				}
			}

			if (S.invertXPlusZ) {
				for (size_t wi = 0; wi < nw; wi++) {
					pred[wi + S.xzIndex] = xzfm[wi];
					if (S.CompInfo[1].Use) pred[wi + S.CompInfo[1].dataindex] = yfm[wi];
				}
			}
			else {
				for (size_t wi = 0; wi < nw; wi++) {
					if (S.CompInfo[0].Use) pred[wi + S.CompInfo[0].dataindex] = xfm[wi];
					if (S.CompInfo[1].Use) pred[wi + S.CompInfo[1].dataindex] = yfm[wi];
					if (S.CompInfo[2].Use) pred[wi + S.CompInfo[2].dataindex] = zfm[wi];
				}				
			}
						
			if (computederivatives) {				
				std::vector<double> xdrv(nw);
				std::vector<double> ydrv(nw);
				std::vector<double> zdrv(nw);

				if (solve_conductivity) {
					for (size_t li = 0; li < nlayers; li++) {
						size_t pindex = li + cIndex;
						T.LEM.calculation_type = CT_CONDUCTIVITYDERIVATIVE;
						T.LEM.derivative_layer = li;
						T.setprimaryfields();
						T.setsecondaryfields();

						fillDerivativeVectors(S, xdrv, ydrv, zdrv);
						//multiply by natural log(10) as parameters are in logbase10 units
						double sf = log(10.0)*e.conductivity[li];
						xdrv *= sf; ydrv *= sf; zdrv *= sf;
						fillMatrixColumn(M, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
					}
				}

				if (solve_thickness) {
					for (size_t li = 0; li < nlayers - 1; li++) {
						size_t pindex = li + tIndex;
						T.LEM.calculation_type = CT_THICKNESSDERIVATIVE;
						T.LEM.derivative_layer = li;
						T.setprimaryfields();
						T.setsecondaryfields();
						fillDerivativeVectors(S, xdrv, ydrv, zdrv);
						//multiply by natural log(10) as parameters are in logbase10 units
						double sf = log(10.0)*e.thickness[li];
						xdrv *= sf; ydrv *= sf; zdrv *= sf;
						fillMatrixColumn(M, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
					}
				}

				if (solve_tx_height) {
					size_t pindex = tx_heightIndex;
					T.LEM.calculation_type = CT_HDERIVATIVE;
					T.LEM.derivative_layer = INT_MAX;
					T.setprimaryfields();
					T.setsecondaryfields();
					fillDerivativeVectors(S, xdrv, ydrv, zdrv);
					fillMatrixColumn(M, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}

				if (solve_txrx_dx) {
					size_t pindex = txrx_dxIndex;
					T.LEM.calculation_type = CT_XDERIVATIVE;
					T.LEM.derivative_layer = INT_MAX;
					T.setprimaryfields();
					T.setsecondaryfields();
					fillDerivativeVectors(S, xdrv, ydrv, zdrv);
					fillMatrixColumn(M, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}

				if (solve_txrx_dy) {
					size_t pindex = txrx_dyIndex;
					T.LEM.calculation_type = CT_YDERIVATIVE;
					T.LEM.derivative_layer = INT_MAX;
					T.setprimaryfields();
					T.setsecondaryfields();
					fillDerivativeVectors(S, xdrv, ydrv, zdrv);
					fillMatrixColumn(M, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}

				if (solve_txrx_dz) {
					size_t pindex = txrx_dzIndex;
					T.LEM.calculation_type = CT_ZDERIVATIVE;
					T.LEM.derivative_layer = INT_MAX;
					T.setprimaryfields();
					T.setsecondaryfields();
					fillDerivativeVectors(S, xdrv, ydrv, zdrv);
					fillMatrixColumn(M, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}

				if (solve_rx_pitch) {
					size_t pindex = rx_pitchIndex;
					T.drx_pitch(xfm, zfm, g.rx_pitch, xdrv, zdrv);
					ydrv *= 0.0;
					fillMatrixColumn(M, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}

				if (solve_rx_roll) {
					size_t pindex = rx_rollIndex;
					T.drx_roll(yfm, zfm, g.rx_roll, ydrv, zdrv);
					xdrv *= 0.0;
					fillMatrixColumn(M, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}								
			}
		}
				
		predicted = cull(pred);		
		if (computederivatives)	J = cull(M);		
	}

	void fillDerivativeVectors(cTDEmSystemInfo& S, std::vector<double>& xdrv, std::vector<double>& ydrv, std::vector<double>& zdrv)
	{
		cTDEmSystem& T = S.T;
		xdrv = T.X;
		ydrv = T.Y;
		zdrv = T.Z;
		if (S.invertPrimaryPlusSecondary) {
			xdrv += T.PrimaryX;
			ydrv += T.PrimaryY;
			zdrv += T.PrimaryZ;
		}
	}

	void fillMatrixColumn(MatrixDouble& M, cTDEmSystemInfo& S, const size_t& pindex, const std::vector<double>& xfm, const std::vector<double>& yfm, const std::vector<double>& zfm, const std::vector<double>& xzfm, const std::vector<double>& xdrv, const std::vector<double>& ydrv, const std::vector<double>& zdrv)
	{
		const size_t nw = S.T.NumberOfWindows;
		if (S.invertXPlusZ) {
			for (size_t w = 0; w < nw; w++) {				
				M(w + S.xzIndex, pindex) = (xfm[w] * xdrv[w] + zfm[w] * zdrv[w]) / xzfm[w];
				if (S.CompInfo[1].Use)M(w + S.CompInfo[1].dataindex,pindex) = xdrv[w];
			}
		}
		else {
			for (size_t w = 0; w < nw; w++) {				
				if (S.CompInfo[0].Use)M(w + S.CompInfo[0].dataindex,pindex) = xdrv[w];
				if (S.CompInfo[1].Use)M(w + S.CompInfo[1].dataindex,pindex) = ydrv[w];
				if (S.CompInfo[2].Use)M(w + S.CompInfo[2].dataindex,pindex) = zdrv[w];				
			}
		}
	}

	VectorDouble compute_parameter_sensitivity()
	{
		VectorDouble s = VectorDouble::Zero(nparam);
		for (size_t pi = 0; pi < nparam; pi++) {
			for (size_t di = 0; di < ndata_inv; di++) {
				s[pi] += (std::fabs(J(di,pi)) * std::sqrt((double)ndata_inv*Wd(di,di)));
			}
		}		
		return s;
	}

	VectorDouble compute_parameter_uncertainty()
	{
		MatrixDouble JtWdJ = J.transpose()*Wd*J;
		MatrixDouble iCm = MatrixDouble::Zero(nparam, nparam);
		for (size_t i = 0; i < nparam; i++) {
			iCm(i,i) = 1.0 / (RefParamStd[i] * RefParamStd[i]);
		}
		MatrixDouble X = (double)ndata_inv*JtWdJ + iCm;		
		MatrixDouble pinvX = pseudoInverse(X);
		VectorDouble s(nparam);
		for (size_t i = 0; i < nparam; i++) {
			s[i] = std::sqrt(pinvX(i,i));
		}
		return s;
	}

	bool invert()
	{
		OutputMessage = "";
		if (parserecord() == false) {
			OutputMessage += ", Skipping - could not parse record";
			return false;
		}

		if (initialise_sample() == false) {
			return false;
		}

		iterate(); 		
		return true;
	}

	void iterate()
	{
		double percentchange = 100.0;
		bool   keepiterating = true;

		while (keepiterating == true) {

			forwardmodel(Param, Pred, true);
			LastPhiD = phiData(Pred);
			TargetPhiD = std::max(LastPhiD * 0.7, MinimumPhiD);
			if (OO.Dump) {
				makedirectorydeep(dumppath());
				writetofile(Obs,dumppath() + "d.dat");
				writetofile(Err,dumppath() + "e.dat");
				writetofile(Pred,dumppath() + "g.dat");
				EM = get_earth(Param);
				GM = get_geometry(Param);
				EM.write(dumppath()+"earth_inv.dat");
				GM.write(dumppath()+"geometry_inv.dat");
				save_iteration_file();
			}

			if (LastIteration >= MaxIterations) {
				keepiterating = false;
				TerminationReason = "Too many iterations";
			}
			else if (LastPhiD <= MinimumPhiD) {
				keepiterating = false;
				TerminationReason = "Reached minimum";
			}
			else if (percentchange < MinimumImprovement) {
				keepiterating = false;
				TerminationReason = "Small % improvement";
			}
			else {
				cTrial t = targetsearch(LastLambda, TargetPhiD);				
				VectorDouble dm = parameterchange(t.lambda);
				VectorDouble mtemp = Param + (t.stepfactor*dm);
				VectorDouble gtemp(ndata_inv);
				forwardmodel(mtemp, gtemp, false);
				double phidtemp = phiData(gtemp);
				percentchange = 100.0*(LastPhiD - phidtemp) / (LastPhiD);

				if (phidtemp < LastPhiD) {
					Param = mtemp;
					Pred = gtemp;
					LastPhiD = phidtemp;
					LastLambda = t.lambda;
					LastPhiM = phiModel(Param, LastPhiC, LastPhiT, LastPhiG, LastPhiS);
				}
			}
			LastIteration++;
		}

		EM = get_earth(Param);
		GM = get_geometry(Param);
		forwardmodel(Param, Pred, false);
		set_predicted();

		forwardmodel(Param, Pred, true);
		ParameterSensitivity = compute_parameter_sensitivity();
		ParameterUncertainty = compute_parameter_uncertainty();
	}

	void save_iteration_file() {
		FILE* fp = fileopen(dumppath() + "iteration.dat", "w");
		fprintf(fp, "Iteration\t%zu\n", LastIteration);
		fprintf(fp, "TargetPhiD\t%lf\n", TargetPhiD);
		fprintf(fp, "PhiD\t%lf\n", LastPhiD);
		fprintf(fp, "Lambda\t%lf\n", LastLambda);
		fprintf(fp, "\n");
		fclose(fp);
	};

	cTrial targetsearch(const double currentlambda, const double target)
	{
		cTrialCache T;
		T.target = target;
		eBracketResult b = brackettarget(T, target, currentlambda);
		cTrial t{};
		if (b == eBracketResult::BRACKETED) {
			//bracketed target - find with Brents Method
			double newphid = DBL_MIN;
			double lambda = brentsmethod(T, target, newphid);
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
			glog.errormsg(_SRC_, "targetsearch(): Unknown value %d returned from target brackettarget()\n", b);
		}
		//T.print(t.lambda, t.phid);
		return t;		
	}

	bool istargetbraketed(cTrialCache& T)
	{
		double target = T.target;
		T.sort_lambda();
		for (size_t i = T.trial.size() - 1; i >= 1; i--) {
			if (T.trial[i].phid >= target && T.trial[i - 1].phid <= target) {
				return true;
			}
			if (T.trial[i].phid <= target && T.trial[i - 1].phid >= target) {
				return true;
			}
		}
		return false;
	}

	bool isminbraketed(cTrialCache& T)
	{
		size_t index = T.minphidindex();
		if (index == 0 || index == T.trial.size() - 1) {
			return false;
		}

		double fa = T.trial[index - 1].phid;
		double fb = T.trial[index].phid;
		double fc = T.trial[index + 1].phid;
		if ((fb < fa) && (fb < fc)) {
			return true;
		}
		return false;
	}

	eBracketResult brackettarget(cTrialCache& T, const double target, const double currentlambda)
	{
		double startx = log10(currentlambda);
		if (LastIteration == 0) {
			std::vector<double> x;
			x.push_back(8); x.push_back(6);
			x.push_back(4); x.push_back(2);
			x.push_back(1); x.push_back(0);
			for (size_t k = 0; k < x.size(); k++) {
				trialfunction(T, pow10(x[k]));
				bool tarbrak = istargetbraketed(T);
				if (tarbrak) {
					return eBracketResult::BRACKETED;//target bracketed		
				}
			}

			double minv = DBL_MAX;
			for (size_t k = 0; k < T.trial.size(); k++) {
				if (fabs(T.trial[k].phid - target) < minv) {
					minv = fabs(T.trial[k].phid - target);
					startx = log10(T.trial[k].lambda);
				}
			}
		}
		else {
			trialfunction(T, pow10(startx));
		}
		//T.print(LastLambda,LastPhiD);
		//prompttocontinue();

		std::vector<double> x;
		x.push_back(+1); x.push_back(-1);
		x.push_back(+2); x.push_back(-2);
		x.push_back(+3); x.push_back(-3);
		for (size_t k = 0; k < x.size(); k++) {
			//T.print(LastLambda, LastPhiD);
			trialfunction(T, pow10(startx + x[k]));
			bool tarbrak = istargetbraketed(T);
			if (tarbrak) {
				return eBracketResult::BRACKETED;//target bracketed		
			}
		}		
		//T.print(LastLambda, LastPhiD);
		//prompttocontinue();

		if (T.maxphid() < target) {
			return eBracketResult::ALLBELOW;//all below target	
		}

		bool minbrak = isminbraketed(T);
		if (minbrak)return eBracketResult::MINBRACKETED;//min bracketed											

		return eBracketResult::ALLABOVE;//all above target
	}

	double brentsmethod(cTrialCache& T, const double target, double& newphid)
	{
		T.sort_lambda();
		size_t index = T.trial.size() - 1;
		for (size_t i = T.trial.size() - 1; i >= 1; i--) {
			double f1 = T.trial[i].phid - target;
			double f2 = T.trial[i - 1].phid - target;
			if (f1 * f2 <= 0.0) {
				index = i;
				break;
			}
		}

		//Adapted from http://en.wikipedia.org/wiki/Brent's_method
		double xerrorTol = 0.01;
		double yerrorTol = target * 0.1;//10% accuracy is good enough

		double a = log10(T.trial[index - 1].lambda);
		double b = log10(T.trial[index].lambda);
		double fa = T.trial[index - 1].phid - target;
		double fb = T.trial[index].phid - target;
		if (fa * fb >= 0.0) {
			glog.warningmsg(_SRC_, "Target must be bracketed\n");
		}

		double c = 0;
		double d = DBL_MAX;
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
			fs = trialfunction(T, pow10(s)) - target;
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
				//T.print(LastLambda, LastPhiD);				
				glog.warningmsg(_SRC_, "Too many bisections\n");
				newphid = fb + target;
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

	double goldensearch(double a, double b, double c, double xtol, const double lambda, const VectorDouble& m, const VectorDouble& dm, VectorDouble& g, cTrialCache& cache)
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
			VectorDouble p = m + x * dm;
			forwardmodel(p, g, false);
			fx = phiData(g);
			t.stepfactor = x;
			t.phid = fx;
			t.phim = phiModel(p);
			t.order = cache.trial.size();
			t.lambda = lambda;
			cache.trial.push_back(t);
			//printf("%lf %lf\n",x,fx);
		}

		double fb = cache.sfsearch(b);
		if (fb < 0) {
			cTrial t;
			VectorDouble p = m + b * dm;
			forwardmodel(p, g, false);
			fb = phiData(g);
			t.stepfactor = b;
			t.phid = fb;
			t.phim = phiModel(p);
			t.order = cache.trial.size();
			t.lambda = lambda;
			cache.trial.push_back(t);
			//printf("%lf %lf\n",b,fb);
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

	double trialfunction(cTrialCache& T, const double triallambda)
	{
		VectorDouble dm(nparam);
		VectorDouble p(nparam);
		VectorDouble g(ndata_inv);
		dm = parameterchange(triallambda);
		cTrialCache cache;
		cTrial t0;
		t0.phid = LastPhiD;
		t0.phim = LastPhiM;
		t0.stepfactor = 0.0;
		t0.lambda = triallambda;
		t0.order = cache.trial.size();
		cache.trial.push_back(t0);

		cTrial t1;
		p = Param + dm;
		forwardmodel(p, g, false);
		t1.phid = phiData(g);
		t1.phim = phiModel(p);
		t1.stepfactor = 1.0;
		t1.lambda = triallambda;
		t1.order = cache.trial.size();
		cache.trial.push_back(t1);

		double pcdiff = 100 * (t1.phid - t0.phid) / t0.phid;
		if (pcdiff > 0.0 || pcdiff < -1.0) {
			//ie dont do not do golden search
			//if only tiny improvement				
			double xtol = 0.1;
			double gsf = goldensearch(0.0, 0.38196601125010510, 1.0, xtol, triallambda, Param, dm, g, cache);
			cTrial t3;
			p = Param + gsf * dm;
			forwardmodel(p, g, false);
			t3.phid = phiData(g);
			t3.phim = phiModel(p);
			t3.stepfactor = gsf;
			t3.lambda = triallambda;
			t3.order = cache.trial.size();
			cache.trial.push_back(t3);
			//double gsf = goldensearch(0.0,0.5,1.0,tau,Param,dm,g,cache);
			//printf("gsf=%lf\n",gsf);				
		}
		//printtrials(cache);
		//prompttocontinue();
		size_t minindex = cache.minphidindex();

		cTrial t = cache.trial[minindex];
		t.order = T.trial.size();
		T.trial.push_back(t);
		return t.phid;
	}

	void writeresult(const int& pointindex)
	{		
		const int& pi = pointindex;
		OM->begin_point_output();
		
		//Id		
		OM->writefield(pi, Id.uniqueid, "uniqueid", "Inversion sequence number", UNITLESS, 1, NC_UINT, DN_NONE, 'I', 12, 0);
		OM->writefield(pi, Id.surveynumber, "survey", "Survey number", UNITLESS, 1, NC_UINT, DN_NONE, 'I', 12, 0);
		OM->writefield(pi, Id.daynumber, "date", "Date number", UNITLESS, 1, NC_UINT, DN_NONE, 'I', 12, 0);
		//TODO fix the example file
		//OM->writefield(pi, Id.flightnumber, "flight", "Flight number", UNITLESS, 1, NC_UINT, DN_NONE, 'I', 12, 0);
		OM->writefield(pi, Id.linenumber, "line", "Line number", UNITLESS, 1, NC_UINT, DN_NONE, 'I', 12, 0);
		OM->writefield(pi, Id.fidnumber, "fiducial", "Fiducial number", UNITLESS, 1, NC_DOUBLE, DN_NONE, 'F', 12, 2);

		//Location
		OM->writefield(pi, Location.x, "easting", "UTM Easting", "m", 1, NC_DOUBLE, DN_NONE, 'F', 10, 1);
		OM->writefield(pi, Location.y, "northing", "UTM Northing", "m", 1, NC_DOUBLE, DN_NONE, 'F', 10, 1);			
		OM->writefield(pi, Location.groundelevation, "elevation", "Ground elevation relative to sea-level", "m", 1, NC_FLOAT, DN_NONE, 'F', 10, 2);
		
		//Geometry Input
		bool invertedfieldsonly = false;
		for (size_t i = 0; i < GI.size(); i++) {
			if (invertedfieldsonly && solvegeometryindex(i) == false)continue;
			OM->writefield(pi, GI[i], "input_" + GI.fname(i), "Input " + GI.description(i), GI.units(i), 1, NC_FLOAT, DN_NONE, 'F', 9, 2);
		}

		//Geometry Modelled		
		invertedfieldsonly = true;
		for (size_t i = 0; i < GM.size(); i++) {
			if (invertedfieldsonly && solvegeometryindex(i) == false)continue;
			OM->writefield(pi, GM[i], "inverted_" + GM.fname(i), "Inverted " + GI.description(i), GI.units(i), 1, NC_FLOAT, DN_NONE, 'F', 9, 2);
		}
				
		//ndata
		OM->writefield(pi,
			ndata_inv, "ndata", "Number of data in inversion", UNITLESS,
			1, NC_UINT, DN_NONE, 'I', 4, 0);

		//Earth	
		OM->writefield(pi,
			nlayers,"nlayers","Number of layers ", UNITLESS,
			1, NC_UINT, DN_NONE, 'I', 4, 0);
		
		OM->writefield(pi,
			EM.conductivity, "conductivity", "Layer conductivity", "S/m",
			EM.conductivity.size(), NC_FLOAT, DN_LAYER, 'E', 15, 6);
		
		double bottomlayerthickness = 100.0;
		if (solve_thickness == false && nlayers > 1) {
			bottomlayerthickness = EM.thickness[nlayers - 2];
		}
		std::vector<double> thickness = EM.thickness;
		thickness.push_back(bottomlayerthickness);

		OM->writefield(pi,
			thickness, "thickness", "Layer thickness", "m",
			thickness.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
					
				
		if (OO.PositiveLayerTopDepths) {			
			std::vector<double> dtop = EM.layer_top_depth();
			OM->writefield(pi,
				dtop, "depth_top", "Depth to top of layer", "m",
				dtop.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
		}

		if (OO.NegativeLayerTopDepths) {
			std::vector<double> ndtop = -1.0*EM.layer_top_depth();
			OM->writefield(pi,
				ndtop, "depth_top_negative", "Negative of depth to top of layer", "m",
				ndtop.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
		}
		
		if (OO.PositiveLayerBottomDepths) {
			std::vector<double> dbot = EM.layer_bottom_depth();
			OM->writefield(pi,
				dbot, "depth_bottom", "Depth to bottom of layer", "m",
				dbot.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
		}

		if (OO.NegativeLayerBottomDepths) {
			std::vector<double> ndbot = -1.0 * EM.layer_bottom_depth();
			OM->writefield(pi,
				ndbot, "depth_bottom_negative", "Negative of depth to bottom of layer", "m",
				ndbot.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
		}

		if (OO.InterfaceElevations) {			
			std::vector<double> etop = EM.layer_top_depth();
			etop += Location.groundelevation;
			OM->writefield(pi,
				etop, "elevation_interface", "Elevation of interface", "m",
				etop.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
		}
				
		if (OO.ParameterSensitivity) {
			std::vector<double> ps = copy(ParameterSensitivity);
			if (solve_conductivity) {
				std::vector<double> v(ps.begin() + cIndex, ps.begin() + cIndex + nlayers);
				OM->writefield(pi,
					v, "conductivity_sensitivity", "Conductivity parameter sensitivity", UNITLESS,
					v.size(), NC_FLOAT, DN_LAYER, 'E', 15, 6);
			}
			
			if (solve_thickness) {
				std::vector<double> v(ps.begin() + tIndex, ps.begin() + tIndex + nlayers-1);
				v.push_back(0.0);//halfspace layer not a parameter
				OM->writefield(pi,
					v, "thickness_sensitivity", "Thickness parameter sensitivity", UNITLESS,
					v.size(), NC_FLOAT, DN_LAYER, 'E', 15, 6);
			}

			size_t k = 0;
			for (size_t gi = 0; gi < GI.size(); gi++) {				
				if (solvegeometryindex(gi) == false)continue;
				std::string name = "inverted_" + GI.fname(gi) + "_sensitivity";
				std::string desc = GI.description(gi) + " parameter sensitivity";
				OM->writefield(pi, 
					ps[gIndex+k], name, desc, UNITLESS,
					1, NC_FLOAT, DN_NONE, 'E', 15, 6);
				k++;
			}
		}

		if (OO.ParameterUncertainty) {
			std::vector<double> pu = copy(ParameterUncertainty);
			if (solve_conductivity) {
				std::vector<double> v(pu.begin() + cIndex, pu.begin() + cIndex + nlayers);
				OM->writefield(pi,
					v, "conductivity_uncertainty", "Conductivity parameter uncertainty", "log10(S/m)",
					v.size(), NC_FLOAT, DN_LAYER, 'E', 15, 6);
			}

			if (solve_thickness) {
				std::vector<double> v(pu.begin() + tIndex, pu.begin() + tIndex + nlayers - 1);
				v.push_back(0.0);//halfspace layer not a parameter
				OM->writefield(pi,
					v, "thickness_uncertainty", "Thickness parameter uncertainty", "log10(m)",
					v.size(), NC_FLOAT, DN_LAYER, 'E', 15, 6);
			}

			size_t k = 0;
			for (size_t gi = 0; gi < GI.size(); gi++) {
				if (solvegeometryindex(gi) == false)continue;
				std::string name = "inverted_" + GI.fname(gi) + "_uncertainty";
				std::string desc = GI.description(gi) + " parameter uncertainty";
				OM->writefield(pi,
					pu[gIndex + k], name, desc, GI.units(gi),
					1, NC_FLOAT, DN_NONE, 'E', 15, 6);
				k++;
			}
		}

				
		//ObservedData
		if (OO.ObservedData) {			
			for (size_t si = 0; si < nsystems; si++) {
				cTDEmSystemInfo& S = SV[si];
				for (size_t ci = 0; ci < 3; ci++) {
					if (S.CompInfo[ci].Use) writeresult_emdata(pi,
						si, S.CompInfo[ci].Name,
						"observed", "Observed",
						'E', 15, 6, S.CompInfo[ci].oP, S.CompInfo[ci].oS, S.invertPrimaryPlusSecondary);
				}
			}
		}

		
		//Noise Estimates
		if (OO.NoiseEstimates) {
			for (size_t si = 0; si < nsystems; si++) {
				cTDEmSystemInfo& S = SV[si];
				for (size_t ci = 0; ci < 3; ci++) {
					if (S.CompInfo[ci].Use) writeresult_emdata(pi,
						si, S.CompInfo[ci].Name,
						"noise", "Estimated noise",						
						'E', 15, 6, 0.0, S.CompInfo[ci].oE, false);
				}
			}
		}
		
		//PredictedData
		if (OO.PredictedData) {
			for (size_t si = 0; si < nsystems; si++) {
				cTDEmSystemInfo& S = SV[si];
				for (size_t ci = 0; ci < 3; ci++) {
					if (S.CompInfo[ci].Use) writeresult_emdata(pi,
						si, S.CompInfo[ci].Name, "predicted", "Predicted", 'E', 15, 6,
						S.predicted.component(ci).Primary,
						S.predicted.component(ci).Secondary,
						S.invertPrimaryPlusSecondary);
				}
			}
		}
		

		//Inversion parameters and norms
		OM->writefield(pi, AlphaC, "AlphaC", "AlphaC inversion parameter", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, AlphaT, "AlphaT", "AlphaT inversion parameter", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, AlphaG, "AlphaG", "AlphaG inversion parameter", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, AlphaS, "AlphaS", "AlphaS inversion parameter", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);		
		OM->writefield(pi, LastPhiD, "PhiD", "Normalised data misfit", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, LastPhiM, "PhiM", "Combined model norm", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, LastPhiC, "PhiC", "Conductivity model norm", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, LastPhiT, "PhiT", "Thickness model norm", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, LastPhiG, "PhiG", "Geometry model norm", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, LastPhiS, "PhiS", "Smoothness model norm", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, LastLambda, "Lambda", "Lambda regularization parameter", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);		
		OM->writefield(pi, LastIteration, "Iterations", "Number of iterations", UNITLESS, 1, NC_UINT, DN_NONE, 'I', 4, 0);
				
		//End of record book keeping
		OM->end_point_output();		
		if (pointsoutput == 0) {			
			OM->end_first_record();//only do this once		
		}
		pointsoutput++;
	};

	void writeresult_emdata(const int& pointindex, const size_t& sysnum, const std::string& comp, const std::string& nameprefix, const std::string& descprefix, const char& form, const int& width, const int& decimals, const double& p, std::vector<double>& s, const bool& includeprimary)
	{
		std::string DN_WINDOW = "em_window";
		std::string sysname = nameprefix + strprint("_EMSystem_%d_", (int)sysnum + 1);
		std::string sysdesc = descprefix + strprint(" EMSystem %d ", (int)sysnum + 1);
		if (includeprimary) {
			std::string name = sysname + comp + "P";
			std::string desc = sysdesc + comp + "-component primary field";			
			OM->writefield(pointindex,
				p, name, desc, UNITLESS,
				1, NC_FLOAT, DN_NONE, form, width, decimals);			
		}

		{
			std::string name = sysname + comp + "S";
			std::string desc = sysdesc + comp + "-component secondary field";
			OM->writefield(pointindex,
				s, name, desc, UNITLESS,
				s.size(), NC_FLOAT, DN_WINDOW, form, width, decimals);
		}
	}

	bool solvegeometryindex(const size_t index)
	{
		eGeometryElementType getype = cTDEmGeometry::elementtype(index);

		switch (getype) {
		case GE_TX_HEIGHT: return solve_tx_height; break;
		case GE_TX_ROLL:   return solve_tx_roll; break;
		case GE_TX_PITCH:  return solve_tx_pitch; break;
		case GE_TX_YAW:    return solve_tx_yaw; break;
		case GE_TXRX_DX:   return solve_txrx_dx; break;
		case GE_TXRX_DY:   return solve_txrx_dy; break;
		case GE_TXRX_DZ:   return solve_txrx_dz; break;
		case GE_RX_ROLL:   return solve_rx_roll; break;
		case GE_RX_PITCH:  return solve_rx_pitch; break;
		case GE_RX_YAW:    return solve_rx_yaw; break;
		default:
			glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
			return false;
		}
	}

	bool readgeometry(const std::vector<cFieldDefinition>& gfd, cTDEmGeometry& g)
	{		
		bool status = true;
		for (size_t i = 0; i < g.size(); i++) {
			bool istatus = IM->read(gfd[i], g[i]);
			if (istatus == false) {
				status = false;
			}
		}
		return status;
	}
	
	void write(const VectorDouble& v, std::string path) const 
	{
		FILE* fp = fileopen(path, "w");
		for (auto  i = 0; i < v.rows(); i++) {
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


	std::string dumppath() const
	{
		std::string s = OO.DumpPath(IM->record(),LastIteration);
		return s;
	};
	

};

#endif
