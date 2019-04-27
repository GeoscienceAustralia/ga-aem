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

#include "blocklanguage.h"
#include "fielddefinition.h"
#include "matrix_ops.h"
#include "airborne_types.h"
#include "tdemsystem.h"
#include "geophysics_netcdf.h"

#define VERSION "1.0"

#if defined _OPENMP
	#include <omp.h>
	//This external thread lock must be set when fftw is being initialised
	extern omp_lock_t fftw_thread_lock;
#endif

enum eNormType { L1, L2 };
enum eSmoothnessMethod { SM_1ST_DERIVATIVE, SM_2ND_DERIVATIVE };
enum eBracketResult { BR_BRACKETED, BR_MINBRACKETED, BR_ALLABOVE, BR_ALLBELOW };

class cTrial{

public:
	size_t order;
	double lambda;
	double phid;
	double phim;
	double stepfactor;

	static bool lambda_compare(const cTrial& a, cTrial& b)
	{				
		if (a.lambda < b.lambda) return true;
		else if (a.lambda > b.lambda) return false;
		else{
			if (a.stepfactor <= b.stepfactor) return true;				
			else return false;
		}		
	}

	static bool phid_compare(const cTrial& a, cTrial& b)
	{
		if (a.phid < b.phid) return true;
		else if (a.phid > b.phid) return false;
		else {
			if (a.stepfactor <= b.stepfactor) return true;
			else return false;
		}
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
};

enum IOType { ASCII, NETCDF };

class cInputManager {

private:		

	IOType IoType = ASCII;
	cGeophysicsNcFile NC;	

	FILE*  Pointer = (FILE*)NULL;
	std::string FileName;
	size_t HeaderLines=0;
	size_t Subsample=1;

	size_t Record=0;
	std::string RecordString;
	std::vector<std::string> FieldStrings;

public:

	cInputManager() {};

	void initialise(const cBlock& b)
	{		
		FileName = b.getstringvalue("DataFile");
		fixseparator(FileName);		
		std::string ext = extractfileextension(FileName);
		glog.logmsg(0,"Opening Input DataFile %s\n", FileName.c_str());
		if (strcasecmp(ext, ".nc") == 0){			
			IoType = NETCDF;
 			//NC.open(FileName,netCDF::NcFile::FileMode::read);			
			glog.logmsg("Debug 1\n");
			NC.open(FileName, netCDF::NcFile::FileMode::read);
			glog.logmsg("Debug 2\n");			
		}
		else {
			IoType = ASCII;
			HeaderLines = b.getsizetvalue("Headerlines");
			if (!isdefined(HeaderLines)) {
				HeaderLines = 0;
			}					
			Pointer = fileopen(FileName, "r");			
		}

		Subsample = b.getsizetvalue("Subsample");
		if (!isdefined(Subsample)) { Subsample = 1; }		
	}

	IOType iotype() const { return IoType; }

	bool readnextrecord()
	{
		if (iotype() == NETCDF) {
			if (Record == 0)Record++;
			else Record += Subsample;	

			if (Record > NC.ntotalsamples())return false;
		}
		else {
			if (Record == 0) {
				//Skip header lines
				for (size_t i = 0; i < HeaderLines; i++) {
					bool status = filegetline(Pointer, RecordString);
					if (status == false)return status;
					Record++;
				}
			}
			else {
				//Skip lines for subsampling
				for (size_t i = 0; i < Subsample - 1; i++) {
					bool status = filegetline(Pointer, RecordString);
					if (status == false)return status;
					Record++;
				}
			}
			bool status = filegetline(Pointer, RecordString);
			if (status == false) return status;
			Record++;
		}
		return true;
	}

	bool parsefieldstrings() {
		if (iotype() == ASCII) {
			FieldStrings = fieldparsestring(RecordString.c_str(), " ,\t\r\n");
			if (FieldStrings.size() <= 1) return false;
			return true;
		}	
		return true;
	}

	static bool contains_non_numeric_characters(const std::string& str)
	{
		size_t pos = str.find_first_not_of("0123456789.+-eE ,\t\r\n");
		if (pos == std::string::npos) return false;
		else return true;
	}

	const std::string& filename() { return FileName; }

	const size_t& record() { return Record;	}
	
	const std::string& recordstring() const { return RecordString; }

	const std::vector<std::string>& fields() const { return FieldStrings; }

	template<typename T>
	bool netcdf_read(const std::string varname, std::vector<T>& v)
	{
		NC.getDataByPointIndex(varname, Record-1, v);
		return true;
	}

	template<typename T>
	bool read(const cFieldDefinition& cd, T& v) 
	{		
		if (cd.definitiontype() == UNAVAILABLE) {
			v = undefinedvalue(v);
			return false;
		}
		else if (cd.definitiontype() == NUMERIC){
			v = (T) cd.numericvalue[0];
			return true;
		}

		if (iotype() == ASCII) {
			cd.getvalue(fields(),v);
			return true;
		}
		else {	
			std::vector<T> vec;
			netcdf_read(cd.varname, vec);
			v = vec[0];
			if (cd.flip) { v = -1*v; }
			cd.applyoperator(v);
			return true;
		}	
	}

	template<typename T>
	bool read(const cFieldDefinition& cd, std::vector<T>& vec, const size_t n)
	{
		vec.resize(n);
		if (cd.definitiontype() == NUMERIC){
			size_t deflen = cd.numericvalue.size();
			for (size_t i = 0; i < n; i++) {
				if(deflen == 1)vec[i] = (T) cd.numericvalue[0];
				else           vec[i] = (T) cd.numericvalue[i];
			}
			return true;
		}

		if (iotype() == ASCII) {
			cd.getvalue(fields(), vec, n);
			return true;
		}
		else {
			netcdf_read(cd.varname, vec);
			for (size_t i = 0; i < vec.size(); i++) {
				if (cd.flip) { vec[i] = -vec[i]; }
				cd.applyoperator(vec[i]);
			}
			return true;
		}
		return false;
	}
};

class cOutputOptions {

private:
	std::string DumpBasePath;

public:
	std::string DataFile;
	std::string Logfile;
	bool PositiveLayerBottomDepths = false;
	bool NegativeLayerBottomDepths = false;
	bool InterfaceElevations = false;
	bool ParameterSensitivity = false;
	bool ParameterUncertainty = false;
	bool ObservedData = false;
	bool NoiseEstimates = false;
	bool PredictedData = false;
	bool Dump = false;

	std::string DumpPath(const size_t datafilerecord, const size_t iteration){
		return DumpBasePath + pathseparatorstring() + 
			strprint("si%07d", (int)datafilerecord) + pathseparatorstring() +
			strprint("it%03d", (int)iteration) + pathseparatorstring();
	};

	cOutputOptions(){};

	cOutputOptions(const cBlock& b) {
		DataFile = b.getstringvalue("DataFile");
		fixseparator(DataFile);

		Logfile = b.getstringvalue("LogFile");				
		fixseparator(Logfile);

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

	void readdata(cInputManager& IM)
	{
		if (Use == false) return;
		IM.read(fd_oP, oP);
		IM.read(fd_oS, oS, nw());
		oE.resize(nw());
		if (EstimateNoiseFromModel) {
			for (size_t w = 0; w < nw(); w++) {
				const double v = 0.01*mn[w]*oS[w];
				oE[w] = std::hypot(an[w],v);
			}
		}
		else {
			IM.read(fd_oE, oE, nw());
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
	cComponentInfo Comp[3];
	int CompIndex[3];
	int xzIndex;

	bool invertXPlusZ;
	bool invertPrimaryPlusSecondary;
	bool reconstructPrimary;
	//bool estimateNoise;


	//int xIndex  = -1;
	//int yIndex  = -1;
	//int zIndex  = -1;
	//double oPX,oPY,oPZ;
	//std::vector<double>  oSX,oSY,oSZ;	
	//std::vector<double>  oEX,oEY,oEZ;		
	//cFieldDefinition fd_oPX,fd_oPY,fd_oPZ;
	//cFieldDefinition fd_oSX,fd_oSY,fd_oSZ;
	//cFieldDefinition fd_oEX,fd_oEY,fd_oEZ;		
	//double xmn,ymn,zmn;
	//std::vector<double> xan,yan,zan;
	sTDEmData predicted;

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

		Comp[0] = cComponentInfo(b.findblock("XComponent"), "X", nwindows, invertPrimaryPlusSecondary);
		Comp[1] = cComponentInfo(b.findblock("YComponent"), "Y", nwindows, invertPrimaryPlusSecondary);
		Comp[2] = cComponentInfo(b.findblock("ZComponent"), "Z", nwindows, invertPrimaryPlusSecondary);

		ncomps = 0;
		if (Comp[0].Use) ncomps++;
		if (Comp[1].Use) ncomps++;
		if (Comp[2].Use) ncomps++;

		if (invertXPlusZ) {
			Comp[0].Use = true;
			Comp[2].Use = true;
		}
		nchans = nwindows * ncomps;
	}
};

class cSBSInverter{

	public:
	
	cBlock  Control;
	int Size;
	int Rank;
	bool UsingOpenMP;
	
	std::vector<cTDEmSystemInfo> SV;
		
	size_t nsystems;	

	cInputManager IM;
	cOutputOptions OO;	
	FILE* ofp = (FILE*)NULL;
	size_t Outputrecord; //output record number
	
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

	size_t nlayers;
	size_t ndata;
	size_t nparam;
	size_t ngeomparam;	

	std::vector<double> Obs;
	std::vector<double> Err;	
	std::vector<double> Pred;			

	std::vector<double> Param;
	std::vector<double> RefParam;
	std::vector<double> RefParamStd;

	std::vector<double> ParameterSensitivity;
	std::vector<double> ParameterUncertainty;
	
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
	
	
	size_t cIndex;
    size_t tIndex;
	size_t gIndex;	
	size_t tx_heightIndex;
	size_t txrx_dxIndex;
	size_t txrx_dyIndex;
	size_t txrx_dzIndex;
	size_t rx_pitchIndex;
	size_t rx_rollIndex;
	
	cSBSInverter(const std::string& controlfile, const int& size, const int& rank, const bool& usingopenmp)
	{		
		_GSTPUSH_
		try {			
			Size = size;
			Rank = rank;
			UsingOpenMP = usingopenmp;
			initialise(controlfile);			
		}
		catch (const std::string msg) {
			glog.logmsg(msg);
		}
		catch (const std::runtime_error e) {
			glog.logmsg(std::string(e.what()));
		}
		catch (const std::exception e) {
			glog.logmsg(std::string(e.what()));
		}		
	};

	~cSBSInverter()
	{
		if (ofp)fclose(ofp);
		glog.close();
	};

	int go()
	{
		_GSTITEM_
		size_t record = 0;
		while (IM.readnextrecord()) {
			record++;
			if (((record-1) % Size) != Rank) continue;
			if (IM.iotype() == IOType::ASCII) {
				bool nonnumeric = IM.contains_non_numeric_characters(IM.recordstring());
				if (nonnumeric) {
					glog.logmsg("Skipping non-numeric record at line %zu of Input DataFile %s\n", IM.record(), IM.filename().c_str());
					glog.logmsg("\n%s\n\n", IM.recordstring().c_str());
					continue;
				}
			}

			//if (I.OO.Dump){
			//	FILE* fp = fileopen(I.dumppath() + "record.dat", "w");
			//	fprintf(fp, "Record\t%zu", I.DataFileRecord);
			//	fclose(fp);
			//}

			double t1 = gettime();
			invert();

			double t2 = gettime();
			double etime = t2 - t1;
			writeresult();
			glog.logmsg("Rec %6zu  %3zu  %5zu  %10lf  Its=%3zu  PhiD=%6.2lf  time=%.1lfs  %s\n", IM.record(), Id.flightnumber, Id.linenumber, Id.fidnumber, LastIteration, LastPhiD, etime, TerminationReason.c_str());
		}
		glog.close();
		return 0;
	}

	void initialise(const std::string& controlfile)
	{
		_GSTITEM_
		loadcontrolfile(controlfile);
		parsecolumns();
		setup_data();
		setup_parameters();
		resize_matrices();
		go();		
	}

	void loadcontrolfile(const std::string& filename)
	{
		glog.logmsg(0, "Loading control file %s\n", filename.c_str());
		Control = cBlock(filename);
		OO = cOutputOptions(Control.findblock("Output"));
		std::string suffix = stringvalue(Rank, ".%04d");
		OO.Logfile = insert_after_filename(OO.Logfile, suffix);
		OO.DataFile = insert_after_filename(OO.DataFile, suffix);
		openlogfile(); //load this first to get outputlogfile opened	

		//Load control file
		parseoptions();
		initialisesystems();

		IM.initialise(Control.findblock("Input"));

		glog.logmsg(0, "Opening Output DataFile %s\n", OO.DataFile.c_str());
		ofp = fileopen(OO.DataFile, "w");
		Outputrecord = 1;
	}

	void openlogfile()
	{
		glog.logmsg(0, "Opening log file %s\n", OO.Logfile.c_str());
		glog.open(OO.Logfile);
		glog.logmsg(0, "Control file %s\n", Control.Filename.c_str());
		glog.logmsg(0, "Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
		glog.logmsg(0, "Working directory %s\n", getcurrentdirectory().c_str());
		glog.logmsg(0, "Processes=%d\tRank=%d\n", Size, Rank);
		glog.log(Control.get_as_string());
		glog.flush();
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


		NormType = L2;//default
		std::string nt = b.getstringvalue("NormType");
		if (!isdefined(nt)) {
			NormType = L2;
		}
		else if (strcasecmp(nt, "L1") == 0) {
			NormType = L1;
		}
		else if (strcasecmp(nt, "L2") == 0) {
			NormType = L2;
		}
		else {
			glog.errormsg("Unknown NormType %s\n", nt.c_str());
		}


		SmoothnessMethod = SM_2ND_DERIVATIVE;//default
		std::string sm = b.getstringvalue("SmoothnessMethod");
		if (!isdefined(sm)) {
			SmoothnessMethod = SM_2ND_DERIVATIVE;
		}
		else if (strcasecmp(sm, "Minimise1stDerivatives") == 0) {
			SmoothnessMethod = SM_1ST_DERIVATIVE;
		}
		else if (strcasecmp(sm, "Minimize1stDerivatives") == 0) {
			SmoothnessMethod = SM_1ST_DERIVATIVE;
		}
		else if (strcasecmp(sm, "Minimise2ndDerivatives") == 0) {
			SmoothnessMethod = SM_2ND_DERIVATIVE;
		}
		else if (strcasecmp(sm, "Minimize2ndDerivatives") == 0) {
			SmoothnessMethod = SM_2ND_DERIVATIVE;
		}
		else {
			glog.errormsg(_SRC_, "Unknown SmoothnessMethod %s\n", sm.c_str());
		}


		MaxIterations = b.getsizetvalue("MaximumIterations");
		MinimumPhiD = b.getdoublevalue("MinimumPhiD");
		MinimumImprovement = b.getdoublevalue("MinimumPercentageImprovement");
	}

	void parsecolumns()
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

		fd_GI = parsegeometry(b);

		cBlock rm = b.findblock("ReferenceModel");
		fd_GR = parsegeometry(rm);
		fd_ERc.initialise(rm, "Conductivity");
		fd_ERt.initialise(rm, "Thickness");

		cBlock sd = b.findblock("StdDevReferenceModel");
		fd_GS = parsegeometry(sd);
		fd_ESc.initialise(sd, "Conductivity");
		fd_ESt.initialise(sd, "Thickness");

		cBlock tfr = b.findblock("TotalFieldReconstruction");
		fd_GTFR = parsegeometry(tfr);
	}

	std::vector<cFieldDefinition> parsegeometry(const cBlock& b)
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
		ndata = 0;
		for (size_t si = 0; si < nsystems; si++) {
			cTDEmSystemInfo& S = SV[si];
			if (S.invertXPlusZ) {
				S.xzIndex = (int)ndata;
				S.CompIndex[0] = -1;
				S.CompIndex[2] = -1;
				ndata += S.nwindows;

				if (S.Comp[1].Use) {
					S.CompIndex[1] = (int)ndata;
					ndata += S.nwindows;
				}
			}
			else {
				for (size_t i = 0; i < 3; i++) {
					if (S.Comp[i].Use) {
						S.CompIndex[i] = (int)ndata;
						ndata += S.nwindows;
					}
				}
			}
		}
		Obs.resize(ndata);
		Err.resize(ndata);
		Pred.resize(ndata);
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

	void resize_matrices()
	{
		J = MatrixDouble(ndata, nparam);
	}

	bool parserecord()
	{
		if (IM.parsefieldstrings() == false) return false;

		Id.uniqueid = IM.record();
		IM.read(sn, Id.surveynumber);
		IM.read(dn, Id.daynumber);
		IM.read(fn, Id.flightnumber);
		IM.read(ln, Id.linenumber);
		IM.read(fidn, Id.fidnumber);
		IM.read(xord, Location.x);
		IM.read(yord, Location.y);
		IM.read(elevation, Location.groundelevation);
		Location.z = ud_double();

		GI = readgeometry(fd_GI);
		GR = readgeometry(fd_GR);
		GR.fillundefined(GI);
		GTFR = readgeometry(fd_GTFR);
		GTFR.fillundefined(GI);

		GS = readgeometry(fd_GS);

		IM.read(fd_ERc, ER.conductivity, nlayers);
		IM.read(fd_ERt, ER.thickness, nlayers - 1);

		if (solve_conductivity) {
			IM.read(fd_ESc, ES.conductivity, nlayers);
		}
		if (solve_thickness) {
			IM.read(fd_ESt, ES.thickness, nlayers - 1);
		}

		for (size_t si = 0; si < nsystems; si++) {
			readsystemdata(si);
		}
		return true;
	}

	void readsystemdata(size_t sysindex)
	{
		cTDEmSystemInfo& S = SV[sysindex];
		size_t nw = S.nwindows;

		S.Comp[0].readdata(IM);
		S.Comp[1].readdata(IM);
		S.Comp[2].readdata(IM);

		/*
		if (S.useX){
			IM.read(S.fd_oPX, S.oPX);
			IM.read(S.fd_oSX, S.oSX, nw);
			S.oEX.resize(nw);
			if (S.estimateNoise){
				for (size_t w = 0; w < nw; w++){
					const double an = S.xan[w];
					const double mn = 0.01*S.xmn*S.oSX[w];
					S.oEX[w] = sqrt(an*an + mn*mn);
				}
			}
			else{
				IM.read(S.fd_oEX, S.oEX, nw);
			}
		}

		if (S.useY){
			IM.read(S.fd_oPY, S.oPY);
			IM.read(S.fd_oSY, S.oSY, nw);
			S.oEY.resize(nw);
			if (S.estimateNoise){
				for (size_t w = 0; w < nw; w++){
					const double an = S.yan[w];
					const double mn = 0.01*S.ymn*S.oSY[w];
					S.oEY[w] = sqrt(an*an + mn*mn);
				}
			}
			else{
				IM.read(S.fd_oEY, S.oEY, nw);
			}
		}

		if (S.useZ){
			IM.read(S.fd_oPZ, S.oPZ);
			IM.read(S.fd_oSZ, S.oSZ, nw);
			S.oEZ.resize(nw);
			if (S.estimateNoise){
				for (size_t w = 0; w < nw; w++){
					const double an = S.zan[w];
					const double mn = 0.01*S.zmn*S.oSZ[w];
					S.oEZ[w] = sqrt(an*an + mn*mn);
				}
			}
			else{
				IM.read(S.fd_oEZ, S.oEZ, nw);
			}
		}
		*/
	}

	void initialise_sample()
	{
		LastIteration = 0;
		LastLambda = 1e8;

		initialise_data();
		initialise_parameters();
		initialise_Wd();
		initialise_Wc();
		initialise_Wt();
		initialise_Wg();

		if (SmoothnessMethod == SM_1ST_DERIVATIVE) {
			initialise_L_Ws_1st_derivative();
		}
		else if (SmoothnessMethod == SM_2ND_DERIVATIVE) {
			initialise_L_Ws_2nd_derivative();
		}

		initialise_Wr_Wm();

		Param = RefParam;
		LastPhiM = phiModel(Param, LastPhiC, LastPhiT, LastPhiG, LastPhiS);
		TerminationReason = "Has not terminated";

		if (OO.Dump) {
			dumptofile(GR, "geometry_start.dat");
			dumptofile(ER, "earth_start.dat");

			dumptofile(GR, "geometry_ref.dat");
			dumptofile(ER, "earth_ref.dat");

			dumptofile(GS, "geometry_std.dat");
			dumptofile(ES, "earth_std.dat");

			FILE* fp = fileopen(dumppath() + "Id.dat", "w");
			fprintf(fp, "%zu\t%zu\t%zu\t%zu\t%zu\t%lf\t%lf\t%lf\t%lf\t%lf", Id.uniqueid, Id.surveynumber, Id.daynumber, Id.flightnumber, Id.linenumber, Id.fidnumber, Location.x, Location.y, Location.groundelevation, Location.z);
			fclose(fp);
		}
	}

	void initialise_data()
	{
		for (size_t si = 0; si < nsystems; si++) {
			cTDEmSystemInfo& S = SV[si];
			cTDEmSystem& T = S.T;
			if (S.reconstructPrimary) {
				T.setgeometry(GTFR);
				T.LEM.calculation_type = CT_FORWARDMODEL;
				T.LEM.derivative_layer = INT_MAX;
				T.setprimaryfields();
				S.Comp[0].oP = T.PrimaryX;
				S.Comp[1].oP = T.PrimaryY;
				S.Comp[2].oP = T.PrimaryZ;
			}

			if (S.invertXPlusZ) {
				for (size_t wi = 0; wi < S.nwindows; wi++) {

					//X+Z Comp
					size_t di = wi + S.xzIndex;
					if (S.invertPrimaryPlusSecondary) {
						Obs[di] = std::hypot(S.Comp[0].oS[wi] + S.Comp[0].oP,
							S.Comp[2].oS[wi] + S.Comp[2].oP);
					}
					else {
						Obs[di] = std::hypot(S.Comp[0].oS[wi],
							S.Comp[2].oS[wi]);
					}
					Err[di] = std::hypot(S.Comp[0].oE[wi],
						S.Comp[2].oE[wi]);

					//Y Comp
					if (S.Comp[1].Use) {
						di = S.CompIndex[1] + wi;
						Obs[di] = S.Comp[1].oS[wi];
						if (S.invertPrimaryPlusSecondary) {
							Obs[di] += S.Comp[1].oP;
						}
						Err[di] = S.Comp[1].oE[wi];
					}
				}
			}
			else {
				for (size_t ci = 0; ci < 3; ci++) {
					for (size_t wi = 0; wi < S.nwindows; wi++) {
						size_t di = S.CompIndex[ci] + wi;
						Obs[di] = S.Comp[ci].oS[wi];
						if (S.invertPrimaryPlusSecondary) {
							Obs[di] += S.Comp[ci].oP;
						}
						Err[di] = S.Comp[ci].oE[wi];
					}
				}
			}
		}

		if (OO.Dump) {
			dumptofile(Obs, "observed.dat");
			dumptofile(Err, "observed_std.dat");
		}
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
		Wd = MatrixDouble(ndata, ndata, 0.0);
		double s = 1.0 / (double)ndata;
		for (size_t i = 0; i < ndata; i++) {
			Wd[i][i] = s / (Err[i] * Err[i]);
		}
		if (OO.Dump) writetofile(Wd, dumppath() + "Wd.dat");
	}

	void initialise_Wc()
	{
		Wc = MatrixDouble(nparam, nparam, 0.0);
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
			Wc[p][p] = s * (t[i] / tavg) / (RefParamStd[p] * RefParamStd[p]);
		}
	}

	void initialise_Wt()
	{
		Wt = MatrixDouble(nparam, nparam, 0.0);
		if (solve_thickness == false)return;

		double s = AlphaT / (double)(nlayers - 1);
		for (size_t i = 0; i < nlayers - 1; i++) {
			size_t p = i + tIndex;
			Wt[p][p] = s / (RefParamStd[p] * RefParamStd[p]);
		}

	}

	void initialise_Wg()
	{
		Wg = MatrixDouble(nparam, nparam, 0.0);
		if (ngeomparam <= 0)return;

		double s = AlphaG / (double)ngeomparam;
		for (size_t i = 0; i < ngeomparam; i++) {
			size_t p = i + gIndex;
			Wg[p][p] = s / (RefParamStd[p] * RefParamStd[p]);
		}
	}

	void initialise_L_Ws_1st_derivative()
	{
		Ws = MatrixDouble(nparam, nparam, 0.0);
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

		MatrixDouble L = MatrixDouble(nlayers - 1, nparam, 0.0);
		size_t neqn = 0;
		for (size_t li = 1; li < nlayers; li++) {
			size_t pindex = cIndex + li;
			double t1 = t[li - 1];
			double t2 = t[li];
			double d12 = (t1 + t2) / 2.0;
			double s = sqrt(t2 / tavg);//sqrt because it gets squared in L'L		
			L[neqn][pindex - 1] = -s / d12;
			L[neqn][pindex] = s / d12;
			neqn++;
		}
		Ws = transpose_mult(L, L);
		Ws *= (AlphaS / (double)(nlayers - 1));
	}

	void initialise_L_Ws_2nd_derivative()
	{
		Ws = MatrixDouble(nparam, nparam, 0.0);
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

		MatrixDouble L = MatrixDouble(nlayers - 2, nparam, 0.0);
		size_t neqn = 0;
		for (size_t li = 1; li < nlayers - 1; li++) {
			size_t pindex = cIndex + li;
			double t1 = t[li - 1];
			double t2 = t[li];
			double t3 = t[li + 1];
			double d12 = (t1 + t2) / 2.0;
			double d23 = (t2 + t3) / 2.0;
			double s = sqrt(t2 / tavg);//sqrt because it gets squared in L'L		
			L[neqn][pindex - 1] = s / d12;
			L[neqn][pindex] = -s / d12 - s / d23;
			L[neqn][pindex + 1] = s / d23;
			neqn++;
		}
		Ws = transpose_mult(L, L);
		Ws *= (AlphaS / (double)(nlayers - 2));
	}

	void initialise_Wr_Wm()
	{
		Wr = MatrixDouble(nparam, nparam, 0.0);
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

	std::vector<double> parameterchange(const double lambda)
	{
		std::vector<double> x = solve(lambda);
		std::vector<double> dm = x - Param;

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

	std::vector<double> solve(const double lambda)
	{
		// Phi = (d-g(m)+Jm) Wd (d-g(m)+Jm) + lambda ( (m-m0)' Wr (m-m0) + m' Ws m) )
		//Ax = b
		//A = [J'WdJ + lambda (Wr + Ws)]
		//x = m(n+1)
		//b = J'Wd(d - g(m) + Jm) + lambda*Wr*m0
		//dm = m(n+1) - m = x - m

		const std::vector<double>& m = Param;
		const std::vector<double>& d = Obs;
		const std::vector<double>& g = Pred;
		const std::vector<double>& e = Err;
		const std::vector<double>& m0 = RefParam;


		MatrixDouble V = Wd;
		if (NormType == L2) {

		}
		else {
			for (size_t i = 0; i < ndata; i++) {
				const double r = (d[i] - g[i]) / e[i];
				V[i][i] *= 1.0 / std::abs(r);
			}
		}

		MatrixDouble JtV = transpose_mult(J, V);
		MatrixDouble JtVJ = JtV * J;

		std::vector<double> b = JtV * (d - g + J * m) + lambda * Wr*m0;
		MatrixDouble        A = JtVJ + lambda * Wm;
		std::vector<double> x = pseudoinverse_od(A)*b;

		return x;
	}

	double l1_norm(const std::vector<double>& g)
	{
		double l1 = 0.0;
		for (size_t i = 0; i < ndata; i++) {
			l1 += std::abs(Obs[i] - g[i]) / Err[i];
		}
		return l1 / ndata;
	}

	double l2_norm(const std::vector<double>& g)
	{
		std::vector<double> v = Obs - g;
		double l2 = mtDm(v, Wd);
		return l2;
	}

	double phiData(const std::vector<double>& g)
	{
		double phid;
		if (NormType == L1) {
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

	double phiModel(const std::vector<double>& p)
	{
		double phic, phit, phig, phis;
		return phiModel(p, phic, phit, phig, phis);
	}

	double phiModel(const std::vector<double>& p, double& phic, double& phit, double& phig, double& phis)
	{
		phic = phiC(p);
		phit = phiT(p);
		phig = phiG(p);
		phis = phiS(p);

		double v = phic + phit + phig + phis;
		return v;
	}

	double phiC(const std::vector<double>& p)
	{
		if (AlphaC == 0.0)return 0.0;
		if (solve_conductivity == false)return 0.0;
		std::vector<double> v = p - RefParam;
		return mtDm(v, Wc);
	}

	double phiT(const std::vector<double>& p)
	{
		if (AlphaT == 0.0)return 0.0;
		if (solve_thickness == false)return 0.0;
		std::vector<double> v = p - RefParam;
		return mtDm(v, Wt);
	}

	double phiG(const std::vector<double>& p)
	{
		if (AlphaG == 0.0)return 0.0;
		if (ngeomparam == 0)return 0.0;
		std::vector<double> v = p - RefParam;
		return mtDm(v, Wg);
	}

	double phiS(const std::vector<double>& p)
	{
		if (AlphaS == 0)return 0.0;
		else return mtAm(p, Ws);
	}

	cEarth1D get_earth(const std::vector<double>& parameters)
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

	cTDEmGeometry get_geometry(const std::vector<double>& parameters)
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

			sTDEmData& d = S.predicted;
			d.xcomponent.Primary = T.PrimaryX;
			d.ycomponent.Primary = T.PrimaryY;
			d.zcomponent.Primary = T.PrimaryZ;
			d.xcomponent.Secondary = T.X;
			d.ycomponent.Secondary = T.Y;
			d.zcomponent.Secondary = T.Z;
		}
	}

	void forwardmodel(const std::vector<double>& parameters, std::vector<double>& predicted, bool computederivatives)
	{
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
					predicted[wi + S.xzIndex] = xzfm[wi];
					if (S.Comp[1].Use) predicted[wi + S.CompIndex[1]] = yfm[wi];
				}
			}
			else {
				for (size_t wi = 0; wi < nw; wi++) {
					if (S.Comp[0].Use) predicted[wi + S.CompIndex[0]] = xfm[wi];
					if (S.Comp[1].Use) predicted[wi + S.CompIndex[1]] = yfm[wi];
					if (S.Comp[2].Use) predicted[wi + S.CompIndex[2]] = zfm[wi];
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
						fillJacobianColumn(S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
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
						fillJacobianColumn(S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
					}
				}

				if (solve_tx_height) {
					size_t pindex = tx_heightIndex;
					T.LEM.calculation_type = CT_HDERIVATIVE;
					T.LEM.derivative_layer = INT_MAX;
					T.setprimaryfields();
					T.setsecondaryfields();
					fillDerivativeVectors(S, xdrv, ydrv, zdrv);
					fillJacobianColumn(S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}

				if (solve_txrx_dx) {
					size_t pindex = txrx_dxIndex;
					T.LEM.calculation_type = CT_XDERIVATIVE;
					T.LEM.derivative_layer = INT_MAX;
					T.setprimaryfields();
					T.setsecondaryfields();
					fillDerivativeVectors(S, xdrv, ydrv, zdrv);
					fillJacobianColumn(S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}

				if (solve_txrx_dy) {
					size_t pindex = txrx_dyIndex;
					T.LEM.calculation_type = CT_YDERIVATIVE;
					T.LEM.derivative_layer = INT_MAX;
					T.setprimaryfields();
					T.setsecondaryfields();
					fillDerivativeVectors(S, xdrv, ydrv, zdrv);
					fillJacobianColumn(S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}

				if (solve_txrx_dz) {
					size_t pindex = txrx_dzIndex;
					T.LEM.calculation_type = CT_ZDERIVATIVE;
					T.LEM.derivative_layer = INT_MAX;
					T.setprimaryfields();
					T.setsecondaryfields();
					fillDerivativeVectors(S, xdrv, ydrv, zdrv);
					fillJacobianColumn(S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}

				if (solve_rx_pitch) {
					size_t pindex = rx_pitchIndex;
					T.drx_pitch(xfm, zfm, g.rx_pitch, xdrv, zdrv);
					ydrv *= 0.0;
					fillJacobianColumn(S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}

				if (solve_rx_roll) {
					size_t pindex = rx_rollIndex;
					T.drx_roll(yfm, zfm, g.rx_roll, ydrv, zdrv);
					xdrv *= 0.0;
					fillJacobianColumn(S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
				}
			}
		}
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

	void fillJacobianColumn(cTDEmSystemInfo& S, const size_t& pindex, const std::vector<double>& xfm, const std::vector<double>& yfm, const std::vector<double>& zfm, const std::vector<double>& xzfm, const std::vector<double>& xdrv, const std::vector<double>& ydrv, const std::vector<double>& zdrv)
	{
		const size_t nw = S.T.NumberOfWindows;
		if (S.invertXPlusZ) {
			for (size_t w = 0; w < nw; w++) {
				J[w + S.xzIndex][pindex] = (xfm[w] * xdrv[w] + zfm[w] * zdrv[w]) / xzfm[w];
				if (S.Comp[1].Use)J[w + S.CompIndex[1]][pindex] = xdrv[w];
			}
		}
		else {
			for (size_t w = 0; w < nw; w++) {
				for (size_t c = 0; c < 3; c++) {
					if (S.Comp[c].Use)J[w + S.CompIndex[c]][pindex] = xdrv[w];
				}
			}
		}
	}

	std::vector<double> compute_parameter_sensitivity()
	{
		std::vector<double> s(nparam, 0.0);
		for (size_t pi = 0; pi < nparam; pi++) {
			for (size_t di = 0; di < ndata; di++) {
				s[pi] += (fabs(J[di][pi]) * sqrt((double)ndata*Wd[di][di]));
			}
		}

		//if (OO.Dump){
		//	dumptofile(s, "layer_sensitivity.dat");
		//	writetofile(JtWdJ, dumppath() + "JtWdJ.dat");
		//}
		return s;
	}

	std::vector<double> compute_parameter_uncertainty()
	{
		MatrixDouble JtWdJ = transpose_mult(J, Wd)*J;
		MatrixDouble iCm(nparam, nparam, 0.0);
		for (size_t i = 0; i < nparam; i++) iCm[i][i] = 1.0 / (RefParamStd[i] * RefParamStd[i]);
		MatrixDouble X = (double)ndata*JtWdJ + iCm;
		//MatrixDouble X = (double)ndata*JtWdJ + iCm + Ws;	
		//MatrixDouble X = (double)ndata*JtWdJ + Ws;
		//MatrixDouble X = (double)ndata*JtWdJ;
		MatrixDouble pinvX = pseudoinverse(X);
		std::vector<double> s(nparam);
		for (size_t i = 0; i < nparam; i++) {
			s[i] = sqrt(pinvX[i][i]);
		}
		return s;
	}

	void invert()
	{
		parserecord();
		initialise_sample();
		iterate();
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
				writetofile(Obs, dumppath() + "d.dat");
				writetofile(Err, dumppath() + "e.dat");
				writetofile(Pred, dumppath() + "g.dat");
				EM = get_earth(Param);
				GM = get_geometry(Param);
				dumptofile(EM, "earth_inv.dat");
				dumptofile(GM, "geometry_inv.dat");
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
				std::vector<double> dm = parameterchange(t.lambda);
				std::vector<double> mtemp = Param + (t.stepfactor*dm);
				std::vector<double> gtemp(ndata);
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

		if (b == BR_BRACKETED) {
			//bracketed target - find with Brents Method
			double newphid = DBL_MIN;
			double lambda = brentsmethod(T, target, newphid);
			return T.findlambda(lambda);
		}
		else if (b == BR_MINBRACKETED) {
			//bracketed minimum but above target - take smallest phid
			return T.minphidtrial();
		}
		else if (b == BR_ALLBELOW) {
			//all below target	- take the largest lambda			
			return T.maxlambdatrial();
		}
		else if (b == BR_ALLABOVE) {
			//all above target - take smallest phid
			return T.minphidtrial();
		}
		else {
			glog.errormsg(_SRC_, "targetsearch(): Unknown value %d returned from target brackettarget()\n", b);
		}
		return T.minphidtrial();
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
					return BR_BRACKETED;//target bracketed		
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

		std::vector<double> x;
		x.push_back(+1); x.push_back(-1);
		x.push_back(+2); x.push_back(-2);
		x.push_back(+3); x.push_back(-3);
		for (size_t k = 0; k < x.size(); k++) {
			trialfunction(T, pow10(startx + x[k]));
			bool tarbrak = istargetbraketed(T);
			if (tarbrak) {
				return BR_BRACKETED;//target bracketed		
			}
		}
		//printtrials(T);
		//prompttocontinue();

		if (T.maxphid() < target) {
			return BR_ALLBELOW;//all below target	
		}

		bool minbrak = isminbraketed(T);
		if (minbrak)return BR_MINBRACKETED;//min bracketed											

		return BR_ALLABOVE;//all above target
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
				printtrials(T);
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

	void printtrials(cTrialCache T)
	{
		T.sort_lambda();
		printf("\n");
		printf("CurrentLambda = %lf CurrentPhid = %lf    Target = %lf\n", LastLambda, LastPhiD, T.target);
		printf("N    Stepfactor       Lambda          Phid\n");
		for (size_t i = 0; i < T.trial.size(); i++) {
			printf("%2zu %12g %12g %12g\n", T.trial[i].order, T.trial[i].stepfactor, T.trial[i].lambda, T.trial[i].phid);
		}
		printf("\n");
	}

	double goldensearch(double a, double b, double c, double xtol, const double lambda, const std::vector<double>& m, const std::vector<double>& dm, std::vector<double>& g, cTrialCache& cache)
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
			std::vector<double> p = m + x * dm;
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
			std::vector<double> p = m + b * dm;
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
		std::vector<double> dm(nparam);
		std::vector<double> p(nparam);
		std::vector<double> g(ndata);
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

	void writeresult()
	{
		cOutputFileInfo OI;
		std::string buf;

		//Id		
		OI.addfield("uniqueid", 'I', 12, 0);
		OI.setcomment("Inversion sequence number");
		buf += strprint("%12lu", Id.uniqueid);

		OI.addfield("survey", 'I', 12, 0);
		OI.setcomment("Survey number");
		buf += strprint("%12lu", Id.surveynumber);

		OI.addfield("date", 'I', 12, 0);
		OI.setcomment("Date number");
		buf += strprint("%12lu", Id.daynumber);

		OI.addfield("flight", 'I', 12, 0);
		OI.setcomment("Flight number, IntrepidFlightNumber");
		buf += strprint("%12lu", Id.flightnumber);

		OI.addfield("line", 'I', 12, 0);
		OI.setcomment("Line number, IntrepidLineNumber");
		buf += strprint("%12lu", Id.linenumber);

		OI.addfield("fiducial", 'F', 12, 2);
		OI.setcomment("Fiducial number, IntrepidFiducial");
		buf += strprint("%12.2lf", Id.fidnumber);

		//Location
		OI.addfield("easting", 'F', 9, 1);
		OI.setunits("m"); OI.setcomment("IntrepidX");
		buf += strprint("%9.1lf", Location.x);

		OI.addfield("northing", 'F', 10, 1);
		OI.setunits("m"); OI.setcomment("IntrepidY");
		buf += strprint("%10.1lf", Location.y);

		OI.addfield("elevation", 'F', 10, 2);
		OI.setunits("m"); OI.setcomment("Ground elevation relative to sea-level");
		buf += strprint("%10.2lf", Location.groundelevation);

		//Geometry	
		writeresult_geometry(buf, OI, GI, "", "Input ", false);
		writeresult_geometry(buf, OI, GM, "inverted_", "Inverted ", true);

		//Earth	
		OI.addfield("nlayers", 'I', 4, 0);
		OI.setcomment("Number of layers");
		buf += strprint("%4lu", nlayers);

		OI.addfield("conductivity", 'E', 15, 6, nlayers);
		OI.setunits("S/m"); OI.setcomment("Layer conductivity");
		for (size_t i = 0; i < nlayers; i++) {
			buf += strprint("%15.6le", EM.conductivity[i]);
		}

		double bottomlayerthickness = 100.0;
		if (solve_thickness == false && nlayers > 1) {
			bottomlayerthickness = EM.thickness[nlayers - 2];
		}
		std::vector<double> thickness = EM.thickness;
		thickness.push_back(bottomlayerthickness);

		OI.addfield("thickness", 'F', 9, 2, nlayers);
		OI.setunits("m"); OI.setcomment("Layer thickness");
		for (size_t i = 0; i < nlayers; i++) {
			buf += strprint("%9.2lf", thickness[i]);
		}

		if (OO.PositiveLayerBottomDepths) {
			OI.addfield("depth_bottom", 'F', 9, 2, nlayers);
			OI.setunits("m"); OI.setcomment("Depth to bottom of layer");
			double tsum = 0.0;
			for (size_t i = 0; i < nlayers; i++) {
				buf += strprint("%9.2lf", tsum);
				tsum += thickness[i];
			}
		}

		if (OO.NegativeLayerBottomDepths) {
			OI.addfield("depth_bottom_negative", 'F', 9, 2, nlayers);
			OI.setunits("m"); OI.setcomment("Negative of depth to bottom of layer");
			double tsum = 0.0;
			for (size_t i = 0; i < nlayers; i++) {
				tsum += thickness[i];
				buf += strprint("%9.2lf", -tsum);
			}
		}

		if (OO.InterfaceElevations) {
			OI.addfield("elevation_interfaces", 'F', 9, 2, nlayers + 1);
			OI.setunits("m"); OI.setcomment("Elevation of interfaces");
			double etop = Location.groundelevation;
			for (size_t i = 0; i < nlayers; i++) {
				buf += strprint("%9.2lf", etop);
				etop -= thickness[i];
			}
			buf += strprint("%9.2lf", etop);
		}

		if (OO.ParameterSensitivity) {
			if (solve_conductivity) {
				OI.addfield("conductivity_sensitivity", 'E', 15, 6, nlayers);
				for (size_t i = 0; i < nlayers; i++) {
					buf += strprint("%15.6le", ParameterSensitivity[cIndex + i]);
				}
			}
			if (solve_thickness) {
				OI.addfield("thickness_sensitivity", 'E', 15, 6, nlayers - 1);
				for (size_t i = 0; i < nlayers - 1; i++) {
					buf += strprint("%15.6le", ParameterSensitivity[tIndex + i]);
				}
			}
			if (solve_tx_height) {
				OI.addfield("tx_height_sensitivity", 'E', 15, 6);
				buf += strprint("%15.6le", ParameterSensitivity[tx_heightIndex]);
			}
			if (solve_txrx_dx) {
				OI.addfield("txrx_dx_sensitivity", 'E', 15, 6);
				buf += strprint("%15.6le", ParameterSensitivity[txrx_dxIndex]);
			}
			if (solve_txrx_dz) {
				OI.addfield("txrx_dz_sensitivity", 'E', 15, 6);
				buf += strprint("%15.6le", ParameterSensitivity[txrx_dzIndex]);
			}
			if (solve_rx_pitch) {
				OI.addfield("rx_pitch_sensitivity", 'E', 15, 6);
				buf += strprint("%15.6le", ParameterSensitivity[rx_pitchIndex]);
			}
		}

		if (OO.ParameterUncertainty) {
			if (solve_conductivity) {
				OI.addfield("conductivity_uncertainty", 'E', 15, 6, nlayers);
				OI.setunits("log10(S/m)");
				for (size_t i = 0; i < nlayers; i++) {
					buf += strprint("%15.6le", ParameterUncertainty[cIndex + i]);
				}
			}
			if (solve_thickness) {
				OI.addfield("thickness_uncertainty", 'E', 15, 6, nlayers - 1);
				OI.setunits("log10(m)");
				for (size_t i = 0; i < nlayers - 1; i++) {
					buf += strprint("%15.6le", ParameterUncertainty[tIndex + i]);
				}
			}
			if (solve_tx_height) {
				OI.addfield("tx_height_uncertainty", 'E', 15, 6);
				OI.setunits("m");
				buf += strprint("%15.6le", ParameterUncertainty[tx_heightIndex]);
			}
			if (solve_txrx_dx) {
				OI.addfield("txrx_dx_uncertainty", 'E', 15, 6);
				OI.setunits("m");
				buf += strprint("%15.6le", ParameterUncertainty[txrx_dxIndex]);
			}
			if (solve_txrx_dz) {
				OI.addfield("txrx_dz_uncertainty", 'E', 15, 6);
				OI.setunits("m");
				buf += strprint("%15.6le", ParameterUncertainty[txrx_dzIndex]);
			}
			if (solve_rx_pitch) {
				OI.addfield("rx_pitch_uncertainty", 'E', 15, 6);
				OI.setunits("degrees");
				buf += strprint("%15.6le", ParameterUncertainty[rx_pitchIndex]);
			}
		}

		//ObservedData
		if (OO.ObservedData) {
			for (size_t si = 0; si < nsystems; si++) {
				cTDEmSystemInfo& S = SV[si];
				for (auto ci = 0; ci < 3; ci++) {
					if (S.Comp[ci].Use) writeresult_component(buf, OI, si, S.Comp[ci].Name, "observed", "Observed", 'E', 15, 6, S.Comp[ci].oP, S.Comp[ci].oS, S.invertPrimaryPlusSecondary);
				}
			}
		}

		//Noise Estimates
		if (OO.NoiseEstimates) {
			for (size_t si = 0; si < nsystems; si++) {
				cTDEmSystemInfo& S = SV[si];
				for (auto ci = 0; ci < 3; ci++) {
					if (S.Comp[ci].Use) writeresult_component(buf, OI, si, S.Comp[ci].Name, "noise", "Estimated noise", 'E', 15, 6, 0.0, S.Comp[ci].oE, false);
				}
			}
		}

		//PredictedData
		if (OO.PredictedData) {
			for (size_t si = 0; si < nsystems; si++) {
				cTDEmSystemInfo& S = SV[si];
				if (S.Comp[0].Use) writeresult_component(buf, OI, si, "X", "predicted", "Predicted", 'E', 15, 6, S.predicted.xcomponent.Primary, S.predicted.xcomponent.Secondary, S.invertPrimaryPlusSecondary);
				if (S.Comp[1].Use) writeresult_component(buf, OI, si, "Y", "predicted", "Predicted", 'E', 15, 6, S.predicted.ycomponent.Primary, S.predicted.ycomponent.Secondary, S.invertPrimaryPlusSecondary);
				if (S.Comp[2].Use) writeresult_component(buf, OI, si, "Z", "predicted", "Predicted", 'E', 15, 6, S.predicted.zcomponent.Primary, S.predicted.zcomponent.Secondary, S.invertPrimaryPlusSecondary);
			}
		}

		OI.addfield("AlphaC", 'E', 15, 6);
		OI.setcomment("AlphaC inversion parameter");
		buf += strprint("%15.6le", AlphaC);

		OI.addfield("AlphaT", 'E', 15, 6);
		OI.setcomment("AlphaT inversion parameter");
		buf += strprint("%15.6le", AlphaT);

		OI.addfield("AlphaG", 'E', 15, 6);
		OI.setcomment("AlphaG inversion parameter");
		buf += strprint("%15.6le", AlphaG);

		OI.addfield("AlphaS", 'E', 15, 6);
		OI.setcomment("AlphaS inversion parameter");
		buf += strprint("%15.6le", AlphaS);

		OI.addfield("PhiD", 'E', 15, 6);
		OI.setcomment("Normalised data misfit");
		buf += strprint("%15.6le", LastPhiD);

		OI.addfield("PhiM", 'E', 15, 6);
		OI.setcomment("Combined model norm");
		buf += strprint("%15.6le", LastPhiM);

		OI.addfield("PhiC", 'E', 15, 6);
		OI.setcomment("Conductivity model norm");
		buf += strprint("%15.6le", LastPhiC);

		OI.addfield("PhiT", 'E', 15, 6);
		OI.setcomment("Thickness model norm");
		buf += strprint("%15.6le", LastPhiT);

		OI.addfield("PhiG", 'E', 15, 6);
		OI.setcomment("Geometry model norm");
		buf += strprint("%15.6le", LastPhiG);

		OI.addfield("PhiS", 'E', 15, 6);
		OI.setcomment("Smoothness model norm");
		buf += strprint("%15.6le", LastPhiS);

		OI.addfield("Lambda", 'E', 15, 6);
		OI.setcomment("Lambda regularization parameter");
		buf += strprint("%15.6le", LastLambda);

		OI.addfield("Iterations", 'I', 4, 0);
		OI.setcomment("Number of iterations");
		buf += strprint("%4lu", LastIteration);

		//Carriage return		
		buf += strprint("\n");
		fprintf(ofp, buf.c_str());
		fflush(ofp);

		OI.lockfields();
		if (Outputrecord == 1) {
			sFilePathParts fpp = getfilepathparts(OO.DataFile);

			std::string hdrfile = fpp.directory + fpp.prefix + ".hdr";
			OI.write_simple_header(hdrfile);

			std::string aseggdffile = fpp.directory + fpp.prefix + ".dfn";
			OI.write_aseggdf_header(aseggdffile);
		}
		Outputrecord++;
	};

	void writeresult_geometry(std::string& buf, cOutputFileInfo& OI, const cTDEmGeometry& g, const std::string& fieldnameprefix, const std::string& commentprefix, const bool invertedfieldsonly)
	{
		for (size_t i = 0; i < g.size(); i++) {
			if (invertedfieldsonly && solvegeometryindex(i) == false)continue;
			OI.addfield(fieldnameprefix + g.fname(i), 'F', 9, 2);
			OI.setunits(g.units(i));
			OI.setcomment(commentprefix + g.description(i));
			buf += strprint("%9.2lf", g[i]);
		}
	}

	void writeresult_component(std::string& buf, cOutputFileInfo& OI, const size_t& sysnum, const std::string& comp, const std::string& nameprefix, const std::string& commprefix, const char& form, const size_t& width, const size_t& decimals, const double& p, std::vector<double>& s, const bool& includeprimary)
	{
		std::string sysfield = nameprefix + strprint("_EMSystem_%d_", (int)sysnum + 1);
		std::string syscomm = commprefix + strprint(" EMSystem %d ", (int)sysnum + 1);

		std::string fmt;
		if (form == 'F') fmt = strprint("%%%d.%dlf", width, decimals);
		else if (form == 'E') fmt = strprint("%%%d.%dle", width, decimals);
		else {
			glog.errormsg(_SRC_, "Invalid output format %c\n", form);
		}

		if (includeprimary) {
			OI.addfield(sysfield + comp + "P", form, width, decimals);
			OI.setcomment(syscomm + comp + "-component primary field");
			buf += strprint(fmt.c_str(), p);
		}
		OI.addfield(sysfield + comp + "S", form, width, decimals, s.size());
		OI.setcomment(syscomm + comp + "-component secondary field windows");
		for (size_t w = 0; w < s.size(); w++) {
			buf += strprint(fmt.c_str(), s[w]);
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

	cTDEmGeometry readgeometry(const std::vector<cFieldDefinition>& gfd)
	{
		cTDEmGeometry g;
		for (size_t i = 0; i < g.size(); i++) {
			IM.read(gfd[i], g[i]);
		}
		return g;
	}

	void dumptofile(const std::vector<double>& v, std::string path)
	{
		FILE* fp = fileopen(dumppath() + path, "w");
		for (size_t i = 0; i < v.size(); i++) {
			fprintf(fp, "%le\n", v[i]);
		}
		fclose(fp);
	}

	void dumptofile(const cEarth1D& e, std::string path)
	{
		FILE* fp = fileopen(dumppath() + path, "w");
		size_t nl = e.conductivity.size();
		for (size_t i = 0; i < nl; i++) {
			if (i < e.thickness.size()) {
				fprintf(fp, "%e\t%e\n", e.conductivity[i], e.thickness[i]);
			}
			else {
				fprintf(fp, "%e\tInf\n", e.conductivity[i]);
			}
		}
		fclose(fp);
	}

	void dumptofile(const cTDEmGeometry& g, std::string path)
	{
		FILE* fp = fileopen(dumppath() + path, "w");
		fprintf(fp, "tx_height\t%lf\n", g.tx_height);
		fprintf(fp, "tx_roll\t%lf\n", g.tx_roll);
		fprintf(fp, "tx_pitch\t%lf\n", g.tx_pitch);
		fprintf(fp, "tx_yaw\t%lf\n", g.tx_yaw);
		fprintf(fp, "txrx_dx\t%lf\n", g.txrx_dx);
		fprintf(fp, "txrx_dy\t%lf\n", g.txrx_dy);
		fprintf(fp, "txrx_dz\t%lf\n", g.txrx_dz);
		fprintf(fp, "rx_roll\t%lf\n", g.rx_roll);
		fprintf(fp, "rx_pitch\t%lf\n", g.rx_pitch);
		fprintf(fp, "rx_yaw\t%lf\n", g.rx_yaw);
		fclose(fp);
	}
	
	std::string dumppath(){
		return OO.DumpPath(IM.record(),LastIteration);
	};
	

};

#endif
