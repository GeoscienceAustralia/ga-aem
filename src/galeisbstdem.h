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

#include "blocklanguage.h"
#include "fielddefinition.h"
#include "matrix_ops.h"
#include "airborne_types.h"
#include "tdemsystem.h"
#include "geophysics_netcdf.h"

enum eNormType { L1, L2 };
enum eSmoothnessMethod { SM_1ST_DERIVATIVE, SM_2ND_DERIVATIVE };
enum eBracketResult { BR_BRACKETED, BR_MINBRACKETED, BR_ALLABOVE, BR_ALLBELOW };

class sTrial{

public:
	size_t order;
	double lambda;
	double phid;
	double phim;
	double stepfactor;
};

class cTrialCache{	

	public:
	double target;	
	std::vector<sTrial> trial;
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

	sTrial findlambda(const double lambda){
		for(size_t i=0; i<trial.size(); i++){
			if(trial[i].lambda == lambda){
			  return trial[i];			  
			}
		}
		return trial[0];
	}

	double minphid(){return trial[minphidindex()].phid;}
	double maxphid(){return trial[maxphidindex()].phid;}
	sTrial minphidtrial(){return trial[minphidindex()];}
	sTrial maxphidtrial(){return trial[maxphidindex()];}
	sTrial maxlambdatrial(){return trial[maxlambdaindex()];}

};

class cTDEmSystemInfo{
	
	public:

	cTDEmSystem T;
	std::string SystemFile;
	size_t nwindows;	
	size_t ncomps;
	size_t nchans;

	bool useX=false;
	bool useY=false;
	bool useZ=false;
	bool invertXPlusZ;
	bool invertPrimaryPlusSecondary;	
	bool reconstructPrimary;
	bool estimateNoise;

	int xzIndex = -1;
	int xIndex  = -1;
	int yIndex  = -1;
	int zIndex  = -1;
	double oPX,oPY,oPZ;
	std::vector<double>  oSX,oSY,oSZ;	
	std::vector<double>  oEX,oEY,oEZ;		
	cFieldDefinition fd_oPX,fd_oPY,fd_oPZ;
	cFieldDefinition fd_oSX,fd_oSY,fd_oSZ;
	cFieldDefinition fd_oEX,fd_oEY,fd_oEZ;		
	double xmn,ymn,zmn;
	std::vector<double> xan,yan,zan;
	sTDEmData predicted;
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
 			NC.open(FileName,netCDF::NcFile::FileMode::read);			
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

	const IOType iotype() const { return IoType; }

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
		}
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

class cSBSInverter{

public:
	
	//Constructors
	cSBSInverter(size_t size, size_t rank);				
	~cSBSInverter();
	
	void initialise(const std::string controlfile);
	void loadcontrolfile(const std::string controlfile);		
	void initialisesystems();
	void openlogfile();
	void parseoptions();
	void parsecolumns();
	void check(const cFieldDefinition& cd);
	bool parserecord();

	std::vector<cFieldDefinition> parsegeometry(const cBlock& b);		
		
	void setup_data();
	void setup_parameters();
	void resize_matrices();	
	
	//void replacegeometry(const sTDEmGeometry& a, sTDEmGeometry& b);	
	void initialise_sample();
	void initialise_data();
	void initialise_parameters();
	
	cTDEmGeometry readgeometry(const std::vector<cFieldDefinition>& gfd);	
	void readsystemdata(size_t sysindex);
	bool solvegeometryindex(const size_t index);

	void writeresult();	
	void writeresult_geometry(std::string& buf, cOutputFileInfo& OI, const cTDEmGeometry& g, const std::string& fieldnameprefix, const std::string& commentprefix, const bool invertedfieldsonly);	
	void writeresult_component(std::string& buf, cOutputFileInfo& OI, const size_t& sysnum, const std::string& comp, const std::string& nameprefix, const std::string& commprefix, const char& form, const size_t& width, const size_t& decimals, const double& p, std::vector<double>& s, const bool& includeprimary);
		
	//Members
	size_t mRank;
	size_t mSize;		
	//FILE*   fp_log = (FILE*)NULL;
	cBlock  Control;		
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
		
	
	void initialise_Wd();
	void initialise_Wc();
	void initialise_Wt();
	void initialise_Wg();
	void initialise_L_Ws_1st_derivative();
	void initialise_L_Ws_2nd_derivative();
	void initialise_Wr_Wm();
	std::vector<double> solve(const double lambda);
			
	double l1_norm(const std::vector<double>& g);
	double l2_norm(const std::vector<double>& g);
	double phiData(const std::vector<double>& g); 	
	double phiModel(const std::vector<double>& p, double& phic, double& phit, double& phig, double& phis);
	double phiModel(const std::vector<double>& p);
	double phiC(const std::vector<double>& p);
	double phiT(const std::vector<double>& p);
	double phiG(const std::vector<double>& p);
	double phiS(const std::vector<double>& p);

	std::string dumppath(){
		return OO.DumpPath(IM.record(),LastIteration);
	};

	void save_iteration_file();
	void dumptofile(const std::vector<double>& v, std::string path);
	void dumptofile(const cEarth1D& e, std::string path);
	void dumptofile(const cTDEmGeometry& g, std::string path);
		
	cEarth1D get_earth(const std::vector<double>& parameters);
	cTDEmGeometry get_geometry(const std::vector<double>& parameters);	
	void set_predicted();

	void forwardmodel(const std::vector<double>& parameters, std::vector<double>& predicted, bool computederivatives);
	void fillDerivativeVectors(cTDEmSystemInfo& S, std::vector<double>& xdrv, std::vector<double>& ydrv, std::vector<double>& zdrv);
	void fillJacobianColumn(cTDEmSystemInfo& S, const size_t& pindex, const std::vector<double>& xfm, const std::vector<double>& yfm, const std::vector<double>& zfm, const std::vector<double>& xzfm, const std::vector<double>& xdrv, const std::vector<double>& ydrv, const std::vector<double>& zdrv);
	std::vector<double> parameterchange(const double lambda);
	void invert();	
	void iterate();		
	void iterate_old();
	std::vector<double> compute_parameter_sensitivity();
	std::vector<double> compute_parameter_uncertainty();	
	
	sTrial targetsearch(const double currentlambda, const double target);
	double trialfunction(cTrialCache& T, const double triallambda);	
	void   printtrials(cTrialCache T);
	bool   istargetbraketed(cTrialCache& T);
	bool   isminbraketed(cTrialCache& T);
	bool   findhighedge(cTrialCache& T, const double xstart, const double deltax, double& xhighedge);
	eBracketResult brackettarget(cTrialCache& T, const double target, const double currentlambda);
	double bisectionsearch(cTrialCache& T, const double target, double& newphid);
	double brentsmethod(cTrialCache& T, const double target, double& newphid);	
	double goldensearch(double a, double b, double c, double xtol, const double lambda, const std::vector<double>& m, const std::vector<double>& dm, std::vector<double>& g, cTrialCache& cache);
};

#endif
