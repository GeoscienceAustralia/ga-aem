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
#include "tdemsystem.h"

enum eSmoothnessMethod { SM_1ST_DERIVATIVE, SM_2ND_DERIVATIVE };
enum eBracketResult { BR_BRACKETED, BR_MINBRACKETED, BR_ALLABOVE, BR_ALLBELOW };

struct sTrial{
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

struct cTDEmSystemInfo{
	
	public:

	cTDEmSystem T;
	std::string SystemFile;
	size_t nwindows;	
	size_t ncomps;
	size_t nchans;

	bool useX;
	bool useY;
	bool useZ;
	bool useTotal;	
	bool reconstructPrimary;
	bool estimateNoise;


	size_t xIndex,yIndex,zIndex;	
	double          oPX,oPY,oPZ;
	std::vector<double>  oSX,oSY,oSZ;	
	std::vector<double>  oEX,oEY,oEZ;		
	FieldDefinition fd_oPX,fd_oPY,fd_oPZ;
	FieldDefinition fd_oSX,fd_oSY,fd_oSZ;
	FieldDefinition fd_oEX,fd_oEY,fd_oEZ;		
	double xmn,ymn,zmn;
	std::vector<double> xan,yan,zan;
	sTDEmData predicted;
};


class cSBSInverter{

public:
	
	//Constructors
	cSBSInverter(size_t size, size_t rank);				
	~cSBSInverter();
	
	void initialise(const std::string controlfile);
	void loadcontrolfile(const std::string controlfile);	
	void initialisesystems();
	void parseoutputs();
	void parseoptions();
	void parsecolumns();

	bool readnextrecord();
	bool parserecord();

	std::vector<FieldDefinition> parsegeometry(const cBlock& b);		
		
	void setup_data();
	void setup_parameters();
	void resize_matrices();	
		
	
	void replacegeometry(const sTDEmGeometry& a, sTDEmGeometry& b);	
	void initialise_sample();
	void initialise_data();
	void initialise_parameters();
	
	int intvalue(const FieldDefinition& coldef);
	double doublevalue(const FieldDefinition& coldef);
	std::vector<double> doublevector(const FieldDefinition& coldef, const size_t& n);
	std::vector<int> intvector(const FieldDefinition& coldef, const size_t& n);	
	sTDEmGeometry readgeometry(const std::vector<FieldDefinition>& gfd);	
	void readsystemdata(size_t sysindex);
			
	void writeresult();	
		
	//Members
	size_t mRank;
	size_t mSize;		
	FILE*  fp_log;
	cBlock  Control;	

	std::string DataFileName;	
	size_t DataFileHeaderLines;	
	size_t DataFileSubsample;
	FILE*  DataFilePointer;
	size_t DataFileRecord;
	std::string DataFileRecordString;
	std::vector<std::string> DataFileFieldStrings;

	std::vector<cTDEmSystemInfo> SV;
		
	size_t nsystems;	
	std::string InputFile;
	
	std::string OutputDataFile;	
	std::string OutputLogfile;	
	bool OutputPositiveLayerBottomDepths;
	bool OutputNegativeLayerBottomDepths;
	bool OutputInterfaceElevations;
	bool OutputParameterSensitivity;
	bool OutputParameterUncertainty;	
	bool OutputPredictedData;
	bool Dump;
	std::string DumpPath;
	

	FILE* ofp;	
	size_t Outputrecord; //output record number
	
	//column definitions
	FieldDefinition sn, dn, fn, ln, fidn;
	FieldDefinition xord, yord, elevation, altimeter;	
	std::vector<FieldDefinition> fd_GI;
	std::vector<FieldDefinition> fd_GR;
	std::vector<FieldDefinition> fd_GS;	
	std::vector<FieldDefinition> fd_GTFR;	
	FieldDefinition fd_ERc;
	FieldDefinition fd_ERt;
	FieldDefinition fd_ESc;
	FieldDefinition fd_ESt;		
			
	sAirborneSampleId Id;
	sAirborneSampleLocation Location; 	
	sTDEmGeometry GI;
	sTDEmGeometry GM;
	sTDEmGeometry GR;
	sTDEmGeometry GS;
	sTDEmGeometry GTFR;
	cEarth1D EM;
	cEarth1D ER;	
	cEarth1D ES;	

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


	
	size_t nlayers;
	size_t ndata;
	size_t nparam;
	size_t ngeomparam;	

	std::vector<double> vObs;
	std::vector<double> vErr;	
	std::vector<double> vPred;			

	std::vector<double> vParam;
	std::vector<double> vRefParam;
	std::vector<double> vRefParamStd;

	std::vector<double> ParameterSensitivity;
	std::vector<double> ParameterUncertainty;
	
	MatrixDouble J;	
	MatrixDouble JtWd;
	MatrixDouble JtWdJ;	
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
			
	double phiData(const std::vector<double>& g); 	
	double phiModel(const std::vector<double>& p, double& phic, double& phit, double& phig, double& phis);
	double phiModel(const std::vector<double>& p);
	double phiC(const std::vector<double>& p);
	double phiT(const std::vector<double>& p);
	double phiG(const std::vector<double>& p);
	double phiS(const std::vector<double>& p);

	void dumptofile(const std::vector<double>& v, std::string path);
	void dumptofile(const cEarth1D& e, std::string path);
	void dumptofile(const sTDEmGeometry& g, std::string path);
		
	cEarth1D get_earth(const std::vector<double>& parameters);
	sTDEmGeometry get_geometry(const std::vector<double>& parameters);	
	void set_predicted();

	
	void forwardmodel(const std::vector<double>& parameters, std::vector<double>& predicted, bool computederivatives);
	std::vector<double> parameterchange(const double lambda);
	void invert();	
	void iterate();		
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
