/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef galeisbsfdemH
#define galeisbsfdemH

#include "tnt.h"
#include "matrix_ops.h"
#include "asciicolumnfile.h"
#include "fielddefinition.h"
#include "airborne_types.h"
#include "fdemsystem.h"

using namespace TNT;

enum  bracket_t { BRACKETED, MINBRACKETED, ALLABOVE, ALLBELOW };
class cTrial{

	public:
		size_t order;
		double lambda;
		double phid;
		double phim;
		double stepfactor;

		static int lambdacompare(const void *pa, const void *pb)
		{
			const cTrial& a = *((cTrial*)pa);
			const cTrial& b = *((cTrial*)pb);
			if(a.lambda < b.lambda)return -1;	
			else if(a.lambda > b.lambda)return 1;	
			else if(a.lambda == b.lambda){
				if(a.stepfactor < b.stepfactor)return -1;
				else if(a.stepfactor > b.stepfactor)return  1;
			else return 0;
			}
			else return 0;
		}

		static int phidcompare(const void *pa, const void *pb)
		{	
			const cTrial& a = *((cTrial*)pa);	
			const cTrial& b = *((cTrial*)pb);	
			if(a.phid < b.phid)return -1;	
			else if(a.phid == b.phid)return 0;
			else return 1;
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
		qsort(&(trial[0]),trial.size(),sizeof(cTrial),cTrial::lambdacompare);
	}
	
	void sort_phid()
	{				
		qsort(&(trial[0]),trial.size(),sizeof(cTrial),cTrial::phidcompare);
	}
};

class FDEmSampleInverter{

	private:

		class cFDEmSystem EMSystem;		
		size_t mRank;
		size_t mSize;	
		std::string mControlFilename;
		cBlock mControl;

		std::string Logfile;
				
		std::string DumpPath;	
		bool   Dump;	

		std::string DataFileName;
		size_t DataFileHeaderLines;
		size_t DataFileSubsample;
		FILE*  DataFilePointer;
		size_t DataFileRecord;
		std::string DataFileRecordString;
		std::vector<std::string> DataFileFieldStrings;
		

		std::string OutputDataFile;
		std::string OutputHeaderFile;
		bool OutputConductivity;
		bool OutputThickness;
		bool OutputPositiveLayerTopDepths;
		bool OutputNegativeLayerTopDepths;
		bool OutputLayerTopElevations;
		bool OutputObservedData;
		bool OutputPredictedData;
		bool OutputPredictedBirdHeight;
		FILE* ofp;
		FILE* hfp;
		size_t Outputrecord; //output record number
		size_t Outputcolumn; //output column number

		size_t nCoilsets;
		size_t nChannels;		
		size_t nlayers;
		size_t ndata;
		size_t nparam;
		size_t ngeomparam;	

		///
		//column definitions
		cFieldDefinition sn, dn, fn, ln, fidn, time;
		cFieldDefinition xord, yord, elevation, altimeter;	
		cFieldDefinition birdheight;	
		cFieldDefinition birdroll;	
		cFieldDefinition birdpitch;	
		cFieldDefinition birdyaw;	
	
		cFieldDefinition emchannels;
		cFieldDefinition std_emchannels;
	
		cFieldDefinition ref_conductivity;
		cFieldDefinition ref_thickness;
		cFieldDefinition ref_birdheight;	
		cFieldDefinition std_conductivity;
		cFieldDefinition std_thickness;
		cFieldDefinition std_birdheight;	

		///	
		double AlphaC;
		double AlphaT;
		double AlphaG;
		double AlphaS;
		bool thickness_weight_referenceconductivity;
		bool thickness_weight_smoothness;

		bool SolveConductivity;
		bool SolveThickness;
		bool SolveBirdHeight;

		double MinimumPhiD;//overall
		double TargetPhiD;//for an iteration
		double MinimumImprovement;//			
		size_t MaxIterations;
		std::string TerminationReason;

		sAirborneSampleId Id;
		sAirborneSampleLocation Location; 	
		sFDEmData DI;
		sFDEmData DE;
		sFDEmData DM;

		sFDEmGeometry GI;
		sFDEmGeometry GM;
		sFDEmGeometry GR;
		sFDEmGeometry GS;	
		cEarth1D EM;
		cEarth1D ER;	
		cEarth1D ES;	

		size_t cIndex;
		size_t tIndex;	
		size_t gIndex;	
		size_t birdheightIndex;	
		size_t emIndex;
	    
		void initialise();	
		void initialise_emsystems();
		void getoptions();
		void getinputs();
		void getoutputs();
		void parsecolumns();	
		void writeresult();
		void addhdrstring(std::string& hstr, const char* s, size_t nband=1);
		void addhdrstring(std::string& hstr, const std::string s, size_t nband=1);

		cFieldDefinition gcdv(const std::string& fieldname);
		int intvalue(const cFieldDefinition& coldef);
		double doublevalue(const cFieldDefinition& coldef);
		std::vector<double> doublevector(const cFieldDefinition& coldef, const size_t& n);
		std::vector<int> intvector(const cFieldDefinition& coldef, const size_t& n);	
	
		std::vector<double> vObs;
		std::vector<double> vErr;	
		std::vector<double> vPred;	
		std::vector<double> vParam;
		std::vector<double> vRefParam;
		std::vector<double> vRefParamStd;		
		Matrix<double> A;	
		std::vector<double> x;
		std::vector<double> b;
		Matrix<double> J;	
		Matrix<double> JtWd;
		Matrix<double> JtWdJ;	
		Matrix<double> Wd;
		Matrix<double> Wc;
		Matrix<double> Wt;
		Matrix<double> Wg;
		Matrix<double> Wm;
		Matrix<double> L;
		Matrix<double> LtL;
	
		double LastPhiD;//for previous iteration
		double LastPhiM;
		double LastPhiC;
		double LastPhiT;
		double LastPhiG;
		double LastPhiS;
		double LastLambda;
		size_t LastIteration;
		
		void setupdata();
		void setupparameters();
		void resizematrices();
		void fillsample();
		void filldata();
		void fillparameters();
		void fillWd();
		void fillWc();
		void fillWt();
		void fillWg();
		void fillLtL();
		void fillAb(const double lambda);
		void solveAxb();

		double phiData(const std::vector<double>& g); 	
		double phiModel(const std::vector<double>& p, double& phic, double& phit, double& phig, double& phis);
		double phiModel(const std::vector<double>& p);
		double phiC(const std::vector<double>& p);
		double phiT(const std::vector<double>& p);
		double phiG(const std::vector<double>& p);
		double phiS(const std::vector<double>& p);
	
		void dumptofile(const std::vector<double>& v, std::string path);
		void dumptofile(const cEarth1D& e, std::string path);
		void dumptofile(const sFDEmGeometry& g, std::string path);
	
		cEarth1D      get_earth(const std::vector<double>& parameters);
		sFDEmGeometry get_geometry(const std::vector<double>& parameters);
		sFDEmData     get_data(const std::vector<double>& data);
	
		void forwardmodel(const std::vector<double>& parameters, std::vector<double>& predicted, bool computederivatives);    
		void iterate();
	
		std::vector<double> parameterchange(const double lambda);

		cTrial targetsearch(const double currentlambda, const double target);
		double trialfunction(cTrialCache& T, const double triallambda);	
		void   printtrials(cTrialCache T);
		bool   istargetbraketed(cTrialCache& T);
		bool   isminbraketed(cTrialCache& T);
		bool   findhighedge(cTrialCache& T, const double xstart, const double deltax, double& xhighedge);
		bracket_t brackettarget(cTrialCache& T, const double target, const double currentlambda);
		double bisectionsearch(cTrialCache& T, const double target, double& newphid);
		double brentsmethod(cTrialCache& T, const double target, double& newphid);	
		double goldensearch(double a, double b, double c, double xtol, const double lambda, const std::vector<double>& m, const std::vector<double>& dm, std::vector<double>& g, cTrialCache& cache);
		
	public:	
		FDEmSampleInverter(std::string controlfilename, size_t size=1, size_t rank=0);
		bool readnextrecord();
		bool parserecord();	
		void invert();
		void finish();
};

#endif
