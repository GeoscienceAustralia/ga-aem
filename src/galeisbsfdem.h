/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef galeisbsfdemH
#define galeisbsfdemH

#include "eigen_utils.h"
#include "asciicolumnfile.h"
#include "fielddefinition.h"
#include "airborne_types.h"
#include "fdemsystem.h"

enum  bracket_t { BRACKETED, MINBRACKETED, ALLABOVE, ALLBELOW };
class cTrial {

public:
	size_t order;
	double lambda;
	double phid;
	double phim;
	double stepfactor;

	static int lambdacompare(const void* pa, const void* pb)
	{
		const cTrial& a = *((cTrial*)pa);
		const cTrial& b = *((cTrial*)pb);
		if (a.lambda < b.lambda)return -1;
		else if (a.lambda > b.lambda)return 1;
		else if (a.lambda == b.lambda) {
			if (a.stepfactor < b.stepfactor)return -1;
			else if (a.stepfactor > b.stepfactor)return  1;
			else return 0;
		}
		else return 0;
	}

	static int phidcompare(const void* pa, const void* pb)
	{
		const cTrial& a = *((cTrial*)pa);
		const cTrial& b = *((cTrial*)pb);
		if (a.phid < b.phid)return -1;
		else if (a.phid == b.phid)return 0;
		else return 1;
	}
};
class cTrialCache {

public:
	double target;
	std::vector<cTrial> trial;
	double sfsearch(const double& x) {
		for (size_t k = 0; k < trial.size(); ++k) {
			if (x == trial[k].stepfactor)return trial[k].phid;
		}
		return -1.0;
	}

	size_t minphidindex() {
		size_t ind = 0;
		for (size_t k = 1; k < trial.size(); ++k) {
			if (trial[k].phid < trial[ind].phid) {
				ind = k;
			}
		}
		return ind;
	}

	size_t maxphidindex() {
		size_t ind = 0;
		for (size_t k = 1; k < trial.size(); ++k) {
			if (trial[k].phid > trial[ind].phid) {
				ind = k;
			}
		}
		return ind;
	}

	size_t maxlambdaindex() {
		size_t ind = 0;
		for (size_t k = 1; k < trial.size(); ++k) {
			if (trial[k].lambda > trial[ind].lambda) {
				ind = k;
			}
		}
		return ind;
	}

	cTrial findlambda(const double lambda) {
		for (size_t i = 0; i < trial.size(); i++) {
			if (trial[i].lambda == lambda) {
				return trial[i];
			}
		}
		return trial[0];
	}

	double minphid() { return trial[minphidindex()].phid; }
	double maxphid() { return trial[maxphidindex()].phid; }
	cTrial minphidtrial() { return trial[minphidindex()]; }
	cTrial maxphidtrial() { return trial[maxphidindex()]; }
	cTrial maxlambdatrial() { return trial[maxlambdaindex()]; }

	void sort_lambda()
	{
		qsort(&(trial[0]), trial.size(), sizeof(cTrial), cTrial::lambdacompare);
	}

	void sort_phid()
	{
		qsort(&(trial[0]), trial.size(), sizeof(cTrial), cTrial::phidcompare);
	}
};

class FDEmSampleInverter {

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
	FILE* DataFilePointer;
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

	//void addhdrstring(std::string& hstr, const char* s, size_t nband = 1);
	//void addhdrstring(std::string& hstr, const std::string s, size_t nband = 1);

	//cFieldDefinition gcdv(const std::string& fieldname);
	//int intvalue(const cFieldDefinition& coldef);
	//double doublevalue(const cFieldDefinition& coldef);
	//std::vector<double> doublevector(const cFieldDefinition& coldef, const size_t& n);
	//std::vector<int> intvector(const cFieldDefinition& coldef, const size_t& n);

	VectorDouble vObs;
	VectorDouble vErr;
	VectorDouble vPred;
	VectorDouble vParam;
	VectorDouble vRefParam;
	VectorDouble vRefParamStd;
	MatrixDouble A;
	VectorDouble x;
	VectorDouble b;
	MatrixDouble J;
	MatrixDouble JtWd;
	MatrixDouble JtWdJ;
	MatrixDouble Wd;
	MatrixDouble Wc;
	MatrixDouble Wt;
	MatrixDouble Wg;
	MatrixDouble Wm;
	MatrixDouble L;
	MatrixDouble LtL;

	double LastPhiD;//for previous iteration
	double LastPhiM;
	double LastPhiC;
	double LastPhiT;
	double LastPhiG;
	double LastPhiS;
	double LastLambda;
	size_t LastIteration;

	void initialise()
	{
		Logfile = mControl.getstringvalue("LogFile");
		if (mSize != 1) {
			std::string tmp = stringvalue(mRank, ".%04d");
			Logfile = insert_after_filename(Logfile, tmp);
		}

		glog.logmsg("Opening log file %s\n", Logfile.c_str());
		glog.open(Logfile);
		glog.logmsg(0, "Logfile opened on %s\n", timestamp().c_str());
		glog.logmsg(0, "Control file %s\n", mControlFilename.c_str());
		glog.logmsg(0, "Version %s Compiled at %s on %s\n", GAAEM_VERSION, __TIME__, __DATE__);
		glog.logmsg(0, "Working directory %s\n", getcurrentdirectory().c_str());
		glog.logmsg(0, "Processes=%d\tRank=%d\n", mSize, mRank);
		glog.log(mControl.get_as_string());

		//Load control file			
		Dump = mControl.getboolvalue("Dump");
		if (Dump) {
			DumpPath = mControl.getstringvalue("DumpPath");
			fixseparator(DumpPath);
			if (DumpPath[DumpPath.length() - 1] != pathseparator()) {
				DumpPath.append(pathseparatorstring());
			}
			makedirectory(DumpPath.c_str());
		}

		initialise_emsystems();
		getoptions();
		getinputs();
		parsecolumns();
		getoutputs();

		setupdata();
		setupparameters();
		resizematrices();
	}
	void initialise_emsystems()
	{
		cBlock b = mControl.findblock("EMSystem");
		std::string SystemFile = b.getstringvalue("SystemFile");
		EMSystem.readsystemdescriptorfile(SystemFile);
		nCoilsets = EMSystem.NumberOfCoilSets();
		nChannels = nCoilsets * 2;
		glog.log(EMSystem.STM.get_as_string());
	}
	void getoptions()
	{
		cBlock b = mControl.findblock("Earth");
		nlayers = (size_t)b.getintvalue("NumberOfLayers");

		//Options
		b = mControl.findblock("Options");
		SolveConductivity = b.getboolvalue("SolveConductivity");
		SolveThickness = b.getboolvalue("SolveThickness");
		SolveBirdHeight = b.getboolvalue("SolveBirdHeight");

		AlphaC = b.getdoublevalue("AlphaConductivity");
		AlphaT = b.getdoublevalue("AlphaThickness");
		AlphaG = b.getdoublevalue("AlphaGeometry");
		AlphaS = b.getdoublevalue("AlphaSmoothness");

		thickness_weight_referenceconductivity = true;
		thickness_weight_smoothness = true;

		MaxIterations = (size_t)b.getintvalue("MaximumIterations");
		MinimumPhiD = b.getdoublevalue("MinimumPhiD");
		MinimumImprovement = b.getdoublevalue("MinimumPercentageImprovement");
	}
	void getinputs()
	{
		//Input
		cBlock b = mControl.findblock("Input");

		DataFileName = b.getstringvalue("DataFile");

		DataFileHeaderLines = b.getsizetvalue("Headerlines");
		if (!isdefined(DataFileHeaderLines)) {
			DataFileHeaderLines = 0;
		}

		DataFileSubsample = b.getsizetvalue("Subsample");
		if (!isdefined(DataFileSubsample)) {
			DataFileSubsample = 1;
		}

		glog.logmsg(0, "Opening Input DataFile %s\n", DataFileName.c_str());
		DataFilePointer = fileopen(DataFileName, "r");
		DataFileRecord = 0;

		glog.logmsg(0, "Opening Output DataFile %s\n", OutputDataFile.c_str());
		ofp = fileopen(OutputDataFile, "w");
		Outputrecord = 1;
	}

	void getoutputs()
	{
		//Output
		cBlock b = mControl.findblock("Output");
		OutputDataFile = b.getstringvalue("DataFile");
		OutputHeaderFile = b.getstringvalue("HeaderFile");
		OutputConductivity = b.getboolvalue("Conductivity");
		OutputThickness = b.getboolvalue("Thickness");
		OutputPositiveLayerTopDepths = b.getboolvalue("PositiveLayerTopDepths");
		OutputNegativeLayerTopDepths = b.getboolvalue("NegativeLayerTopDepths");
		OutputLayerTopElevations = b.getboolvalue("LayerTopElevations");
		OutputObservedData = b.getboolvalue("ObservedData");
		OutputPredictedData = b.getboolvalue("PredictedData");
		OutputPredictedBirdHeight = b.getboolvalue("PredictedBirdHeight");

		if (mSize > 1) {
			std::string tmp = stringvalue(mRank, ".%04d");
			OutputDataFile = insert_after_filename(OutputDataFile, tmp);
			OutputHeaderFile = insert_after_filename(OutputHeaderFile, tmp);
		}
		glog.logmsg("Opening OutputDataFile file %s\n", OutputDataFile.c_str());
		ofp = fileopen(OutputDataFile, "w");

		glog.logmsg("Opening OutputHeaderFile file %s\n", OutputHeaderFile.c_str());
		hfp = fileopen(OutputHeaderFile, "w");
		Outputrecord = 1;



	}
	void parsecolumns()
	{
		cBlock b = mControl.findblock("Input.Columns");
		sn.initialise(b, "SurveyNumber");
		dn.initialise(b, "DateNumber");
		fn.initialise(b, "FlightNumber");
		ln.initialise(b, "LineNumber");
		fidn.initialise(b, "FidNumber");
		xord.initialise(b, "Easting");
		yord.initialise(b, "Northing");
		elevation.initialise(b, "GroundElevation");
		altimeter.initialise(b, "Altimeter");

		birdheight.initialise(b, "BirdHeight");
		birdroll.initialise(b, "BirdRoll");
		birdpitch.initialise(b, "BirdPitch");
		birdyaw.initialise(b, "BirdYaw");

		emchannels.initialise(b, "EMChannels");
		std_emchannels.initialise(b, "StdDevEMChannels");

		cBlock rm = b.findblock("ReferenceModel");
		ref_conductivity.initialise(rm, "Conductivity");
		ref_thickness.initialise(rm, "Thickness");
		ref_birdheight.initialise(rm, "BirdHeight");

		cBlock sd = b.findblock("StdDevReferenceModel");
		std_conductivity.initialise(sd, "Conductivity");
		std_thickness.initialise(sd, "Thickness");
		std_birdheight.initialise(sd, "BirdHeight");
	}

	void setupdata()
	{
		size_t dindex = 0;
		ndata = 0;

		ndata += nChannels;
		emIndex = dindex;
		dindex += nChannels;

		vObs.resize(ndata);
		vErr.resize(ndata);
		vPred.resize(ndata);
	}

	void setupparameters()
	{
		nparam = 0;
		ngeomparam = 0;

		size_t pindex = 0;
		if (SolveConductivity) {
			cIndex = pindex;
			nparam += nlayers;
			pindex += nlayers;
		}

		if (SolveThickness) {
			tIndex = pindex;
			nparam += nlayers - 1;
			pindex += nlayers - 1;
		}

		if (SolveBirdHeight) {
			gIndex = pindex;
			birdheightIndex = pindex;
			pindex++;
			nparam++;
			ngeomparam++;
		}

		vParam.resize(nparam);
		vRefParam.resize(nparam);
		vRefParamStd.resize(nparam);
	};

	void resizematrices()
	{
		Wc.resize(nparam, nparam);
		Wt.resize(nparam, nparam);
		Wg.resize(nparam, nparam);
		Wm.resize(nparam, nparam);

		Wd.resize(ndata, ndata);
		J.resize(ndata, nparam);

		JtWd.resize(nparam, ndata);
		JtWdJ.resize(nparam, nparam);
		A.resize(nparam, nparam);

		x.resize(nparam);
		b.resize(nparam);
	};

	void writeresult()
	{
		Outputcolumn = 1;
		std::string buf;
		std::string hdr;

		//Id	
		addhdrstring(hdr, "uniqueid");
		buf += strprint("%7d ", Id.uniqueid);

		addhdrstring(hdr, "survey");
		buf += strprint("%4d ", Id.surveynumber);

		addhdrstring(hdr, "date");
		buf += strprint("%8d ", Id.daynumber);

		addhdrstring(hdr, "flight");
		buf += strprint("%4d ", Id.flightnumber);

		addhdrstring(hdr, "line");
		buf += strprint("%10d ", Id.linenumber);

		addhdrstring(hdr, "fid");
		buf += strprint("%10.2lf ", Id.fidnumber);

		//Location
		addhdrstring(hdr, "easting");
		buf += strprint("%8.1lf ", Location.x);

		addhdrstring(hdr, "northing");
		buf += strprint("%9.1lf ", Location.y);

		addhdrstring(hdr, "elevation");
		buf += strprint("%7.2lf ", Location.groundelevation);

		addhdrstring(hdr, "altimeter");
		buf += strprint("%7.2lf ", Location.z);

		//Observed Geometry
		addhdrstring(hdr, "birdheight");
		buf += strprint("%7.2lf ", GI.birdheight);
		addhdrstring(hdr, "birdroll");
		buf += strprint("%7.2lf ", GI.birdroll);
		addhdrstring(hdr, "birdpitch");
		buf += strprint("%7.2lf ", GI.birdpitch);
		addhdrstring(hdr, "birdyaw");
		buf += strprint("%7.2lf ", GI.birdyaw);

		//Predicted TX height
		if (OutputPredictedBirdHeight) {
			addhdrstring(hdr, "predicted_birdheight");
			buf += strprint("%7.2lf ", GM.birdheight);
		}

		//Earth
		addhdrstring(hdr, "nlayers");
		buf += strprint("%3d ", nlayers);

		addhdrstring(hdr, "conductivity", nlayers);
		for (size_t i = 0; i < nlayers; i++) {
			buf += strprint("%14.6le ", EM.conductivity[i]);
		}

		addhdrstring(hdr, "thickness", nlayers - 1);
		for (size_t i = 0; i < nlayers - 1; i++) {
			buf += strprint("%8.2lf ", EM.thickness[i]);
		}

		if (OutputPositiveLayerTopDepths) {
			addhdrstring(hdr, "depth_top", nlayers - 1);
			double tsum = 0.0;
			for (size_t i = 0; i < nlayers - 1; i++) {
				buf += strprint("%8.2lf ", tsum);
				tsum += EM.thickness[i];
			}
		}

		if (OutputNegativeLayerTopDepths) {
			addhdrstring(hdr, "depth_top_negative", nlayers - 1);
			double tsum = 0.0;
			for (size_t i = 0; i < nlayers - 1; i++) {
				buf += strprint("%8.2lf ", -tsum);
				tsum += EM.thickness[i];
			}
		}

		if (OutputLayerTopElevations) {
			addhdrstring(hdr, "elevation_top", nlayers - 1);
			double etop = Location.groundelevation;
			for (size_t i = 0; i < nlayers - 1; i++) {
				buf += strprint("%8.2lf ", etop);
				etop -= EM.thickness[i];
			}
		}

		//Observed
		if (OutputObservedData) {
			addhdrstring(hdr, "observed_data", nChannels);
			for (size_t i = 0; i < nCoilsets; i++) {
				buf += strprint("%8.2lf ", DI.inphase[i]);
				buf += strprint("%8.2lf ", DI.quadrature[i]);
			}
		}

		if (OutputObservedData) {
			addhdrstring(hdr, "observed_noise", nChannels);
			for (size_t i = 0; i < nCoilsets; i++) {
				buf += strprint("%8.2lf ", DE.inphase[i]);
				buf += strprint("%8.2lf ", DE.quadrature[i]);
			}
		}

		//Predicted
		if (OutputPredictedData) {
			addhdrstring(hdr, "predicted_data", nChannels);
			for (size_t i = 0; i < nCoilsets; i++) {
				buf += strprint("%8.2lf ", DM.inphase[i]);
				buf += strprint("%8.2lf ", DM.quadrature[i]);
			}
		}

		addhdrstring(hdr, "AlphaC");
		buf += strprint("%14.6le ", AlphaC);

		addhdrstring(hdr, "AlphaT");
		buf += strprint("%14.6le ", AlphaT);

		addhdrstring(hdr, "AlphaG");
		buf += strprint("%14.6le ", AlphaG);

		addhdrstring(hdr, "AlphaS");
		buf += strprint("%14.6le ", AlphaS);

		addhdrstring(hdr, "PhiD");
		buf += strprint("%14.6le ", LastPhiD);

		addhdrstring(hdr, "PhiM");
		buf += strprint("%14.6le ", LastPhiM);

		addhdrstring(hdr, "PhiC");
		buf += strprint("%14.6le ", LastPhiC);

		addhdrstring(hdr, "PhiT");
		buf += strprint("%14.6le ", LastPhiT);

		addhdrstring(hdr, "PhiG");
		buf += strprint("%14.6le ", LastPhiG);

		addhdrstring(hdr, "PhiS");
		buf += strprint("%14.6le ", LastPhiS);

		addhdrstring(hdr, "Lambda");
		buf += strprint("%14.6le ", LastLambda);

		addhdrstring(hdr, "Iterations");
		buf += strprint("%3d", LastIteration);

		//Carriage return		
		buf += strprint("\n");
		fprintf(ofp, buf.c_str());
		fflush(ofp);

		if (Outputrecord == 1) {
			fprintf(hfp, hdr.c_str());
			fclose(hfp);
		}
		Outputrecord++;

	};
	
	void addhdrstring(std::string& hstr, const std::string s, size_t nband=1)
	{
		addhdrstring(hstr, s.c_str(), nband);
	}

	void addhdrstring(std::string& hstr, const char* s, size_t nband)
	{
		if (nband == 1) {
			hstr += strprint("%d\t%s\n", Outputcolumn, s);
			Outputcolumn++;
		}
		else {
			hstr += strprint("%d-%d\t%s\n", Outputcolumn, Outputcolumn + nband - 1, s);
			Outputcolumn += nband;
		}
	};

	void fillsample()
	{
		filldata();
		fillparameters();
		fillWd();
		fillWc();
		fillWt();
		fillWg();
		fillLtL();

		if (Dump) {
			dumptofile(GI, "geometry_in.dat");
			dumptofile(ER, "earth_in.dat");
		}

		if (Dump) {
			/*
			dumptofile(GR,"geometry_start.dat");
			dumptofile(ER,"earth_start.dat");

			dumptofile(GR,"geometry_ref.dat");
			dumptofile(ER,"earth_ref.dat");

			dumptofile(ES,"earth_std.dat");

			FILE* fp = fileopen(DumpPath+"Id.dat","w");
			fprintf(fp,"%d\t%d\t%d\t%d\t%d\t%lf\t%lf\t%lf\t%lf\t%lf",Id.uniqueid,Id.surveynumber,Id.daynumber,Id.flightnumber,Id.linenumber,Id.fidnumber,Location.x,Location.y,Location.groundelevation,Location.z);
			fclose(fp);
			*/
		}
	}
	
	void filldata()
	{
		//Data  	
		for (size_t i = 0; i < nCoilsets; i++) {
			vObs[emIndex + i * 2] = DI.inphase[i];
			vObs[emIndex + i * 2 + 1] = DI.quadrature[i];
			vErr[emIndex + i * 2] = DE.inphase[i];
			vErr[emIndex + i * 2 + 1] = DE.quadrature[i];
		}

		if (Dump) {
			dumptofile(vObs, "observed.dat");
			dumptofile(vErr, "observed_std.dat");
		}

		if (Dump) {
			dumptofile(vObs, "observed.dat");
			dumptofile(vErr, "observed_std.dat");
		}
	}
	
	void fillparameters()
	{
		if (SolveConductivity) {
			for (size_t i = 0; i < nlayers; i++) {
				vRefParam[i + cIndex] = log10(ER.conductivity[i]);
				vRefParamStd[i + cIndex] = ES.conductivity[i];
			}
		}

		if (SolveThickness) {
			for (size_t i = 0; i < nlayers - 1; i++) {
				vRefParam[i + tIndex] = log10(ER.thickness[i]);
				vRefParamStd[i + tIndex] = ES.thickness[i];
			}
		}

		if (SolveBirdHeight) {
			vRefParam[birdheightIndex] = GI.birdheight;
			vRefParamStd[birdheightIndex] = GS.birdheight;
		}
	}
	
	void fillWd()
	{
		Wd.setZero();
		for (size_t i = 0; i < ndata; i++) {
			Wd(i, i) = 1.0 / (ndata * vErr[i] * vErr[i]);
		}

		//if(Dump) writetofile(Wd,DumpPath+"Wd.dat");
	}
	
	void fillWc()
	{
		Wc.setZero();
		if (SolveConductivity == false)return;

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
			t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3]) * t[nlayers - 2];
		}


		double tsum = 0;
		for (size_t i = 0; i < nlayers; i++)tsum += t[i];

		for (size_t i = 0; i < nlayers; i++) {
			size_t p = i + cIndex;
			if (nlayers > 1 && thickness_weight_referenceconductivity) {
				Wc(p, p) = (double)nlayers * (t[i] / tsum) / (vRefParamStd[p] * vRefParamStd[p]);
			}
			else {
				Wc(p, p) = 1.0 / ((double)nlayers * vRefParamStd[p] * vRefParamStd[p]);
			}
		}
		//if(Dump)writetofile(Wc,DumpPath+"Wc.dat");		
	}
	
	void fillWt()
	{
		Wt.setZero();
		if (SolveThickness == false)return;

		for (size_t i = 0; i < nlayers - 1; i++) {
			size_t p = i + tIndex;
			Wt(p, p) = 1.0 / ((double)(nlayers - 1) * vRefParamStd[p] * vRefParamStd[p]);
		}
		//if(Dump)writetofile(Wt,DumpPath+"Wt.dat");
	}
	
	void fillWg()
	{
		Wg.setZero();
		if (ngeomparam <= 0)return;

		for (size_t i = 0; i < ngeomparam; i++) {
			size_t p = i + gIndex;
			Wg(p, p) = 1.0 / ((double)ngeomparam * vRefParamStd[p] * vRefParamStd[p]);
		}
		//if(Dump)writetofile(Wg,DumpPath+"Wg.dat");	
	}
	
	void fillLtL()
	{
		if (AlphaS == 0 || nlayers < 3) return;
		L = MatrixDouble(nlayers - 2, nparam);
		LtL = MatrixDouble(nparam, nparam);

		L.setZero();

		std::vector<double> t(nlayers);
		if (nlayers == 1) {
			t[0] = 1.0;
		}
		else {
			for (size_t i = 0; i < (nlayers - 1); i++) {
				t[i] = ER.thickness[i];
			}
			t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3]) * t[nlayers - 2];
		}

		if (nlayers > 2 && SolveConductivity) {
			size_t neqn = 0;
			if (thickness_weight_smoothness) {
				for (size_t li = 1; li < nlayers - 1; li++) {
					double t1, t2, t3;
					size_t pindex = cIndex + li;
					t1 = t[li - 1];
					t2 = t[li];
					t3 = t[li + 1];
					L(neqn, pindex - 1) = 1.0 / (0.5 * (t1 + t2));
					L(neqn, pindex) = -1.0 / (0.5 * (t1 + t2)) - 1.0 / (0.5 * (t2 + t3));
					L(neqn, pindex + 1) = 1.0 / (0.5 * (t2 + t3));
					neqn++;
				}
			}
			else {
				for (size_t li = 1; li < nlayers - 1; li++) {
					size_t pindex = cIndex + li;
					L(neqn, pindex - 1) = 1.0;
					L(neqn, pindex) = -2.0;
					L(neqn, pindex + 1) = 1.0;
					neqn++;
				}
			}
		}
		const double s = (1.0 / (double)(nlayers - 2));
		LtL = s * (L.transpose() * L);

		//if(Dump){
		//writetofile(L,DumpPath+"L.dat");		
		//writetofile(LtL,DumpPath+"LtL.dat");
		//}
	}

	void dumptofile(const std::vector<double>& v, std::string path)
	{
		FILE* fp = fileopen(DumpPath + path, "w");
		for (size_t c = 0; c < nCoilsets; c++) {
			fprintf(fp, "%lf\t", v[emIndex + c * 2]);
			fprintf(fp, "%lf\t", v[emIndex + c * 2 + 1]);
			fprintf(fp, "\n");
		}
		fclose(fp);
	}

	void dumptofile(const VectorDouble& v, std::string path)
	{
		FILE* fp = fileopen(DumpPath + path, "w");
		for (size_t c = 0; c < nCoilsets; c++) {
			fprintf(fp, "%lf\t", v[emIndex + c * 2]);
			fprintf(fp, "%lf\t", v[emIndex + c * 2 + 1]);
			fprintf(fp, "\n");
		}
		fclose(fp);
	}

	void dumptofile(const cEarth1D& e, std::string path)
	{
		FILE* fp = fileopen(DumpPath + path, "w");
		size_t nlayers = e.conductivity.size();
		for (size_t i = 0; i < nlayers; i++) {
			if (i < nlayers - 1)fprintf(fp, "%lf\t%lf\n", e.conductivity[i], e.thickness[i]);
			else fprintf(fp, "%lf\t200.0\n", e.conductivity[i]);
		}
		fclose(fp);
	}

	void dumptofile(const sFDEmGeometry& g, std::string path)
	{
		FILE* fp = fileopen(DumpPath + path, "w");
		fprintf(fp, "height\t%lf\n", g.birdheight);
		fclose(fp);
	}

	double phiData(const VectorDouble& g)
	{
		double val;
		double sum = 0.0;
		for (size_t i = 0; i < ndata; i++) {
			val = (vObs[i] - g[i]) / vErr[i];
			sum += (val * val);
		}
		return (sum / ndata);
	}

	double phiModel(const VectorDouble& p)
	{
		double val;
		double sum = 0.0;
		for (size_t i = 0; i < nparam; i++) {
			val = (vRefParam[i] - p[i]) / vRefParamStd[i];
			sum += (val * val);
		}
		return (sum / nparam);
	}

	double phiModel(const VectorDouble& p, double& phic, double& phit, double& phig, double& phis)
	{
		phic = phiC(p);
		phit = phiT(p);
		phig = phiG(p);
		phis = phiS(p);

		double v = 0;
		v += AlphaC * phic;
		v += AlphaT * phit;
		v += AlphaG * phig;
		v += AlphaS * phis;
		return v;
	}
	
	double phiC(const VectorDouble& p)
	{
		if (SolveConductivity == false)return 0.0;
		VectorDouble v = p - vRefParam;
		return mtDm(v, Wc);
	}

	double phiT(const VectorDouble& p)
	{
		if (SolveThickness == false)return 0.0;
		VectorDouble v = p - vRefParam;
		return mtDm(v, Wt);
	}

	double phiG(const VectorDouble& p)
	{
		if (SolveBirdHeight == false)return 0.0;
		VectorDouble v = p - vRefParam;
		return mtDm(v, Wg);
	}
	
	double phiS(const VectorDouble& p)
	{
		if (AlphaS == 0)return 0.0;
		else return mtAm(p, LtL);
	}

	cEarth1D get_earth(const VectorDouble& parameters)
	{
		cEarth1D e = ER;
		if (SolveConductivity) {
			for (size_t li = 0; li < nlayers; li++) {
				e.conductivity[li] = pow10(parameters[li + cIndex]);
			}
		}

		if (SolveThickness) {
			for (size_t li = 0; li < nlayers - 1; li++) {
				e.thickness[li] = pow10(parameters[li + tIndex]);
			}
		}
		return e;
	}
	
	sFDEmGeometry get_geometry(const VectorDouble& parameters)
	{
		sFDEmGeometry g = GR;
		if (SolveBirdHeight) {
			g.birdheight = parameters[birdheightIndex];
		}
		return g;
	}
	
	sFDEmData get_data(const VectorDouble& data)
	{
		sFDEmData d;
		d.inphase.resize(nCoilsets);
		d.quadrature.resize(nCoilsets);
		for (size_t i = 0; i < nCoilsets; i++) {
			d.inphase[i] = data[emIndex + i * 2];
			d.quadrature[i] = data[emIndex + i * 2 + 1];
		}
		return d;
	}

	void forwardmodel(const VectorDouble& parameters, VectorDouble& predicted, bool computederivatives)
	{
		cEarth1D e = get_earth(parameters);
		sFDEmGeometry g = get_geometry(parameters);

		EMSystem.setearth(e);
		EMSystem.setrollpitchyaw(g.birdroll, g.birdpitch, g.birdyaw);
		EMSystem.setheight(g.birdheight);
		EMSystem.setupcomputations();

		//Forwardmodel
		cvector fm;
		fm = EMSystem.ppms();

		for (size_t i = 0; i < nCoilsets; i++) {
			predicted[emIndex + i * 2] = fm[i].real();
			predicted[emIndex + i * 2 + 1] = fm[i].imag();
		}

		if (computederivatives) {
			//Derivatives		
			cvector deriv;

			J.setZero();//initialise to zero

			if (SolveConductivity) {
				for (size_t li = 0; li < nlayers; li++) {
					deriv = EMSystem.dppms(eDC, li);
					size_t pindex = li + cIndex;
					for (size_t i = 0; i < nCoilsets; i++) {
						//multiply by natural log(10) as parameters are in logbase10 units
						J(emIndex + i * 2, pindex) = std::log(10.0) * e.conductivity[li] * deriv[i].real();
						J(emIndex + i * 2 + 1, pindex) = std::log(10.0) * e.conductivity[li] * deriv[i].imag();
					}
				}
			}

			if (SolveThickness) {
				for (size_t li = 0; li < nlayers - 1; li++) {
					deriv = EMSystem.dppms(eDT, li);
					size_t pindex = li + tIndex;
					for (size_t i = 0; i < nCoilsets; i++) {
						//multiply by natural log(10) as parameters are in logbase10 units
						J(emIndex + i * 2, pindex) = std::log(10.0) * e.thickness[li] * deriv[i].real();
						J(emIndex + i * 2 + 1, pindex) = std::log(10.0) * e.thickness[li] * deriv[i].imag();
					}
				}
			}

			if (SolveBirdHeight) {
				deriv = EMSystem.dppms(eDB, 0);
				size_t pindex = birdheightIndex;
				for (size_t i = 0; i < nCoilsets; i++) {
					J(emIndex + i * 2, pindex) = deriv[i].real();
					J(emIndex + i * 2 + 1, pindex) = deriv[i].imag();
				}
			}

			JtWd = J.transpose() * Wd;
			JtWdJ = JtWd * J;

			if (Dump) {
				writetofile(J, DumpPath + "J.dat");
				writetofile(JtWd, DumpPath + "JtWd.dat");
				writetofile(JtWdJ, DumpPath + "JtWdJ.dat");
			}
		}
	}

	void fillAb(const double lambda)
	{
		//Ax = b
		//A = [J'WdJ + lambda*Wm + lanbda*LtL]
		//x = m(n+1)
		//b = J'Wd(d - g(m) + Jm) + labda*Wm*m0
		//dm = m(n+1) - m = x - m

		const VectorDouble& m = vParam;
		const VectorDouble& d = vObs;
		const VectorDouble& g = vPred;
		const VectorDouble& m0 = vRefParam;

		double zc = lambda * AlphaC;
		double zt = lambda * AlphaT;
		double zg = lambda * AlphaG;
		double zs = lambda * AlphaS;

		Wm.setZero();
		if (zc > 0) Wm += zc * Wc;
		if (zt > 0) Wm += zt * Wt;
		if (zg > 0) Wm += zg * Wg;

		b = JtWd * (d - g + J * m) + Wm * m0;
		if (zs > 0) {
			A = JtWdJ + Wm + zs * LtL;
		}
		else {
			A = JtWdJ + Wm;
		}

		if (Dump) {
			writetofile(d, DumpPath + "d.dat");
			writetofile(g, DumpPath + "g.dat");
			writetofile(m, DumpPath + "m.dat");
			writetofile(m0, DumpPath + "m0.dat");
			writetofile(J, DumpPath + "J.dat");
			writetofile(Wc, DumpPath + "Wc.dat");
			writetofile(Wt, DumpPath + "Wt.dat");
			writetofile(A, DumpPath + "A.dat");
			writetofile(b, DumpPath + "b.dat");
			writetofile(A, DumpPath + "A.dat");
		}
	}

	void solveAxb()
	{
		MatrixDouble pinvA = pseudoInverse(A);
		x = pinvA * b;

		//writetofile(A,DumpPath+"A.dat");	
		//writetofile(pinvA,DumpPath+"pinvA.dat");	
		//prompttocontinue();
	}
	
	void iterate()
	{
		VectorDouble dm(nparam);
		VectorDouble gtemp(ndata);
		VectorDouble mtemp(nparam);
		double phidtemp, phimtemp, phictemp, phittemp, phigtemp, phistemp;

		size_t iteration = 0;
		vParam = vRefParam;
		forwardmodel(vParam, vPred, false);
		LastPhiD = phiData(vPred);
		LastPhiM = phiModel(vParam, phictemp, phittemp, phigtemp, phistemp);
		LastPhiC = phictemp;
		LastPhiT = phittemp;
		LastPhiG = phigtemp;
		LastPhiS = phistemp;
		LastLambda = 1e8;
		TargetPhiD = LastPhiD * 0.7;
		LastIteration = 0;


		TerminationReason = "Has not terminated";
		cEarth1D earth = get_earth(vParam);
		sFDEmGeometry geometry = get_geometry(vParam);

		bool keepiterating = true;
		if (LastPhiD < MinimumPhiD) {
			keepiterating = false;
			TerminationReason = "Reached minimum";
		}

		while (keepiterating == true) {
			iteration++;
			TargetPhiD = LastPhiD * 0.7;

			mtemp = vParam;
			gtemp = vPred;

			if (TargetPhiD < MinimumPhiD)TargetPhiD = MinimumPhiD;
			forwardmodel(mtemp, gtemp, true);

			cTrial t = targetsearch(LastLambda, TargetPhiD);

			dm = parameterchange(t.lambda);
			mtemp = vParam + (t.stepfactor * dm);
			forwardmodel(mtemp, gtemp, false);
			phidtemp = phiData(gtemp);
			phimtemp = phiModel(mtemp, phictemp, phittemp, phigtemp, phistemp);
			double percentchange = 100.0 * (LastPhiD - phidtemp) / (LastPhiD);

			if (phidtemp < LastPhiD) {
				vParam = mtemp;
				vPred = gtemp;

				DM = get_data(vPred);
				EM = get_earth(vParam);
				GM = get_geometry(vParam);

				LastPhiD = phidtemp;
				LastPhiM = phimtemp;
				LastPhiC = phictemp;
				LastPhiT = phittemp;
				LastPhiG = phigtemp;
				LastPhiS = phistemp;
				LastLambda = t.lambda;
				LastIteration = iteration;

				if (LastPhiD <= MinimumPhiD) {
					keepiterating = false;
					TerminationReason = "Reached minimum";
				}
				else if (percentchange < MinimumImprovement) {
					keepiterating = false;
					TerminationReason = "Small % improvement";
				}
			}
			else {
				keepiterating = false;
				TerminationReason = "No improvement";
			}

			if (iteration >= MaxIterations) {
				keepiterating = false;
				TerminationReason = "Too many iterations";
			}
		}

		//forwardmodel(vParam,gtemp,true);	
		//LayerSensitivity1 = compute_doi_sensitivity1();	
		//LayerSensitivity2 = compute_doi_sensitivity2();	

		if (Dump) {
			dumptofile(earth, "earth_inv.dat");
			dumptofile(geometry, "geometry_inv.dat");
			dumptofile(vPred, "predicted.dat");
			FILE* fp = fileopen(DumpPath + "iteration.dat", "w");
			fprintf(fp, "Iteration\t%d\n", (int)LastIteration);
			fprintf(fp, "TargetPhiD\t%lf\n", TargetPhiD);
			fprintf(fp, "PhiD\t%lf\n", LastPhiD);
			fprintf(fp, "Lambda\t%lf\n", LastLambda);
			fprintf(fp, "\n");
			fclose(fp);
		}
	}
	
	VectorDouble parameterchange(const double lambda)
	{
		fillAb(lambda);
		solveAxb();
		VectorDouble dm = x - vParam;

		if (SolveConductivity) {
			for (size_t li = 0; li < nlayers; li++) {
				size_t pindex = li + cIndex;
				const double maxcond = 10;
				const double mincond = 1e-6;
				if (vParam[pindex] + dm[pindex] > log10(maxcond)) {
					//printf("upper limit li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
					dm[pindex] = log10(maxcond) - vParam[pindex];
				}
				else if (vParam[pindex] + dm[pindex] < log10(mincond)) {
					//printf("lower limit li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
					dm[pindex] = log10(mincond) - vParam[pindex];
				}
			}
		}

		if (SolveThickness) {
			for (size_t li = 0; li < nlayers - 1; li++) {
				size_t pindex = li + tIndex;
				if (dm[pindex] > 0.5) {
					//printf("li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
					dm[pindex] = 0.5;
				}
				else if (dm[pindex] < -0.5) {
					//printf("li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
					dm[pindex] = -0.5;
				}
			}
		}

		if (SolveBirdHeight) {
			size_t pindex = birdheightIndex;
			if (dm[pindex] > 0.5) {
				//printf("li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] = 0.5;
			}
			else if (dm[pindex] < -0.5) {
				//printf("li=%d pindex=%d dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] = -0.5;
			}

			if (vParam[pindex] + dm[pindex] > 1000) {
				dm[pindex] = 1000 - vParam[pindex];
			}
			else if (vParam[pindex] + dm[pindex] < 10) {
				dm[pindex] = 10 - vParam[pindex];
			}
		}

		//writetofile(x,DumpPath+"x.dat");		
		//writetofile(dm,DumpPath+"dm.dat");		
		return dm;
	}
	
	cTrial targetsearch(const double currentlambda, const double target)
	{
		cTrialCache T;
		T.target = target;
		bracket_t b = brackettarget(T, target, currentlambda);
		//printtrials(T);
		//prompttocontinue();

		if (b == BRACKETED) {
			//bracketed target - find with Brents Method
			double newphid = DBL_MIN;
			double lambda = brentsmethod(T, target, newphid);
			//printtrials(T);
			//prompttocontinue();
			return T.findlambda(lambda);
		}
		else if (b == MINBRACKETED) {
			//bracketed minimum but above target - take smallest phid
			return T.minphidtrial();
		}
		else if (b == ALLBELOW) {
			//all below target	- take the largest lambda			
			return T.maxlambdatrial();
		}
		else if (b == ALLABOVE) {
			//all above target - take smallest phid
			return T.minphidtrial();
		}
		else {
			glog.errormsg(_SRC_, "Error unknown value %d returned from target brackettarget()\n", b);
		}
		return T.minphidtrial();
	}
	
	bool istargetbraketed(cTrialCache& T)
	{
		double target = T.target;
		T.sort_lambda();
		for (size_t i = 0; i < T.trial.size() - 1; i++) {
			if (T.trial[i].phid >= target && T.trial[i + 1].phid <= target) {
				return true;
			}
			if (T.trial[i].phid <= target && T.trial[i + 1].phid >= target) {
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
	
	bracket_t brackettarget(cTrialCache& T, const double target, const double currentlambda)
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
					return BRACKETED;//target bracketed		
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
				return BRACKETED;//target bracketed		
			}
		}
		//printtrials(T);
		//prompttocontinue();

		if (T.maxphid() < target) {
			return ALLBELOW;//all below target	
		}

		bool minbrak = isminbraketed(T);
		if (minbrak)return MINBRACKETED;//min bracketed											

		return ALLABOVE;//all above target
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
			glog.warningmsg(_SRC_, "Target must be bracketed for cSBSInverter::brentsmethod()\n");
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
				glog.warningmsg(_SRC_, "Too many bisections in cSBSInverter::brentsmethod()\n");
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
			printf("%2d %12g %12g %12g\n", (int)T.trial[i].order, T.trial[i].stepfactor, T.trial[i].lambda, T.trial[i].phid);
		}
		printf("\n");
	}
	
	double goldensearch(double a, double b, double c, double xtol, const double lambda, const VectorDouble& m, const VectorDouble& dm, VectorDouble& g, cTrialCache& cache)
	{
		//adapted from http://en.wikipedia.org/wiki/Golden_section_search	
		double resphi = 2 - ((1 + sqrt(5.0)) / 2.0);
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
		VectorDouble g(ndata);
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
		p = vParam + dm;
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
			double gsf = goldensearch(0.0, 0.38196601125010510, 1.0, xtol, triallambda, vParam, dm, g, cache);
			cTrial t3;
			p = vParam + gsf * dm;
			forwardmodel(p, g, false);
			t3.phid = phiData(g);
			t3.phim = phiModel(p);
			t3.stepfactor = gsf;
			t3.lambda = triallambda;
			t3.order = cache.trial.size();
			cache.trial.push_back(t3);
			//double gsf = goldensearch(0.0,0.5,1.0,tau,vParam,dm,g,cache);
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

	int intvalue(const cFieldDefinition& cd)
	{
		int v;
		cd.getvalue(DataFileFieldStrings, v);
		return v;
	}
	
	double doublevalue(const cFieldDefinition& cd)
	{
		double v;
		cd.getvalue(DataFileFieldStrings, v);
		return v;
	}
	
	std::vector<int> intvector(const cFieldDefinition& cd, const size_t& n)
	{
		std::vector<int> v;
		cd.getvalue(DataFileFieldStrings, v, n);
		return v;
	}

	std::vector<double> doublevector(const cFieldDefinition& cd, const size_t& n)
	{
		std::vector<double> v;
		cd.getvalue(DataFileFieldStrings, v, n);
		return v;
	}

public:
	
	FDEmSampleInverter(std::string controlfilename, size_t size, size_t rank)
	{
		mSize = size;
		mRank = rank;
		mControlFilename = controlfilename;
		glog.logmsg("Loading control file %s\n", controlfilename.c_str());
		mControl.loadfromfile(mControlFilename);
		initialise();
	};

	bool readnextrecord() {
		if (DataFileRecord == 0) {
			//Skip header lines
			for (size_t i = 0; i < DataFileHeaderLines; i++) {
				bool status = filegetline(DataFilePointer, DataFileRecordString);
				if (status == false)return status;
				DataFileRecord++;
			}
		}
		else {
			//Skip lines for subsampling
			for (size_t i = 0; i < DataFileSubsample - 1; i++) {
				bool status = filegetline(DataFilePointer, DataFileRecordString);
				if (status == false)return status;
				DataFileRecord++;
			}
		}

		bool status = filegetline(DataFilePointer, DataFileRecordString);
		if (status == false)return status;
		DataFileRecord++;
		return true;
	};

	bool parserecord()
	{
		DataFileFieldStrings = fieldparsestring(DataFileRecordString.c_str(), " ,\t\r\n");
		if (DataFileFieldStrings.size() <= 1)return false;

		Id.uniqueid = DataFileRecord;
		Id.surveynumber = (size_t)intvalue(sn);
		Id.daynumber = (size_t)intvalue(dn);
		Id.flightnumber = (size_t)intvalue(fn);
		Id.linenumber = (size_t)intvalue(ln);
		Id.fidnumber = doublevalue(fidn);
		Location.x = doublevalue(xord);
		Location.y = doublevalue(yord);
		Location.groundelevation = doublevalue(elevation);
		Location.z = doublevalue(altimeter);

		GI.birdheight = doublevalue(birdheight);
		GI.birdroll = doublevalue(birdroll);
		GI.birdpitch = doublevalue(birdpitch);
		GI.birdyaw = doublevalue(birdyaw);

		std::vector<double> chan = doublevector(emchannels, nChannels);
		std::vector<double> stdchan = doublevector(std_emchannels, nChannels);

		DI.inphase.resize(nCoilsets);
		DI.quadrature.resize(nCoilsets);
		DE.inphase.resize(nCoilsets);
		DE.quadrature.resize(nCoilsets);
		for (size_t i = 0; i < nCoilsets; i++) {
			DI.inphase[i] = chan[i * 2];
			DI.quadrature[i] = chan[i * 2 + 1];
			DE.inphase[i] = stdchan[i * 2];
			DE.quadrature[i] = stdchan[i * 2 + 1];
		}

		ER.conductivity = doublevector(ref_conductivity, nlayers);
		if (SolveConductivity) {
			ES.conductivity = doublevector(std_conductivity, nlayers);
		}

		ER.thickness = doublevector(ref_thickness, nlayers - 1);
		if (SolveThickness) {
			ES.thickness = doublevector(std_thickness, nlayers - 1);
		}

		GR = GI;
		if (SolveBirdHeight) {
			GR.birdheight = doublevalue(ref_birdheight);
			GS.birdheight = doublevalue(std_birdheight);
		}

		return true;
	};
	
	void invert()
	{
		if (Dump) {
			FILE* fp = fileopen(DumpPath + "record.dat", "w");
			fprintf(fp, "Record\t%zu", DataFileRecord);
			fclose(fp);
		}

		double t1 = gettime();
		fillsample();
		iterate();
		writeresult();
		double t2 = gettime();

		double etime = t2 - t1;
		glog.logmsg("Rec %6lu\t %3lu\t %5lu\t %10lf ...", DataFileRecord, Id.flightnumber, Id.linenumber, Id.fidnumber);
		glog.logmsg("Its=%3lu\tPhiD=%6.2lf\t%s time=%.1lfs\n", LastIteration, LastPhiD, TerminationReason.c_str(), etime);
	}

	void finish()
	{
		fclose(ofp);
		glog.close();
	}
};

#endif
