/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cstdio>
#include <cfloat>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <vector>
#include <random>



#include "file_utils.h"
#include "file_formats.h"
#include "lem.h"
#include "tdemsystem.h"
#include "matrix_ops.h"
#include "vector_utils.h"
#include "galeisbstdem.h"
#include "stacktrace.h"
#include "large_loop.h"

cStackTrace globalstacktrace;

#define VERSION "1.0"

#if defined _MPI_ENABLED
	#include "mpi_wrapper.h"
#endif

#if defined _OPENMP
	#include <omp.h>
	//This thread lock must be set when fftw is being initialised
	omp_lock_t fftw_thread_lock;
#endif

int triallambdacompare(const void *pa, const void *pb)
{
	const sTrial& a = *((sTrial*)pa);
	const sTrial& b = *((sTrial*)pb);
	if (a.lambda < b.lambda)return -1;
	else if (a.lambda > b.lambda)return 1;
	else if (a.lambda == b.lambda){
		if (a.stepfactor < b.stepfactor)return -1;
		else if (a.stepfactor > b.stepfactor)return  1;
		else return 0;
	}
	else return 0;
}
int trialphidcompare(const void *pa, const void *pb)
{
	const sTrial& a = *((sTrial*)pa);
	const sTrial& b = *((sTrial*)pb);
	if (a.phid < b.phid)return -1;
	else if (a.phid == b.phid)return 0;
	else return 1;
}
void sorttrial_lambda(cTrialCache& x)
{
	qsort(&(x.trial[0]), x.trial.size(), sizeof(sTrial), triallambdacompare);
}
void sorttrial_phid(cTrialCache& x)
{
	qsort(&(x.trial[0]), x.trial.size(), sizeof(sTrial), trialphidcompare);
}

cSBSInverter::cSBSInverter(size_t size, size_t rank)
{	
	mSize = size;
	mRank = rank;
}
cSBSInverter::~cSBSInverter()
{
	fclose(fp_log);
	fclose(ofp);
};
void cSBSInverter::initialise(const std::string controlfile)
{
	_GSTPUSH_
	loadcontrolfile(controlfile);
	parsecolumns();
	setup_data();
	setup_parameters();
	resize_matrices();
	_GSTPOP_
}
void cSBSInverter::loadcontrolfile(const std::string filename)
{
	rootmessage("Loading control file %s\n", filename.c_str());
	Control = cBlock(filename);

	cBlock OB = Control.findblock("Output");
	OO = cOutputOptions(OB);	
	std::string suffix = stringvalue(mRank, ".%04lu");
	OO.Logfile  = insert_after_filename(OO.Logfile, suffix);
	OO.DataFile = insert_after_filename(OO.DataFile, suffix);
	openlogfile(); //load this first to get outputlogfile opened	
	
	//Load control file				
	parseoptions();
	initialisesystems();

	cBlock IB = Control.findblock("Input");	
	DataFileName  = IB.getstringvalue("DataFile");
	fixseparator(DataFileName);

	DataFileHeaderLines = IB.getsizetvalue("Headerlines");		
	if (!isdefined(DataFileHeaderLines)){
		DataFileHeaderLines = 0;
	}
	
	DataFileSubsample   = IB.getsizetvalue("Subsample");
	if (!isdefined(DataFileSubsample)){
		DataFileSubsample = 1;
	}	
	
	rootmessage(fp_log, "Opening Input DataFile %s\n", DataFileName.c_str());
	DataFilePointer = fileopen(DataFileName, "r");
	DataFileRecord  = 0;
	
	rootmessage(fp_log,"Opening Output DataFile %s\n", OO.DataFile.c_str());
	ofp = fileopen(OO.DataFile, "w");	
	Outputrecord = 1;
}
void cSBSInverter::openlogfile()
{			
	rootmessage("Opening log file %s\n", OO.Logfile.c_str());
	fp_log = fileopen(OO.Logfile, "w");
	rootmessage(fp_log, "Logfile opened on %s\n", timestamp().c_str());
	rootmessage(fp_log, "Control file %s\n", Control.Filename.c_str());
	rootmessage(fp_log, "Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
	rootmessage(fp_log, "Working directory %s\n", getcurrentdirectory().c_str());
	rootmessage(fp_log, "Processes=%lu\tRank=%lu\n", mSize, mRank);
	Control.write(fp_log);	
}
void cSBSInverter::parseoptions()
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
	
	SmoothnessMethod = SM_2ND_DERIVATIVE;//default
	std::string sm = b.getstringvalue("SmoothnessMethod");
	if (!isdefined(sm)){
		SmoothnessMethod = SM_2ND_DERIVATIVE;
	}
	else if (strcasecmp(sm,"Minimise1stDerivatives")==0){
		SmoothnessMethod = SM_1ST_DERIVATIVE;
	}
	else if (strcasecmp(sm, "Minimize1stDerivatives") == 0){
		SmoothnessMethod = SM_1ST_DERIVATIVE;
	}
	else if (strcasecmp(sm, "Minimise2ndDerivatives") == 0){
		SmoothnessMethod = SM_2ND_DERIVATIVE;
	}
	else if (strcasecmp(sm, "Minimize2ndDerivatives") == 0){
		SmoothnessMethod = SM_2ND_DERIVATIVE;
	}
	else{
		errormessage(fp_log,"Unknown SmoothnessMethod %s\n", sm.c_str());
	}
		

	MaxIterations = b.getsizetvalue("MaximumIterations");
	MinimumPhiD   = b.getdoublevalue("MinimumPhiD");
	MinimumImprovement = b.getdoublevalue("MinimumPercentageImprovement");
}
void cSBSInverter::parsecolumns()
{	
	cBlock b = Control.findblock("Input.Columns");
	sn.set(b, "SurveyNumber");
	dn.set(b, "DateNumber");
	fn.set(b, "FlightNumber");
	ln.set(b, "LineNumber");
	fidn.set(b, "FidNumber");
	xord.set(b, "Easting");
	yord.set(b, "Northing");
	elevation.set(b, "GroundElevation");	

	fd_GI = parsegeometry(b);

	cBlock rm = b.findblock("ReferenceModel");
	fd_GR = parsegeometry(rm);
	fd_ERc.set(rm, "Conductivity");
	fd_ERt.set(rm, "Thickness");

	cBlock sd = b.findblock("StdDevReferenceModel");
	fd_GS = parsegeometry(sd);
	fd_ESc.set(sd, "Conductivity");
	fd_ESt.set(sd, "Thickness");

	cBlock tfr = b.findblock("TotalFieldReconstruction");
	fd_GTFR = parsegeometry(tfr);
}
std::vector<FieldDefinition> cSBSInverter::parsegeometry(const cBlock& b)
{
	std::vector<FieldDefinition> g(10);
	for (size_t i = 0; i < g.size(); i++) {
		g[i].set(b, cTDEmGeometry::fname(i));
	}	
	return g;
}
void cSBSInverter::initialisesystems()
{
	nsystems = Control.getsizetvalue("NumberOfSystems");
	SV.resize(nsystems);
	for (size_t si = 0; si < nsystems; si++){
		cTDEmSystemInfo& S = SV[si];		
		std::string str = strprint("EMSystem%lu", si + 1);
		cBlock  b = Control.findblock(str);

		cTDEmSystem& T = S.T;
		std::string stmfile = b.getstringvalue("SystemFile");
		fixseparator(stmfile);

		rootmessage("Reading system file %s\n", stmfile.c_str());
		T.readsystemdescriptorfile(stmfile);
		fprintf(fp_log, "==============System file %s\n", stmfile.c_str());
		T.STM.write(fp_log);
		fprintf(fp_log, "==========================================================================\n");
		S.nwindows = T.NumberOfWindows;

		S.useX = b.getboolvalue("UseXComponent");
		S.useY = b.getboolvalue("UseYComponent");
		S.useZ = b.getboolvalue("UseZComponent");
		S.useTotal = b.getboolvalue("InvertTotalField");
		S.reconstructPrimary = b.getboolvalue("ReconstructPrimaryFieldFromInputGeometry");
		S.estimateNoise = b.getboolvalue("EstimateNoiseFromModel");
		S.ncomps = 0;
		if (S.useX){
			S.ncomps++;
			if (S.estimateNoise){
				S.xmn = b.getdoublevalue("XMultiplicativeNoise");
				S.xan = b.getdoublevector("XAdditiveNoise");
			}
			S.oSX.resize(S.nwindows);
			S.oEX.resize(S.nwindows);

			S.fd_oPX.set(b, "XComponentPrimary");
			S.fd_oSX.set(b, "XComponentSecondary");
			S.fd_oEX.set(b, "XComponentNoise");
		}

		if (S.useY){
			S.ncomps++;
			if (S.estimateNoise){
				S.ymn = b.getdoublevalue("YMultiplicativeNoise");
				S.yan = b.getdoublevector("YAdditiveNoise");
			}
			S.oSY.resize(S.nwindows);
			S.oEY.resize(S.nwindows);

			S.fd_oPY.set(b, "YComponentPrimary");
			S.fd_oSY.set(b, "YComponentSecondary");
			S.fd_oEY.set(b, "YComponentNoise");
		}
		if (S.useZ){
			S.ncomps++;
			if (S.estimateNoise){
				S.zmn = b.getdoublevalue("ZMultiplicativeNoise");
				S.zan = b.getdoublevector("ZAdditiveNoise");
			}
			S.oSZ.resize(S.nwindows);
			S.oEZ.resize(S.nwindows);

			S.fd_oPZ.set(b, "ZComponentPrimary");
			S.fd_oSZ.set(b, "ZComponentSecondary");
			S.fd_oEZ.set(b, "ZComponentNoise");
		}
		S.nchans = S.nwindows*S.ncomps;
	}
}
void cSBSInverter::setup_data()
{
	ndata = 0;
	for (size_t si = 0; si < nsystems; si++){
		cTDEmSystemInfo& S = SV[si];
		if (S.useX){
			S.xIndex = ndata;
			ndata += S.nwindows;
		}
		if (S.useY){
			S.yIndex = ndata;
			ndata += S.nwindows;
		}
		if (S.useZ){
			S.zIndex = ndata;
			ndata += S.nwindows;
		}
	}
	vObs.resize(ndata);
	vErr.resize(ndata);
	vPred.resize(ndata);
}
void cSBSInverter::setup_parameters()
{
	nlayers = Control.getsizetvalue("Earth.NumberOfLayers");
	nparam = 0;
	ngeomparam = 0;
	if (solve_conductivity){
		cIndex = nparam;
		nparam += nlayers;
	}

	if (solve_thickness){
		tIndex = nparam;
		nparam += nlayers - 1;
	}


	//Geometry params
	gIndex = nparam;
	if (solve_tx_height){
		tx_heightIndex = nparam;
		nparam++;
		ngeomparam++;
	}

	if (solve_txrx_dx){
		txrx_dxIndex = nparam;
		nparam++;
		ngeomparam++;
	}

	if (solve_txrx_dy){
		txrx_dyIndex = nparam;
		nparam++;
		ngeomparam++;
	}

	if (solve_txrx_dz){
		txrx_dzIndex = nparam;
		nparam++;
		ngeomparam++;
	}

	if (solve_rx_pitch){
		rx_pitchIndex = nparam;
		nparam++;
		ngeomparam++;
	}

	if (solve_rx_roll){
		rx_rollIndex = nparam;
		nparam++;
		ngeomparam++;
	}

	///////////////
	vParam.resize(nparam);
	vRefParam.resize(nparam);
	vRefParamStd.resize(nparam);
}
void cSBSInverter::resize_matrices()
{	
	J     = MatrixDouble(ndata, nparam);
	JtWd  = MatrixDouble(nparam, ndata);
	JtWdJ = MatrixDouble(nparam, nparam);	
}
bool cSBSInverter::readnextrecord(){
	if (DataFileRecord == 0){
		//Skip header lines
		for (size_t i = 0; i < DataFileHeaderLines; i++){
			bool status = filegetline(DataFilePointer, DataFileRecordString);
			if (status == false)return status;
			DataFileRecord++;
		}
	}
	else{
		//Skip lines for subsampling
		for (size_t i = 0; i < DataFileSubsample - 1; i++){
			bool status = filegetline(DataFilePointer, DataFileRecordString);
			if (status == false)return status;
			DataFileRecord++;
		}
	}
		
	bool status = filegetline(DataFilePointer, DataFileRecordString);	
	if (status == false)return status;	
	DataFileRecord++;
	return true;
}
bool cSBSInverter::parserecord()
{				
	DataFileFieldStrings = fieldparsestring(DataFileRecordString.c_str(), " ,\t\r\n");
	if (DataFileFieldStrings.size() <= 1)return false;

	Id.uniqueid = DataFileRecord;

	Id.surveynumber = (size_t)intvalue(sn);
	Id.daynumber    = (size_t)intvalue(dn);
	Id.flightnumber = (size_t)intvalue(fn);
	Id.linenumber   = (size_t)intvalue(ln);
	Id.fidnumber    = doublevalue(fidn);
	Location.x      = doublevalue(xord);
	Location.y      = doublevalue(yord);
	Location.groundelevation = doublevalue(elevation);
	Location.z      = ud_double();

	GI = readgeometry(fd_GI);	
	GR = readgeometry(fd_GR);
	GR.fillundefined(GI);	
	GTFR = readgeometry(fd_GTFR);
	GTFR.fillundefined(GI);

	GS = readgeometry(fd_GS);

	ER.conductivity = doublevector(fd_ERc, nlayers);
	ER.thickness = doublevector(fd_ERt, nlayers - 1);

	if (solve_conductivity){
		ES.conductivity = doublevector(fd_ESc, nlayers);
	}
	if (solve_thickness){
		ES.thickness = doublevector(fd_ESt, nlayers - 1);
	}

	for (size_t si = 0; si < nsystems; si++){
		readsystemdata(si);
	}

	initialise_sample();
	return true;
}
void cSBSInverter::readsystemdata(size_t sysindex)
{
	cTDEmSystemInfo& S = SV[sysindex];
	size_t nw = S.nwindows;

	if (S.useX){
		S.oPX = doublevalue(S.fd_oPX);
		S.oSX = doublevector(S.fd_oSX, nw);
		S.oEX.resize(nw);
		if (S.estimateNoise){
			for (size_t w = 0; w < nw; w++){
				const double an = S.xan[w];
				const double mn = 0.01*S.xmn*S.oSX[w];
				S.oEX[w] = sqrt(an*an + mn*mn);
			}
		}
		else{
			S.oEX = doublevector(S.fd_oEX, nw);
		}
	}

	if (S.useY){
		S.oPY = doublevalue(S.fd_oPY);
		S.oSY = doublevector(S.fd_oSY, nw);
		S.oEY.resize(nw);
		if (S.estimateNoise){
			for (size_t w = 0; w < nw; w++){
				const double an = S.yan[w];
				const double mn = 0.01*S.ymn*S.oSY[w];
				S.oEY[w] = sqrt(an*an + mn*mn);
			}
		}
		else{
			S.oEY = doublevector(S.fd_oEY, nw);
		}
	}

	if (S.useZ){
		S.oPZ = doublevalue(S.fd_oPZ);
		S.oSZ = doublevector(S.fd_oSZ, nw);
		S.oEZ.resize(nw);
		if (S.estimateNoise){
			for (size_t w = 0; w < nw; w++){
				const double an = S.zan[w];
				const double mn = 0.01*S.zmn*S.oSZ[w];
				S.oEZ[w] = sqrt(an*an + mn*mn);
			}
		}
		else{
			S.oEZ = doublevector(S.fd_oEZ, nw);
		}
	}
}
void cSBSInverter::initialise_sample()
{
	initialise_data();
	initialise_parameters();
	initialise_Wd();
	initialise_Wc();
	initialise_Wt();
	initialise_Wg();

	if (SmoothnessMethod == SM_1ST_DERIVATIVE){
		initialise_L_Ws_1st_derivative();
	}
	else if (SmoothnessMethod == SM_2ND_DERIVATIVE){
		initialise_L_Ws_2nd_derivative();
	}

	initialise_Wr_Wm();

	if (OO.Dump){
		dumptofile(GR, "geometry_start.dat");
		dumptofile(ER, "earth_start.dat");

		dumptofile(GR, "geometry_ref.dat");
		dumptofile(ER, "earth_ref.dat");

		dumptofile(ES, "earth_std.dat");

		FILE* fp = fileopen(OO.DumpPath + "Id.dat", "w");
		fprintf(fp, "%lu\t%lu\t%lu\t%lu\t%lu\t%lf\t%lf\t%lf\t%lf\t%lf", Id.uniqueid, Id.surveynumber, Id.daynumber, Id.flightnumber, Id.linenumber, Id.fidnumber, Location.x, Location.y, Location.groundelevation, Location.z);
		fclose(fp);
	}
}
void cSBSInverter::initialise_data()
{
	for (size_t si = 0; si < nsystems; si++){
		cTDEmSystemInfo& S = SV[si];
		cTDEmSystem& T = S.T;
		if (S.reconstructPrimary){
			T.setgeometry(GTFR);
			T.LEM.calculation_type = CT_FORWARDMODEL;
			T.LEM.derivative_layer = INT_MAX;
			T.setprimaryfields();
			S.oPX = T.PrimaryX;
			S.oPY = T.PrimaryY;
			S.oPZ = T.PrimaryZ;
		}

		if (S.useX){
			for (size_t wi = 0; wi < S.nwindows; wi++){
				size_t di = wi + S.xIndex;
				vErr[di] = S.oEX[wi];
				vObs[di] = S.oSX[wi];
				if (S.useTotal)vObs[di] += S.oPX;
			}
		}
		if (S.useY){
			for (size_t wi = 0; wi < S.nwindows; wi++){
				size_t di = wi + S.yIndex;
				vErr[di] = S.oEY[wi];
				vObs[di] = S.oSY[wi];
				if (S.useTotal)vObs[di] += S.oPY;				
			}
		}
		if (S.useZ){
			for (size_t wi = 0; wi < S.nwindows; wi++){
				size_t di = wi + S.zIndex;
				vErr[di] = S.oEZ[wi];
				vObs[di] = S.oSZ[wi];
				if (S.useTotal)vObs[di] += S.oPZ;				
			}
		}
	}

	if (OO.Dump){
		dumptofile(vObs, "observed.dat");
		dumptofile(vErr, "observed_std.dat");
	}
}
void cSBSInverter::initialise_parameters()
{
	if (solve_conductivity){
		for (size_t i = 0; i < nlayers; i++){
			vRefParam[i + cIndex] = log10(ER.conductivity[i]);
			vRefParamStd[i + cIndex] = ES.conductivity[i];
		}
	}

	if (solve_thickness){
		for (size_t i = 0; i < nlayers - 1; i++){
			vRefParam[i + tIndex] = log10(ER.thickness[i]);
			vRefParamStd[i + tIndex] = ES.thickness[i];
		}
	}
	
	if (solve_tx_height){
		vRefParam[tx_heightIndex] = GR.tx_height;
		vRefParamStd[tx_heightIndex] = GS.tx_height;
	}

	if (solve_txrx_dx){
		vRefParam[txrx_dxIndex] = GR.txrx_dx;
		vRefParamStd[txrx_dxIndex] = GS.txrx_dx;
	}
	if (solve_txrx_dy){
		vRefParam[txrx_dyIndex] = GR.txrx_dy;
		vRefParamStd[txrx_dyIndex] = GS.txrx_dy;
	}
	if (solve_txrx_dz){
		vRefParam[txrx_dzIndex] = GR.txrx_dz;
		vRefParamStd[txrx_dzIndex] = GS.txrx_dz;
	}

	if (solve_rx_pitch){
		vRefParam[rx_pitchIndex] = GR.rx_pitch;
		vRefParamStd[rx_pitchIndex] = GS.rx_pitch;
	}

	if (solve_rx_roll){
		vRefParam[rx_rollIndex] = GR.rx_roll;
		vRefParamStd[rx_rollIndex] = GS.rx_roll;
	}

}
void cSBSInverter::initialise_Wd()
{
	Wd = MatrixDouble(ndata, ndata, 0.0);
	double s = 1.0 / (double)ndata;
	for (size_t i = 0; i < ndata; i++){
		Wd[i][i] = s / (vErr[i]*vErr[i]);
	}
	if(OO.Dump) writetofile(Wd,OO.DumpPath+"Wd.dat");
}
void cSBSInverter::initialise_Wc()
{	
	Wc = MatrixDouble(nparam, nparam, 0.0);	
	if (solve_conductivity == false)return;

	std::vector<double> t(nlayers);
	if (nlayers == 1){
		t[0] = 1;
	}
	else if (nlayers == 2){
		t[0] = ER.thickness[0];
		t[1] = ER.thickness[0];
	}
	else{
		for (size_t i = 0; i < (nlayers - 1); i++){
			t[i] = ER.thickness[i];
		}
		t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3])*t[nlayers - 2];
	}


	double tsum = 0.0;
	for (size_t i = 0; i < nlayers; i++)tsum += t[i];
	double tavg = tsum / (double)nlayers;
	
	double s = AlphaC / (double)(nlayers);
	for (size_t i = 0; i < nlayers; i++){
		size_t p = i + cIndex;		
		Wc[p][p] = s * (t[i]/tavg) / (vRefParamStd[p] * vRefParamStd[p]);				
	}			
}
void cSBSInverter::initialise_Wt()
{
	Wt = MatrixDouble(nparam, nparam, 0.0);
	if (solve_thickness == false)return;

	double s = AlphaT / (double)(nlayers - 1);
	for (size_t i = 0; i < nlayers - 1; i++){
		size_t p = i + tIndex;
		Wt[p][p] = s / (vRefParamStd[p] * vRefParamStd[p]);
	}	
	
}
void cSBSInverter::initialise_Wg()
{
	Wg = MatrixDouble(nparam, nparam, 0.0);
	if (ngeomparam <= 0)return;

	double s = AlphaG / (double)ngeomparam;
	for (size_t i = 0; i < ngeomparam; i++){
		size_t p = i + gIndex;
		Wg[p][p] = s / (vRefParamStd[p] * vRefParamStd[p]);
	}	
}
void cSBSInverter::initialise_L_Ws_1st_derivative()
{
	Ws = MatrixDouble(nparam, nparam, 0.0);
	if (AlphaS == 0 || nlayers < 3) return;
	if (solve_conductivity == false) return;

	std::vector<double> t(nlayers);
	for (size_t i = 0; i < (nlayers - 1); i++){
		t[i] = ER.thickness[i];
	}
	t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3])*t[nlayers - 2];


	double tsum = 0.0;
	for (size_t i = 0; i < nlayers; i++)tsum += t[i];
	double tavg = tsum / (double)nlayers;

	MatrixDouble L = MatrixDouble(nlayers - 1, nparam, 0.0);
	size_t neqn = 0;
	for (size_t li = 1; li < nlayers; li++){
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
void cSBSInverter::initialise_L_Ws_2nd_derivative()
{
	Ws = MatrixDouble(nparam, nparam, 0.0);
	if (AlphaS == 0 || nlayers < 3) return;
	if (solve_conductivity == false) return;

	std::vector<double> t(nlayers);
	if (nlayers == 1){
		t[0] = 1.0;
	}
	else{
		for (size_t i = 0; i < (nlayers - 1); i++){
			t[i] = ER.thickness[i];
		}
		t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3])*t[nlayers - 2];
	}
	
	double tsum = 0.0;
	for (size_t i = 0; i < nlayers; i++)tsum += t[i];
	double tavg = tsum / (double)nlayers;

	MatrixDouble L = MatrixDouble(nlayers - 2, nparam, 0.0);
	size_t neqn = 0;
	for (size_t li = 1; li < nlayers - 1; li++){
		size_t pindex = cIndex + li;
		double t1 = t[li - 1];
		double t2 = t[li];
		double t3 = t[li + 1];
		double d12 = (t1 + t2) / 2.0;
		double d23 = (t2 + t3) / 2.0;
		double s = sqrt(t2 / tavg);//sqrt because it gets squared in L'L		
		L[neqn][pindex - 1] =  s / d12;
		L[neqn][pindex]     = -s / d12 - s / d23;
		L[neqn][pindex + 1] =  s / d23;
		neqn++;
	}
	Ws  = transpose_mult(L, L);
	Ws *= (AlphaS / (double)(nlayers - 2));	
}
void cSBSInverter::initialise_Wr_Wm()
{
	Wr = MatrixDouble(nparam, nparam, 0.0);
	if (AlphaC > 0.0) Wr += Wc;
	if (AlphaT > 0.0) Wr += Wt;
	if (AlphaG > 0.0) Wr += Wg;
	
	Wm = Wr + Ws;

	if (OO.Dump){
		writetofile(Wc, OO.DumpPath + "Wc.dat");
		writetofile(Wt, OO.DumpPath + "Wt.dat");
		writetofile(Wg, OO.DumpPath + "Wg.dat");
		writetofile(Wr, OO.DumpPath + "Wr.dat");
		writetofile(Ws, OO.DumpPath + "Ws.dat");
		writetofile(Wm, OO.DumpPath + "Wm.dat");
	}
	
}
std::vector<double> cSBSInverter::parameterchange(const double lambda)
{	
	std::vector<double> x  = solve(lambda);
	std::vector<double> dm = x - vParam;

	if (solve_conductivity){
		for (size_t li = 0; li < nlayers; li++){
			size_t pindex = li + cIndex;
			const double maxcond = 50;
			const double mincond = 1e-6;
			if (vParam[pindex] + dm[pindex] > log10(maxcond)){
				//printf("upper limit li=%lu pindex=%lu dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] = log10(maxcond) - vParam[pindex];
			}
			else if (vParam[pindex] + dm[pindex] < log10(mincond)){
				//printf("lower limit li=%lu pindex=%lu dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] = log10(mincond) - vParam[pindex];
			}
		}
	}

	if (solve_thickness){
		for (size_t li = 0; li<nlayers - 1; li++){
			size_t pindex = li + tIndex;
			if (dm[pindex] > 0.5){
				//printf("li=%lu pindex=%lu dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] = 0.5;
			}
			else if (dm[pindex] < -0.5){
				//printf("li=%lu pindex=%lu dm=%lf\n",li,pindex,dm[pindex]);
				dm[pindex] = -0.5;
			}
		}
	}

	if (solve_tx_height){
		size_t pindex = tx_heightIndex;
		if (dm[pindex] > 0.5){
			//printf("li=%lu pindex=%lu dm=%lf\n",li,pindex,dm[pindex]);
			dm[pindex] = 0.5;
		}
		else if (dm[pindex] < -0.5){
			//printf("li=%lu pindex=%lu dm=%lf\n",li,pindex,dm[pindex]);
			dm[pindex] = -0.5;
		}

		if (vParam[pindex] + dm[pindex] > 1000){
			dm[pindex] = 1000 - vParam[pindex];
		}
		else if (vParam[pindex] + dm[pindex] < 10){
			dm[pindex] = 10 - vParam[pindex];
		}
	}
	return dm;
}
std::vector<double> cSBSInverter::solve(const double lambda)
{	
	// Phi = (d-g(m)+Jm) Wd (d-g(m)+Jm) + lambda ( (m-m0)' Wr (m-m0) + m' Ws m) )
	//Ax = b
	//A = [J'WdJ + lambda (Wr + Ws)]
	//x = m(n+1)
	//b = J'Wd(d - g(m) + Jm) + lambda*Wr*m0
	//dm = m(n+1) - m = x - m

	const std::vector<double>& m = vParam;
	const std::vector<double>& d = vObs;
	const std::vector<double>& g = vPred;
	const std::vector<double>& m0 = vRefParam;

	std::vector<double> b = JtWd*(d - g + J*m) + lambda*Wr*m0;
	MatrixDouble   A = JtWdJ + lambda*Wm;

	std::vector<double> x = pseudoinverse_od(A)*b;

	if (OO.Dump){
		writetofile(d, OO.DumpPath + "d.dat");
		writetofile(g, OO.DumpPath + "g.dat");
		writetofile(m, OO.DumpPath + "m.dat");
		writetofile(m0, OO.DumpPath + "m0.dat");
		writetofile(J, OO.DumpPath + "J.dat");
		writetofile(JtWdJ, OO.DumpPath + "JtWdJ.dat");
		writetofile(b, OO.DumpPath + "b.dat");
		writetofile(A, OO.DumpPath + "A.dat");
	}

	return x;
}
double cSBSInverter::phiData(const std::vector<double>& g)
{
	std::vector<double> v = vObs - g;
	double phid = mtDm(v, Wd);

	//this reports invalid models 
	if (phid < 0.0){
		phid = 1e9;
		warningmessage("Caught invalid PhiD\n");
	}

	return phid;
}
double cSBSInverter::phiModel(const std::vector<double>& p)
{
	double phic,phit,phig,phis;
	return phiModel(p, phic, phit, phig, phis);
}
double cSBSInverter::phiModel(const std::vector<double>& p, double& phic, double& phit, double& phig, double& phis)
{
	phic = phiC(p);
	phit = phiT(p);
	phig = phiG(p);
	phis = phiS(p);

	double v = phic + phit + phig + phis;	
	return v;
}
double cSBSInverter::phiC(const std::vector<double>& p)
{
	if (AlphaC == 0.0)return 0.0;
	if (solve_conductivity == false)return 0.0;
	std::vector<double> v = p - vRefParam;
	return mtDm(v, Wc);
}
double cSBSInverter::phiT(const std::vector<double>& p)
{
	if (AlphaT == 0.0)return 0.0;
	if (solve_thickness == false)return 0.0;
	std::vector<double> v = p - vRefParam;
	return mtDm(v, Wt);
}
double cSBSInverter::phiG(const std::vector<double>& p)
{
	if (AlphaG == 0.0)return 0.0;
	if (ngeomparam == 0)return 0.0;
	std::vector<double> v = p - vRefParam;
	return mtDm(v, Wg);
}
double cSBSInverter::phiS(const std::vector<double>& p)
{
	if (AlphaS == 0)return 0.0;
	else return mtAm(p, Ws);
}
cEarth1D cSBSInverter::get_earth(const std::vector<double>& parameters)
{
	cEarth1D e = ER;
	if (solve_conductivity){
		for (size_t li = 0; li < nlayers; li++){
			e.conductivity[li] = pow10(parameters[li + cIndex]);
		}
	}

	if (solve_thickness){
		for (size_t li = 0; li < nlayers - 1; li++){
			e.thickness[li] = pow10(parameters[li + tIndex]);
		}
	}
	return e;
}
cTDEmGeometry cSBSInverter::get_geometry(const std::vector<double>& parameters)
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
void cSBSInverter::set_predicted()
{
	for (size_t si = 0; si < nsystems; si++){
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
void cSBSInverter::forwardmodel(const std::vector<double>& parameters, std::vector<double>& predicted, bool computederivatives)
{
	cEarth1D      e = get_earth(parameters);
	cTDEmGeometry g = get_geometry(parameters);
	for (size_t si = 0; si < nsystems; si++){
		cTDEmSystemInfo& S = SV[si];
		cTDEmSystem& T = S.T;
		T.setconductivitythickness(e.conductivity, e.thickness);
		T.setgeometry(g);

		//Forwardmodel
		T.LEM.calculation_type = CT_FORWARDMODEL;
		T.LEM.derivative_layer = INT_MAX;
		T.setupcomputations();
		T.setprimaryfields();
		T.setsecondaryfields();

		std::vector<double> x, y, z;
		if (S.useTotal){			
			x = T.X + T.PrimaryX;
			y = T.Y + T.PrimaryY;
			z = T.Z + T.PrimaryZ;
		}
		else{
			x = T.X;
			y = T.Y;
			z = T.Z;
		}


		size_t nw = T.NumberOfWindows;
		for (size_t wi = 0; wi < nw; wi++){
			if (S.useX) predicted[wi + S.xIndex] = x[wi];
			if (S.useY) predicted[wi + S.yIndex] = y[wi];
			if (S.useZ) predicted[wi + S.zIndex] = z[wi];
		}

		if (computederivatives){
			//Receiver orientation derivatives go here while the forward model is still current
			if (solve_rx_pitch){
				size_t pindex = rx_pitchIndex;
				std::vector<double> dxbdp;
				std::vector<double> dzbdp;
				T.drx_pitch(x, z, g.rx_pitch, dxbdp, dzbdp);
				for (size_t i = 0; i < T.NumberOfWindows; i++){
					if (S.useX) J[i + S.xIndex][pindex] = dxbdp[i];
					if (S.useY) J[i + S.yIndex][pindex] = 0.0;
					if (S.useZ) J[i + S.zIndex][pindex] = dzbdp[i];
				}
			}

			if (solve_rx_roll){
				size_t pindex = rx_rollIndex;
				std::vector<double> dybdr;
				std::vector<double> dzbdr;
				T.drx_roll(y, z, g.rx_pitch, dybdr, dzbdr);
				for (size_t i = 0; i < T.NumberOfWindows; i++){
					if (S.useX) J[i + S.xIndex][pindex] = 0.0;
					if (S.useY) J[i + S.yIndex][pindex] = dybdr[i];
					if (S.useZ) J[i + S.zIndex][pindex] = dzbdr[i];
				}
			}

			if (solve_conductivity){
				for (size_t li = 0; li < nlayers; li++){
					size_t pindex = li + cIndex;
					T.LEM.calculation_type = CT_CONDUCTIVITYDERIVATIVE;
					T.LEM.derivative_layer = li;
					T.setprimaryfields();
					T.setsecondaryfields();
					//multiply by natural log(10) as parameters are in logbase10 units
					double c = log(10.0)*e.conductivity[li];
					if (S.useTotal){
						x = c*(T.X + T.PrimaryX); y = c*(T.Y + T.PrimaryY); z = c*(T.Z + T.PrimaryZ);
					}
					else{
						x = c*T.X;
						y = c*T.Y;
						z = c*T.Z;
					}
					for (size_t wi = 0; wi < nw; wi++){
						if (S.useX)J[wi + S.xIndex][pindex] = x[wi];
						if (S.useY)J[wi + S.yIndex][pindex] = y[wi];
						if (S.useZ)J[wi + S.zIndex][pindex] = z[wi];
					}
				}
			}
			if (solve_thickness){
				for (size_t li = 0; li < nlayers - 1; li++){
					size_t pindex = li + tIndex;
					T.LEM.calculation_type = CT_THICKNESSDERIVATIVE;
					T.LEM.derivative_layer = li;
					T.setprimaryfields();
					T.setsecondaryfields();
					//multiply by natural log(10) as parameters are in logbase10 units
					double c = log(10.0)*e.thickness[li];
					if (S.useTotal){
						x = c*(T.X + T.PrimaryX); y = c*(T.Y + T.PrimaryY); z = c*(T.Z + T.PrimaryZ);
					}
					else{
						x = c*T.X; y = c*T.Y; z = c*T.Z;
					}
					for (size_t wi = 0; wi < nw; wi++){
						if (S.useX)J[wi + S.xIndex][pindex] = x[wi];
						if (S.useY)J[wi + S.yIndex][pindex] = y[wi];
						if (S.useZ)J[wi + S.zIndex][pindex] = z[wi];
					}
				}
			}

			if (solve_tx_height){
				size_t pindex = tx_heightIndex;
				T.LEM.calculation_type = CT_HDERIVATIVE;
				T.LEM.derivative_layer = INT_MAX;
				T.setprimaryfields();
				T.setsecondaryfields();
				if (S.useTotal){
					x = (T.X + T.PrimaryX); y = (T.Y + T.PrimaryY); z = (T.Z + T.PrimaryZ);
				}
				else{
					x = T.X; y = T.Y; z = T.Z;
				}
				for (size_t wi = 0; wi < nw; wi++){
					if (S.useX)J[wi + S.xIndex][pindex] = x[wi];
					if (S.useY)J[wi + S.yIndex][pindex] = y[wi];
					if (S.useZ)J[wi + S.zIndex][pindex] = z[wi];
				}
			}

			if (solve_txrx_dx){
				size_t pindex = txrx_dxIndex;
				T.LEM.calculation_type = CT_XDERIVATIVE;
				T.LEM.derivative_layer = INT_MAX;
				T.setprimaryfields();
				T.setsecondaryfields();
				if (S.useTotal){
					x = (T.X + T.PrimaryX); y = (T.Y + T.PrimaryY); z = (T.Z + T.PrimaryZ);
				}
				else{
					x = T.X; y = T.Y; z = T.Z;
				}
				for (size_t wi = 0; wi < nw; wi++){
					if (S.useX)J[wi + S.xIndex][pindex] = x[wi];
					if (S.useY)J[wi + S.yIndex][pindex] = y[wi];
					if (S.useZ)J[wi + S.zIndex][pindex] = z[wi];
				}
			}

			if (solve_txrx_dy){
				size_t pindex = txrx_dyIndex;
				T.LEM.calculation_type = CT_YDERIVATIVE;
				T.LEM.derivative_layer = INT_MAX;
				T.setprimaryfields();
				T.setsecondaryfields();
				if (S.useTotal){
					x = (T.X + T.PrimaryX); y = (T.Y + T.PrimaryY); z = (T.Z + T.PrimaryZ);
				}
				else{
					x = T.X; y = T.Y; z = T.Z;
				}
				for (size_t wi = 0; wi < nw; wi++){
					if (S.useX)J[wi + S.xIndex][pindex] = x[wi];
					if (S.useY)J[wi + S.yIndex][pindex] = y[wi];
					if (S.useZ)J[wi + S.zIndex][pindex] = z[wi];
				}
			}

			if (solve_txrx_dz){
				size_t pindex = txrx_dzIndex;
				T.LEM.calculation_type = CT_ZDERIVATIVE;
				T.LEM.derivative_layer = INT_MAX;
				T.setprimaryfields();
				T.setsecondaryfields();
				if (S.useTotal){
					x = (T.X + T.PrimaryX); y = (T.Y + T.PrimaryY); z = (T.Z + T.PrimaryZ);
				}
				else{
					x = T.X; y = T.Y; z = T.Z;
				}
				for (size_t wi = 0; wi < nw; wi++){
					if (S.useX)J[wi + S.xIndex][pindex] = x[wi];
					if (S.useY)J[wi + S.yIndex][pindex] = y[wi];
					if (S.useZ)J[wi + S.zIndex][pindex] = z[wi];
				}
			}
		}		
		JtWd  = transpose_mult(J, Wd);
		JtWdJ = JtWd*J;


		//if(Dump){		
		//	writetofile(J,DumpPath+"J.dat");			
		//	writetofile(JtWd,DumpPath+"JtWd.dat");			
		//	writetofile(JtWdJ,DumpPath+"JtWdJ.dat");
		//}	
	}
}
std::vector<double> cSBSInverter::compute_parameter_sensitivity()
{
	std::vector<double> s(nparam, 0.0);	
	for (size_t pi = 0; pi < nparam; pi++){
		for (size_t di = 0; di < ndata; di++){			
			s[pi] += (fabs(J[di][pi]) * sqrt((double)ndata*Wd[di][di]));
		}
	}	

	if (OO.Dump){
		dumptofile(s, "layer_sensitivity.dat");
		writetofile(JtWdJ, OO.DumpPath + "JtWdJ.dat");
	}
	return s;
}
std::vector<double> cSBSInverter::compute_parameter_uncertainty()
{			
	MatrixDouble iCm(nparam, nparam, 0.0);
	for (size_t i = 0; i<nparam; i++) iCm[i][i] = 1.0/(vRefParamStd[i]*vRefParamStd[i]);
	MatrixDouble X = (double)ndata*JtWdJ + iCm;
	//MatrixDouble X = (double)ndata*JtWdJ + iCm + Ws;	
	//MatrixDouble X = (double)ndata*JtWdJ + Ws;
	//MatrixDouble X = (double)ndata*JtWdJ;
	MatrixDouble pinvX = pseudoinverse(X);
	std::vector<double> s(nparam);
	for (size_t i = 0; i < nparam; i++){
		s[i] = sqrt(pinvX[i][i]);
	}	
	return s;
}

void cSBSInverter::invert()
{
	parserecord();
	iterate();
}

void cSBSInverter::iterate()
{		
	double percentchange = 100.0;
	vParam = vRefParam;
	LastLambda = 1e8;
	LastIteration = 0;	
	LastPhiM   = phiModel(vParam, LastPhiC, LastPhiT, LastPhiG, LastPhiS);
	TerminationReason = "Has not terminated";
	
	size_t iteration = 0;	
	bool keepiterating = true;		
	do{		
		forwardmodel(vParam, vPred, true);
		LastPhiD = phiData(vPred);
		
		if (iteration >= MaxIterations){
			keepiterating = false;
			TerminationReason = "Too many iterations";			
		}
		else if (LastPhiD <= MinimumPhiD){
			keepiterating = false;
			TerminationReason = "Reached minimum";			
		}				
		else if (percentchange < MinimumImprovement){
			keepiterating = false;
			TerminationReason = "Small % improvement";
		}		
		else{			
			TargetPhiD = std::max(LastPhiD * 0.7, MinimumPhiD);
			sTrial t = targetsearch(LastLambda, TargetPhiD);

			std::vector<double> dm    = parameterchange(t.lambda);
			std::vector<double> mtemp = vParam + (t.stepfactor*dm);
			std::vector<double> gtemp(ndata);
			forwardmodel(mtemp, gtemp, false);
			double phidtemp = phiData(gtemp);
			percentchange = 100.0*(LastPhiD - phidtemp) / (LastPhiD);

			if (phidtemp < LastPhiD){	
				iteration++;
				vParam = mtemp;
				vPred  = gtemp;
				LastPhiD = phidtemp;
				LastLambda = t.lambda;
				LastIteration = iteration;
				LastPhiM = phiModel(vParam, LastPhiC, LastPhiT, LastPhiG, LastPhiS);								
				//printf("%d %lf %lf %lf\n", LastIteration, LastLambda, LastPhiD, LastPhiM);
			}	
		}		
	}
	while (keepiterating == true);

	EM = get_earth(vParam);
	GM = get_geometry(vParam);
	forwardmodel(vParam, vPred, false);	
	set_predicted();
	
	forwardmodel(vParam, vPred, true);
	ParameterSensitivity = compute_parameter_sensitivity();
	ParameterUncertainty = compute_parameter_uncertainty();

	if (OO.Dump){
		dumptofile(EM, "earth_inv.dat");
		dumptofile(GM, "geometry_inv.dat");
		dumptofile(vPred, "predicted.dat");
		FILE* fp = fileopen(OO.DumpPath + "iteration.dat", "w");
		fprintf(fp, "Iteration\t%lu\n", LastIteration);
		fprintf(fp, "TargetPhiD\t%lf\n", TargetPhiD);
		fprintf(fp, "PhiD\t%lf\n", LastPhiD);
		fprintf(fp, "Lambda\t%lf\n", LastLambda);
		fprintf(fp, "\n");
		fclose(fp);
	}
}
sTrial cSBSInverter::targetsearch(const double currentlambda, const double target)
{
	cTrialCache T;
	T.target = target;
	eBracketResult b = brackettarget(T, target, currentlambda);
	
	if (b == BR_BRACKETED){
		//bracketed target - find with Brents Method
		double newphid = DBL_MIN;
		double lambda = brentsmethod(T, target, newphid);		
		return T.findlambda(lambda);
	}
	else if (b == BR_MINBRACKETED){
		//bracketed minimum but above target - take smallest phid
		return T.minphidtrial();
	}
	else if (b == BR_ALLBELOW){
		//all below target	- take the largest lambda			
		return T.maxlambdatrial();
	}
	else if (b == BR_ALLABOVE){
		//all above target - take smallest phid
		return T.minphidtrial();
	}
	else{
		errormessage(fp_log,"targetsearch(): Unknown value %d returned from target brackettarget()\n", b);
	}
	return T.minphidtrial();
}
bool cSBSInverter::istargetbraketed(cTrialCache& T)
{
	double target = T.target;
	sorttrial_lambda(T);
	for (size_t i = T.trial.size() - 1; i >= 1; i--){
		if (T.trial[i].phid >= target && T.trial[i - 1].phid <= target){
			return true;
		}
		if (T.trial[i].phid <= target && T.trial[i - 1].phid >= target){
			return true;
		}
	}
	return false;
}
bool cSBSInverter::isminbraketed(cTrialCache& T)
{
	size_t index = T.minphidindex();
	if (index == 0 || index == T.trial.size() - 1){
		return false;
	}

	double fa = T.trial[index - 1].phid;
	double fb = T.trial[index].phid;
	double fc = T.trial[index + 1].phid;
	if ((fb < fa) && (fb < fc)){
		return true;
	}
	return false;
}
eBracketResult cSBSInverter::brackettarget(cTrialCache& T, const double target, const double currentlambda)
{
	double startx = log10(currentlambda);
	if (LastIteration == 0){
		std::vector<double> x;
		x.push_back(8); x.push_back(6);
		x.push_back(4); x.push_back(2);
		x.push_back(1); x.push_back(0);
		for (size_t k = 0; k < x.size(); k++){
			trialfunction(T, pow10(x[k]));
			bool tarbrak = istargetbraketed(T);
			if (tarbrak){
				return BR_BRACKETED;//target bracketed		
			}
		}

		double minv = DBL_MAX;
		for (size_t k = 0; k < T.trial.size(); k++){
			if (fabs(T.trial[k].phid - target) < minv){
				minv = fabs(T.trial[k].phid - target);
				startx = log10(T.trial[k].lambda);
			}
		}
	}
	else{
		trialfunction(T, pow10(startx));
	}

	std::vector<double> x;
	x.push_back(+1); x.push_back(-1);
	x.push_back(+2); x.push_back(-2);
	x.push_back(+3); x.push_back(-3);
	for (size_t k = 0; k < x.size(); k++){
		trialfunction(T, pow10(startx + x[k]));
		bool tarbrak = istargetbraketed(T);
		if (tarbrak){
			return BR_BRACKETED;//target bracketed		
		}
	}
	//printtrials(T);
	//prompttocontinue();

	if (T.maxphid() < target){
		return BR_ALLBELOW;//all below target	
	}

	bool minbrak = isminbraketed(T);
	if (minbrak)return BR_MINBRACKETED;//min bracketed											

	return BR_ALLABOVE;//all above target
}
double cSBSInverter::brentsmethod(cTrialCache& T, const double target, double& newphid)
{
	sorttrial_lambda(T);
	size_t index = T.trial.size() - 1;
	for (size_t i = T.trial.size() - 1; i >= 1; i--){
		double f1 = T.trial[i].phid - target;
		double f2 = T.trial[i - 1].phid - target;
		if (f1 * f2 <= 0.0){
			index = i;
			break;
		}
	}

	//Adapted from http://en.wikipedia.org/wiki/Brent's_method
	double xerrorTol = 0.01;
	double yerrorTol = target*0.1;//10% accuracy is good enough

	double a = log10(T.trial[index - 1].lambda);
	double b = log10(T.trial[index].lambda);
	double fa = T.trial[index - 1].phid - target;
	double fb = T.trial[index].phid - target;
	if (fa * fb >= 0.0){
		warningmessage("brentsmethod(): Target must be bracketed for cSBSInverter::brentsmethod()\n");
	}

	double c = 0;
	double d = DBL_MAX;
	double fc = 0;
	double s = 0;
	double fs = 0;

	// if f(a) f(b) >= 0 then error-exit
	if (fa * fb >= 0)
	{
		if (fa < fb){
			newphid = fa + target;
			return pow10(a);
		}
		else{
			newphid = fb + target;
			return pow10(b);
		}
	}

	// if |f(a)| < |f(b)| then swap (a,b) end if
	if (fabs(fa) < fabs(fb)){
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
		if (fabs(fs) < yerrorTol){
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
		if (i > 20){
			printtrials(T);
			warningmessage("Too many bisections in cSBSInverter::brentsmethod()\n");
			newphid = fb + target;
			return pow10(b);
		}
	}
	if (fb < fa){
		newphid = fb + target;
		return pow10(b);
	}
	else{
		newphid = fa + target;
		return pow10(a);
	}
}
void cSBSInverter::printtrials(cTrialCache T)
{
	sorttrial_lambda(T);
	printf("\n");
	printf("CurrentLambda = %lf CurrentPhid = %lf    Target = %lf\n", LastLambda, LastPhiD, T.target);
	printf("N    Stepfactor       Lambda          Phid\n");
	for (size_t i = 0; i<T.trial.size(); i++){
		printf("%2lu %12g %12g %12g\n", T.trial[i].order, T.trial[i].stepfactor, T.trial[i].lambda, T.trial[i].phid);
	}
	printf("\n");
}
double cSBSInverter::goldensearch(double a, double b, double c, double xtol, const double lambda, const std::vector<double>& m, const std::vector<double>& dm, std::vector<double>& g, cTrialCache& cache)
{
	//adapted from http://en.wikipedia.org/wiki/Golden_section_search	
	const double resphi = 2 - ((1 + sqrt(5.0)) / 2.0);
	double x;
	if (c - b > b - a){
		x = b + resphi * (c - b);
	}
	else{
		x = b - resphi * (b - a);
	}

	//if(fabs(c - a) < tau * (fabs(b) + fabs(x))){
	//  return (c + a) / 2; 
	//}
	if (fabs(c - a) < xtol){
		return (c + a) / 2;
	}

	double fx = cache.sfsearch(x);
	if (fx < 0){
		sTrial t;
		std::vector<double> p = m + x*dm;
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
	if (fb < 0){
		sTrial t;
		std::vector<double> p = m + b*dm;
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
	if (fx < fb){
		if (c - b > b - a){
			return goldensearch(b, x, c, xtol, lambda, m, dm, g, cache);
		}
		else{
			return goldensearch(a, x, b, xtol, lambda, m, dm, g, cache);
		}
	}
	else{
		if (c - b > b - a){
			return goldensearch(a, b, x, xtol, lambda, m, dm, g, cache);
		}
		else{
			return goldensearch(x, b, c, xtol, lambda, m, dm, g, cache);
		}
	}
}
double cSBSInverter::trialfunction(cTrialCache& T, const double triallambda)
{
	std::vector<double> dm(nparam);
	std::vector<double> p(nparam);
	std::vector<double> g(ndata);
	dm = parameterchange(triallambda);
	cTrialCache cache;
	sTrial t0;
	t0.phid = LastPhiD;
	t0.phim = LastPhiM;
	t0.stepfactor = 0.0;
	t0.lambda = triallambda;
	t0.order = cache.trial.size();
	cache.trial.push_back(t0);

	sTrial t1;
	p = vParam + dm;
	forwardmodel(p, g, false);
	t1.phid = phiData(g);
	t1.phim = phiModel(p);
	t1.stepfactor = 1.0;
	t1.lambda = triallambda;
	t1.order = cache.trial.size();
	cache.trial.push_back(t1);

	double pcdiff = 100 * (t1.phid - t0.phid) / t0.phid;
	if (pcdiff > 0.0 || pcdiff < -1.0){
		//ie dont do not do golden search
		//if only tiny improvement				
		double xtol = 0.1;
		double gsf = goldensearch(0.0, 0.38196601125010510, 1.0, xtol, triallambda, vParam, dm, g, cache);
		sTrial t3;
		p = vParam + gsf*dm;
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

	sTrial t = cache.trial[minindex];
	t.order = T.trial.size();
	T.trial.push_back(t);
	return t.phid;
}
void cSBSInverter::writeresult()
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
	writeresult_geometry(buf, OI, GI, "","Input ", false);
	writeresult_geometry(buf, OI, GM, "inverted_", "Inverted ", true);
			
	//Earth	
	OI.addfield("nlayers", 'I', 4, 0);
	OI.setcomment("Number of layers");
	buf += strprint("%4lu", nlayers);
	
	OI.addfield("conductivity", 'E', 15, 6, nlayers);	
	OI.setunits("S/m"); OI.setcomment("Layer conductivity");
	for (size_t i = 0; i < nlayers; i++){
		buf += strprint("%15.6le", EM.conductivity[i]);
	}
		
	double bottomlayerthickness = 100.0;
	if (solve_thickness==false && nlayers > 1){
		bottomlayerthickness = EM.thickness[nlayers - 2];
	}
	std::vector<double> thickness = EM.thickness;
	thickness.push_back(bottomlayerthickness);

	OI.addfield("thickness", 'F', 9, 2, nlayers);
	OI.setunits("m"); OI.setcomment("Layer thickness");
	for (size_t i = 0; i < nlayers; i++){
		buf += strprint("%9.2lf", thickness[i]);
	}
	
	if (OO.PositiveLayerBottomDepths){
		OI.addfield("depth_bottom", 'F', 9, 2, nlayers);
		OI.setunits("m"); OI.setcomment("Depth to bottom of layer");
		double tsum = 0.0;
		for (size_t i = 0; i < nlayers; i++){
			buf += strprint("%9.2lf", tsum);
			tsum += thickness[i];
		}
	}

	if (OO.NegativeLayerBottomDepths){
		OI.addfield("depth_bottom_negative", 'F', 9, 2, nlayers);
		OI.setunits("m"); OI.setcomment("Negative of depth to bottom of layer");
		double tsum = 0.0;
		for (size_t i = 0; i < nlayers; i++){
			tsum += thickness[i];
			buf += strprint("%9.2lf", -tsum);			
		}
	}

	if (OO.InterfaceElevations){
		OI.addfield("elevation_interfaces", 'F', 9, 2, nlayers+1);
		OI.setunits("m"); OI.setcomment("Elevation of interfaces");
		double etop = Location.groundelevation;
		for (size_t i = 0; i < nlayers; i++){
			buf += strprint("%9.2lf", etop);
			etop -= thickness[i];
		}
		buf += strprint("%9.2lf", etop);
	}

	if (OO.ParameterSensitivity){
		if (solve_conductivity){
			OI.addfield("conductivity_sensitivity", 'E', 15, 6, nlayers);			
			for (size_t i = 0; i < nlayers; i++){
				buf += strprint("%15.6le", ParameterSensitivity[cIndex+i]);
			}
		}
		if (solve_thickness){
			OI.addfield("thickness_sensitivity", 'E', 15, 6, nlayers-1);			
			for (size_t i = 0; i < nlayers-1; i++){
				buf += strprint("%15.6le", ParameterSensitivity[tIndex+i]);
			}
		}
		if (solve_tx_height){
			OI.addfield("tx_height_sensitivity", 'E', 15, 6);
			buf += strprint("%15.6le", ParameterSensitivity[tx_heightIndex]);
		}
		if (solve_txrx_dx){
			OI.addfield("txrx_dx_sensitivity", 'E', 15, 6);
			buf += strprint("%15.6le", ParameterSensitivity[txrx_dxIndex]);
		}
		if (solve_txrx_dz){
			OI.addfield("txrx_dz_sensitivity", 'E', 15, 6);
			buf += strprint("%15.6le", ParameterSensitivity[txrx_dzIndex]);
		}
		if (solve_rx_pitch){
			OI.addfield("rx_pitch_sensitivity", 'E', 15, 6);
			buf += strprint("%15.6le", ParameterSensitivity[rx_pitchIndex]);
		}
	}

	if (OO.ParameterUncertainty){
		if (solve_conductivity){
			OI.addfield("conductivity_uncertainty", 'E', 15, 6, nlayers);			
			OI.setunits("log10(S/m)");
			for (size_t i = 0; i < nlayers; i++){
				buf += strprint("%15.6le", ParameterUncertainty[cIndex+i]);
			}
		}
		if (solve_thickness){
			OI.addfield("thickness_uncertainty", 'E', 15, 6, nlayers - 1);
			OI.setunits("log10(m)");
			for (size_t i = 0; i < nlayers - 1; i++){
				buf += strprint("%15.6le", ParameterUncertainty[tIndex+i]);
			}
		}
		if (solve_tx_height){
			OI.addfield("tx_height_uncertainty", 'E', 15, 6);
			OI.setunits("m");
			buf += strprint("%15.6le", ParameterUncertainty[tx_heightIndex]);
		}
		if (solve_txrx_dx){			
			OI.addfield("txrx_dx_uncertainty", 'E', 15, 6);
			OI.setunits("m");
			buf += strprint("%15.6le", ParameterUncertainty[txrx_dxIndex]);
		}
		if (solve_txrx_dz){			
			OI.addfield("txrx_dz_uncertainty", 'E', 15, 6);
			OI.setunits("m");
			buf += strprint("%15.6le", ParameterUncertainty[txrx_dzIndex]);
		}
		if (solve_rx_pitch){			
			OI.addfield("rx_pitch_uncertainty", 'E', 15, 6);
			OI.setunits("degrees");
			buf += strprint("%15.6le", ParameterUncertainty[rx_pitchIndex]);
		}
	}

	//ObservedData
	if (OO.ObservedData) {
		for (size_t si = 0; si < nsystems; si++) {
			cTDEmSystemInfo& S = SV[si];							
			if (S.useX) writeresult_component(buf, OI, si, "X", "observed", "Observed", 'E', 15, 6, S.oPX, S.oSX, S.useTotal);
			if (S.useY) writeresult_component(buf, OI, si, "Y", "observed", "Observed", 'E', 15, 6, S.oPY, S.oSY, S.useTotal);
			if (S.useZ) writeresult_component(buf, OI, si, "Z", "observed", "Observed", 'E', 15, 6, S.oPZ, S.oSZ, S.useTotal);			
		}
	}

	//Noise Estimates
	if (OO.NoiseEstimates) {
		for (size_t si = 0; si < nsystems; si++) {
			cTDEmSystemInfo& S = SV[si];
			if (S.useX) writeresult_component(buf, OI, si, "X", "noise", "Estimated noise", 'E', 15, 6, 0.0, S.oEX, false);
			if (S.useY) writeresult_component(buf, OI, si, "Y", "noise", "Estimated noise", 'E', 15, 6, 0.0, S.oEY, false);
			if (S.useZ) writeresult_component(buf, OI, si, "Z", "noise", "Estimated noise", 'E', 15, 6, 0.0, S.oEZ, false);
		}
	}

	//PredictedData
	if (OO.PredictedData){
		for (size_t si = 0; si < nsystems; si++){
			cTDEmSystemInfo& S = SV[si];
			if (S.useX) writeresult_component(buf, OI, si, "X", "predicted", "Predicted", 'E', 15, 6, S.predicted.xcomponent.Primary, S.predicted.xcomponent.Secondary, S.useTotal);
			if (S.useY) writeresult_component(buf, OI, si, "Y", "predicted", "Predicted", 'E', 15, 6, S.predicted.ycomponent.Primary, S.predicted.ycomponent.Secondary, S.useTotal);
			if (S.useZ) writeresult_component(buf, OI, si, "Z", "predicted", "Predicted", 'E', 15, 6, S.predicted.zcomponent.Primary, S.predicted.zcomponent.Secondary, S.useTotal);
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
	if (Outputrecord == 1){		
		sFilePathParts fpp = getfilepathparts(OO.DataFile);
		
		std::string hdrfile = fpp.directory + fpp.prefix + ".hdr";
		OI.write_simple_header(hdrfile);

		std::string aseggdffile = fpp.directory + fpp.prefix + ".dfn";
		OI.write_aseggdf_header(aseggdffile);		
	}
	Outputrecord++;
};

void cSBSInverter::writeresult_geometry(std::string& buf, cOutputFileInfo& OI, const cTDEmGeometry& g, const std::string& fieldnameprefix, const std::string& commentprefix, const bool invertedfieldsonly)
{
	for (size_t i = 0; i < g.size(); i++){
		if (invertedfieldsonly && solvegeometryindex(i) == false)continue;
		OI.addfield(fieldnameprefix + g.fname(i), 'F', 9, 2);
		OI.setunits(g.units(i));
		OI.setcomment(commentprefix + g.description(i));
		buf += strprint("%9.2lf", g[i]);		
	}
}

void cSBSInverter::writeresult_component(std::string& buf, cOutputFileInfo& OI, const size_t& sysnum, const std::string& comp, const std::string& nameprefix, const std::string& commprefix, const char& form, const size_t& width, const size_t& decimals, const double& p, std::vector<double>& s, const bool& includeprimary)
{	
	std::string sysfield = nameprefix + strprint("_EMSystem_%d_", (int)sysnum+1);
	std::string syscomm  = commprefix + strprint(" EMSystem %d ", (int)sysnum+1);

	std::string fmt;
	if (form == 'F') fmt = strprint("%%%d.%dlf", width, decimals);
	else if (form == 'E') fmt = strprint("%%%d.%dle", width, decimals);
	else {
		rootmessage("Invalid output format %c\n", form);
		std::string e = strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
		throw(std::runtime_error(e));		
	}
	
	if (includeprimary){
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

bool cSBSInverter::solvegeometryindex(const size_t index)
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
		rootmessage("Geometry index %llu out of range\n", index);
		std::string e = strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
		throw(std::runtime_error(e));
		break;
	}
}
int cSBSInverter::intvalue(const FieldDefinition& cd)
{
	return cd.intvalue(DataFileFieldStrings);
}
double cSBSInverter::doublevalue(const FieldDefinition& cd)
{
	return cd.doublevalue(DataFileFieldStrings);
}
std::vector<int> cSBSInverter::intvector(const FieldDefinition& cd, const size_t& n)
{
	return cd.intvector(DataFileFieldStrings, n);
}
std::vector<double> cSBSInverter::doublevector(const FieldDefinition& cd, const size_t& n)
{
	return cd.doublevector(DataFileFieldStrings, n);
}
cTDEmGeometry cSBSInverter::readgeometry(const std::vector<FieldDefinition>& gfd)
{
	cTDEmGeometry g;
	g.tx_height = doublevalue(gfd[0]);
	g.tx_roll = doublevalue(gfd[1]);
	g.tx_pitch = doublevalue(gfd[2]);
	g.tx_yaw = doublevalue(gfd[3]);
	g.txrx_dx = doublevalue(gfd[4]);
	g.txrx_dy = doublevalue(gfd[5]);
	g.txrx_dz = doublevalue(gfd[6]);
	g.rx_roll = doublevalue(gfd[7]);
	g.rx_pitch = doublevalue(gfd[8]);
	g.rx_yaw = doublevalue(gfd[9]);
	return g;
}
void cSBSInverter::dumptofile(const std::vector<double>& v, std::string path)
{
	FILE* fp = fileopen(OO.DumpPath + path, "w");
	for (size_t i = 0; i < v.size(); i++){
		fprintf(fp, "%le\n", v[i]);
	}
	fclose(fp);
}
void cSBSInverter::dumptofile(const cEarth1D& e, std::string path)
{
	FILE* fp = fileopen(OO.DumpPath + path, "w");
	size_t nl = e.conductivity.size();
	for (size_t i = 0; i < nl; i++){
		if (i < e.thickness.size()){
			fprintf(fp, "%e\t%e\n", e.conductivity[i], e.thickness[i]);
		}
		else{
			fprintf(fp, "%e\tInf\n", e.conductivity[i]);
		}
	}
	fclose(fp);
}
void cSBSInverter::dumptofile(const cTDEmGeometry& g, std::string path)
{
	FILE* fp = fileopen(OO.DumpPath + path, "w");
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

////////////////////////////////////////////////////////////////////
int process(std::string controlfile, size_t Size, size_t Rank, bool usingopenmp)
{
	cSBSInverter I(Size, Rank);	
	
	//If OpenMP is being used set the thread lock while FFTW initialises	
	if (usingopenmp){
		#if defined _OPENMP
			omp_set_lock(&fftw_thread_lock);
		#endif	
	}
	
	I.initialise(controlfile);

	if (usingopenmp){
		#if defined _OPENMP			
			omp_unset_lock(&fftw_thread_lock);
		#endif	
	}
				
	size_t record = 0;
	while (I.readnextrecord()){
		record++;
		if ((record - 1) % (Size) != Rank)continue;
		
		bool nonnumeric = I.contains_non_numeric_characters(I.DataFileRecordString);
		if (nonnumeric){
			rootmessage(I.fp_log, "Skipping non-numeric record at line %lu of Input DataFile %s\n", I.DataFileRecord,I.DataFileName.c_str());
			rootmessage(I.fp_log, "\n%s\n\n", I.DataFileRecordString.c_str());
			continue;
		}

		if (I.OO.Dump){
			FILE* fp = fileopen(I.OO.DumpPath + "record.dat", "w");
			fprintf(fp, "Record\t%lu", I.DataFileRecord);
			fclose(fp);
		}

		double t1 = gettime();
		I.invert();

		double t2 = gettime();
		double etime = t2 - t1;
		I.writeresult();

		message(I.fp_log, "Rec %6lu  %3lu  %5lu  %10lf  ", I.DataFileRecord, I.Id.flightnumber, I.Id.linenumber, I.Id.fidnumber);
		message(I.fp_log, "Its=%3lu  PhiD=%6.2lf  time=%.1lfs %s\n", I.LastIteration, I.LastPhiD, etime, I.TerminationReason.c_str());
	}	
	//my_barrier();	
	rootmessage(I.fp_log, "Logfile closing at %s\n", timestamp().c_str());
	return 0;
};

int main(int argc, char** argv)
{		
	_GSTPUSH_
	int exitstatus;
	_GSTPOP_
	int mpisize = 1;
	int mpirank = 0;
	std::string mpipname = "Standalone";
#if defined _MPI_ENABLED			
	MPI_Init(&argc, &argv);
	mpirank  = cMpiEnv::world_rank();
	mpisize  = cMpiEnv::world_size();
	mpipname = cMpiEnv::processor_name();
	if (mpirank == 0)printf("MPI Started Processes=%d\tRank=%d\tProcessor name = %s\n", mpisize, mpirank, mpipname.c_str());
#endif

	if (mpirank == 0){
		printf("%s\n", commandlinestring(argc, argv).c_str());
		printf("%s\n", versionstring(VERSION, __TIME__, __DATE__).c_str());				
	}
	
	if (argc < 2){				
		printf("Usage: %s control_file_name [number_of_openmp_threads]\n", argv[0]);		
		exitstatus = EXIT_FAILURE;
	}
	else if (argc == 2){		
		bool usingopenmp = false;
		std::string controlfile = string(argv[1]);							
		process(controlfile, (size_t)mpisize, (size_t)mpirank, usingopenmp);
		exitstatus = EXIT_SUCCESS;
	}
	else if (argc == 3){		
		int openmpsize = atoi(argv[2]);
		if (mpisize > 1){
			if (mpirank == 0){
				printf("**Error: You may not use OpenMP with MPI\n");
				printf("**       Do not use [number_of_openmp_threads] when launched with mpiexec or mpirun\n");
			}
			exitstatus = EXIT_FAILURE;
		}
		else if (openmpsize < 1){
			printf("**Error: %d is a silly number of threads\n", openmpsize);
			exitstatus = EXIT_FAILURE;		
		}
		else{
			#if defined _OPENMP				
				bool usingopenmp = true;
				std::string controlfile = string(argv[1]);

				int openmpmaxthreads = omp_get_max_threads();
				if (openmpsize > openmpmaxthreads){
					printf("**Warning: The number of requested threads (%d) is more than the processors available (%d)\n", openmpsize, openmpmaxthreads);
				}

				printf("OpenMP threading Processes=%d\n", openmpsize);
				omp_init_lock(&fftw_thread_lock);
				#pragma omp parallel num_threads(openmpsize)
				{
					int openmprank = omp_get_thread_num();
					process(controlfile, (size_t)openmpsize, (size_t)openmprank, usingopenmp);
				}
				exitstatus = EXIT_SUCCESS;
			#else
				printf("**Error: This executable has not been compiled with OpenMP enabbled\n");
				printf("**       Compile with OpenMP or do not specify [number_of_openmp_threads]\n");
				exitstatus = EXIT_FAILURE;
			#endif
		}
	}
	else{		
		printf("Usage: %s control_file_name [number_of_openmp_threads]\n", argv[0]);
		printf("Too many command line arguments\n");		
		exitstatus = EXIT_FAILURE;
	}	

	#if defined _MPI_ENABLED
		if(mpirank==0)printf("Finalizing MPI\n");		
		MPI_Finalize();
	#endif

	return exitstatus;
}

