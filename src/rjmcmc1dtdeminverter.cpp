/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <cmath>
#include <cfloat>

#include "general_utils.h"
#include "file_utils.h"
#include "file_formats.h"
#include "vector_utils.h"
#include "rjmcmc1dtdeminverter.h"

rjmcmc1dTDEmInverter::rjmcmc1dTDEmInverter(const std::string& executable, const std::string& controlfile, size_t size, size_t rank)
{
	mSize = size;
	mRank = rank;
	initialise(executable, controlfile);
	initialise_sampler();
}
rjmcmc1dTDEmInverter::~rjmcmc1dTDEmInverter()
{
	fclose(fp_indata);
	fclose(fp_log);
}
void rjmcmc1dTDEmInverter::initialise(const std::string& executable, const std::string& controlfile)
{
	rootmessage("Loading control file %s\n", controlfile.c_str());
	Control = cBlock(controlfile);
	
	cBlock OB = Control.findblock("Output");
	std::string OutputDirectory = OB.getstringvalue("Directory");
	addtrailingseparator(OutputDirectory);
	if (exists(OutputDirectory) == false){
		rootmessage("Creating OutputDirectory: %s\n", OutputDirectory.c_str());
		makedirectory(OutputDirectory);
	}
	
	std::string s = OutputDirectory + OB.getstringvalue("LogFile");
	std::string rankstr = stringvalue(mRank, ".%04lu");
	LogFile = insert_after_filename(s, rankstr);

	rootmessage("Opening log file %s\n", LogFile.c_str());
	fp_log = fileopen(LogFile, "w");
	rootmessage(fp_log, "Logfile opened on %s\n", timestamp().c_str());
	rootmessage(fp_log, "Executing %s\n", executable.c_str());
	rootmessage(fp_log, "Control file %s\n", controlfile.c_str());
	rootmessage(fp_log, "Version %s Compiled at %s on %s\n", VERSION, __TIME__, __DATE__);
	rootmessage(fp_log, "Working directory %s\n", getcurrentdirectory().c_str());
	rootmessage(fp_log, "Processes=%lu\tRank=%lu\n", mSize, mRank);
	Control.write(fp_log);

	//Type of sampling Multichain or not
	cBlock SB = Control.findblock("Sampler");
	SaveMaps = SB.getboolvalue("SaveMaps");
	SaveMapsRate = SB.getintvalue("SaveMapsRate");
	if (SaveMapsRate == INT_MIN)SaveMapsRate = 1;

	SaveChains = SB.getboolvalue("SaveChains");
	SaveChainsRate = SB.getintvalue("SaveChainsRate");
	if (SaveChainsRate == INT_MIN)SaveChainsRate = 1;

	cBlock IB = Control.findblock("Input");
	HeaderLines = (size_t)IB.getintvalue("HeaderLines");
	SubSample = (size_t)IB.getintvalue("Subsample");
	FirstRecord = (size_t)IB.getintvalue("FirstRecord");
	LastRecord = (size_t)IB.getintvalue("LastRecord");
	if (FirstRecord <= 0)FirstRecord = 1;
	if (LastRecord <= 0)LastRecord = UINT_MAX;
	CurrentRecord = 0;

	InputDataFile = IB.getstringvalue("DataFile");
	fixseparator(InputDataFile);
	fp_indata = fileopen(InputDataFile, "r");

	s = OutputDirectory + OB.getstringvalue("DataFile");
	fixseparator(s);
	OutputDataFile = insert_after_filename(s, rankstr);
	
	if (SaveMaps){
		MapsDirectory = OB.getstringvalue("MapsDirectory");
		addtrailingseparator(MapsDirectory);
		if (exists(MapsDirectory) == false){
			rootmessage(fp_log, "Creating MapsDirectory: %s\n", MapsDirectory.c_str());
			makedirectory(MapsDirectory);
		}
	}

	if (SaveChains){
		ChainsDirectory = OB.getstringvalue("ChainsDirectory");
		addtrailingseparator(ChainsDirectory);
		if (exists(ChainsDirectory) == false){
			rootmessage(fp_log, "Creating ChainsDirectory: %s\n", ChainsDirectory.c_str());
			makedirectory(ChainsDirectory);
		}		
	}

	//Load stm file		
	initialise_systems();
	getcolumnnumbers();
}
void rjmcmc1dTDEmInverter::initialise_systems()
{
	nsystems = (size_t)Control.getintvalue("NumberOfSystems");
	SV.resize(nsystems);
	ndata = 0;
	for (size_t i = 0; i < nsystems; i++){
		sTDEmSystemInfo& S = SV[i];
		std::string str = strprint("EMSystem%lu", i + 1);
		cBlock b = Control.findblock(str);

		cTDEmSystem& T = S.T;
		std::string stmfile = b.getstringvalue("SystemFile");
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
		S.estimateNoiseFromModel = b.getboolvalue("EstimateNoiseFromModel");
		S.ncomps = 0;
		if (S.useX){
			S.ncomps++;
			ndata += S.nwindows;
			if (S.estimateNoiseFromModel){
				S.x_multiplicativenoise = b.getdoublevalue("XMultiplicativeNoise");
				S.x_additivenoise = b.getdoublevector("XAdditiveNoise");
			}
			S.oSX.resize(S.nwindows);
			S.oEX.resize(S.nwindows);
		}
		if (S.useY){
			S.ncomps++;
			ndata += S.nwindows;
			if (S.estimateNoiseFromModel){
				S.y_multiplicativenoise = b.getdoublevalue("YMultiplicativeNoise");
				S.y_additivenoise = b.getdoublevector("YAdditiveNoise");
			}
			S.oSY.resize(S.nwindows);
			S.oEY.resize(S.nwindows);
		}
		if (S.useZ){
			S.ncomps++;
			ndata += S.nwindows;
			if (S.estimateNoiseFromModel){
				S.z_multiplicativenoise = b.getdoublevalue("ZMultiplicativeNoise");
				S.z_additivenoise = b.getdoublevector("ZAdditiveNoise");
			}
			S.oSZ.resize(S.nwindows);
			S.oEZ.resize(S.nwindows);
		}
		S.nchans = S.nwindows*S.ncomps;
	}
	obs.resize(ndata);
	err.resize(ndata);	
}
void rjmcmc1dTDEmInverter::getcolumnnumbers()
{
	cBlock b = Control.findblock("Input.Columns");

	fd_surveynumber.set(b, "SurveyNumber");
	fd_datenumber.set(b, "DateNumber");
	fd_flightnumber.set(b, "FlightNumber");
	fd_linenumber.set(b, "LineNumber");
	fd_fidnumber.set(b, "FidNumber");
	fd_xord.set(b, "Easting");
	fd_yord.set(b, "Northing");
	fd_elevation.set(b, "GroundElevation");
	fd_altimeter.set(b, "Altimeter");

	fd_geometry[0].set(b, "TX_Height");
	fd_geometry[1].set(b, "TX_Roll");
	fd_geometry[2].set(b, "TX_Pitch");
	fd_geometry[3].set(b, "TX_Yaw");
	fd_geometry[4].set(b, "TXRX_DX");
	fd_geometry[5].set(b, "TXRX_DY");
	fd_geometry[6].set(b, "TXRX_DZ");
	fd_geometry[7].set(b, "RX_Roll");
	fd_geometry[8].set(b, "RX_Pitch");
	fd_geometry[9].set(b, "RX_Yaw");

	for (size_t i = 0; i < nsystems; i++){
		sTDEmSystemInfo& S = SV[i];
		std::string str = strprint("EMSystem%lu", i + 1);

		cBlock c = Control.findblock(str);
		S.fd_oPX.set(c, "XComponentPrimary");
		S.fd_oPY.set(c, "YComponentPrimary");
		S.fd_oPZ.set(c, "ZComponentPrimary");

		S.fd_oSX.set(c, "XComponentSecondary");
		S.fd_oSY.set(c, "YComponentSecondary");
		S.fd_oSZ.set(c, "ZComponentSecondary");

		S.fd_oEX.set(c, "XComponentNoise");
		S.fd_oEY.set(c, "YComponentNoise");
		S.fd_oEZ.set(c, "ZComponentNoise");
	}
}
void rjmcmc1dTDEmInverter::initialise_sampler()
{
	cBlock b = Control.findblock("Sampler");

	description = "";
	reportstats = true;

	nl_min = (size_t)b.getintvalue("NLayersMin");
	nl_max = (size_t)b.getintvalue("NLayersMax");

	pmax = b.getdoublevalue("DepthMax");
	std::string str = b.getstringvalue("DepthScaling");
	if (strcmp("LOG10", str.c_str()) == 0){
		param_position = LOG10;
		pmax = log10(pmax);
	}
	else{
		param_position = LINEAR;
	}
	size_t ndepth = (size_t)b.getintvalue("NDepthCells");

	vmin = b.getdoublevalue("ConductivityMin");
	vmax = b.getdoublevalue("ConductivityMax");
	str = b.getstringvalue("ConductivityScaling");
	if (strcmp("LOG10", str.c_str()) == 0){
		param_value = LOG10;
		vmin = log10(vmin);
		vmax = log10(vmax);
	}
	else{
		param_value = LINEAR;
	}
	size_t ncond = (size_t)b.getintvalue("NConductivityCells");
	pmap.initialise(nl_min, nl_max, pmax, ndepth, vmin, vmax, ncond);

	nchains = (size_t)b.getintvalue("NChains");
	nsamples = (size_t)b.getintvalue("NSamples");
	nburnin = (size_t)b.getintvalue("NBurnIn");
	thinrate = (size_t)b.getintvalue("ThinRate");
	//sd_bd_valuechange = b.getdoublevalue("STDValueChangeBirthDeath");
	//sd_valuechange = b.getdoublevalue("STDValueChange");
	//sd_move = b.getdoublevalue("STDMove");


	for (size_t i = 0; i<NTYPE_UNKNOWN; i++){
		std::string str = strprint("Nuisance%lu", i + 1);
		cBlock c = b.findblock(str);
		if (c.Entries.size() == 0)break;

		rjMcMCNuisance n;
		std::vector<std::string> ninit;
		std::string typestr = c.getstringvalue("Type");
		n.settype(typestr);

		std::string sinit = c.getstringvalue("Initial");
		if (strcasecmp(sinit, "DataFile") == 0){
			n.value = DBL_MIN;
		}
		else if (strncasecmp(sinit, "Column", 6) == 0){
			n.value = DBL_MIN;
		}
		else if (strncasecmp(sinit, "-Column", 7) == 0){
			n.value = DBL_MIN;
		}
		else{
			n.value = c.getdoublevalue("Initial");
		}

		n.min = c.getdoublevalue("Min");
		n.max = c.getdoublevalue("Max");
		if (n.min >= n.max){
			printf("Error: Min >= Max for nuisance");
			exit(1);
		}
		if (n.max == DBL_MIN){
			n.max = DBL_MAX;
		}
		n.sd_valuechange = c.getdoublevalue("STDValueChange");
		ntemplate.push_back(n);
		ninitial.push_back(sinit);
	}

}
bool rjmcmc1dTDEmInverter::readnextrecord()
{
	if (filegetline(fp_indata, CurrentRecordString) == false){
		return false;
	}
	else{
		CurrentRecord++;
		if (CurrentRecord > LastRecord + HeaderLines){
			return false;
		}
		else{
			return true;
		}
	}
}
bool rjmcmc1dTDEmInverter::readnextrecord_thisprocess()
{
	//Skip to next record for this process	
	while (readnextrecord()){
		int  r = ((int)CurrentRecord - (int)HeaderLines - (int)FirstRecord);
		if (r < 0)continue;
		size_t n = (size_t)r%SubSample;
		size_t b = (size_t)floor((double)r / (double)SubSample);
		size_t p = b%mSize;
		if (n == 0 && p == mRank)return true;
	}
	return false;
}
void rjmcmc1dTDEmInverter::parsecurrentrecord()
{
	CurrentRecordFields = fieldparsestring(CurrentRecordString.c_str(), " ,\t\r\n");
	std::vector<std::string>& f = CurrentRecordFields;

	surveynumber = (size_t)fd_surveynumber.intvalue(f);
	datenumber = (size_t)fd_datenumber.intvalue(f);
	flightnumber = (size_t)fd_flightnumber.intvalue(f);
	linenumber = (size_t)fd_linenumber.intvalue(f);
	fidnumber = fd_fidnumber.doublevalue(f);
	xord = fd_xord.doublevalue(f);
	yord = fd_yord.doublevalue(f);
	elevation = fd_elevation.doublevalue(f);
	altimeter = fd_altimeter.doublevalue(f);

	IG.tx_height = fd_geometry[0].doublevalue(f);
	IG.tx_roll = fd_geometry[1].doublevalue(f);
	IG.tx_pitch = fd_geometry[2].doublevalue(f);
	IG.tx_yaw = fd_geometry[3].doublevalue(f);
	IG.txrx_dx = fd_geometry[4].doublevalue(f);
	IG.txrx_dy = fd_geometry[5].doublevalue(f);
	IG.txrx_dz = fd_geometry[6].doublevalue(f);
	IG.rx_roll = fd_geometry[7].doublevalue(f);
	IG.rx_pitch = fd_geometry[8].doublevalue(f);
	IG.rx_yaw = fd_geometry[9].doublevalue(f);

	for (size_t i = 0; i < nsystems; i++){
		sTDEmSystemInfo& S = SV[i];

		if (S.useX){
			S.oPX = S.fd_oPX.doublevalue(f);			
			S.oSX = S.fd_oSX.doublevector(f, S.nwindows);			
			if (S.estimateNoiseFromModel){
				for (size_t wi = 0; wi < S.nwindows; wi++){					
					double an = S.x_additivenoise[wi];
					double mn = 0.01 * S.x_multiplicativenoise * S.oSX[wi];
					S.oEX[wi] = sqrt(an*an + mn*mn);
				}				
			}
			else{
				S.oEX = S.fd_oEX.doublevector(f, S.nwindows);				
			}
		}

		if (S.useY){			
			S.oPY = S.fd_oPY.doublevalue(f);						
			S.oSY = S.fd_oSY.doublevector(f, S.nwindows);			
			if (S.estimateNoiseFromModel){
				for (size_t wi = 0; wi < S.nwindows; wi++){					
					double an = S.y_additivenoise[wi];
					double mn = 0.01 * S.y_multiplicativenoise * S.oSY[wi];
					S.oEY[wi] = sqrt(an*an + mn*mn);				
				}
			}
			else{				
				S.oEY = S.fd_oEY.doublevector(f, S.nwindows);				
			}
		}

		if (S.useZ){			
			S.oPZ = S.fd_oPZ.doublevalue(f);			
			S.oSZ = S.fd_oSZ.doublevector(f, S.nwindows);
			if (S.estimateNoiseFromModel){
				for (size_t wi = 0; wi < S.nwindows; wi++){								
					double an = S.z_additivenoise[wi];
					double mn = 0.01 * S.z_multiplicativenoise * S.oSZ[wi];
					S.oEZ[wi] = sqrt(an*an + mn*mn);
				}
			}
			else{				
				S.oEZ = S.fd_oEZ.doublevector(f, S.nwindows);
			}
		}		
	}
	set_data();
	set_nuisance();
}
void rjmcmc1dTDEmInverter::set_data()
{
	size_t di = 0;
	for (size_t i = 0; i < nsystems; i++){
		sTDEmSystemInfo& S = SV[i];
		cTDEmSystem& T = S.T;

		if (S.reconstructPrimary){
			T.setgeometry(IG);
			T.Earth.calculation_type = CT_FORWARDMODEL;
			T.Earth.derivative_layer = INT_MAX;
			T.setprimaryfields();

			S.oPX = T.PrimaryX;
			S.oPY = T.PrimaryY;
			S.oPZ = T.PrimaryZ;
		}

		if (S.useX){
			for (size_t wi = 0; wi < S.nwindows; wi++){
				err[di] = S.oEX[wi];
				obs[di] = S.oSX[wi];
				if (S.useTotal)obs[di] += S.oPX;
				di++;
			}
		}
		if (S.useY){
			for (size_t wi = 0; wi < S.nwindows; wi++){
				err[di] = S.oEY[wi];
				obs[di] = S.oSY[wi];
				if (S.useTotal)obs[di] += S.oPY;
				di++;
			}
		}
		if (S.useZ){
			for (size_t wi = 0; wi < S.nwindows; wi++){
				err[di] = S.oEZ[wi];
				obs[di] = S.oSZ[wi];
				if (S.useTotal)obs[di] += S.oPZ;
				di++;
			}
		}
	}	
}
void rjmcmc1dTDEmInverter::set_nuisance()
{
	nuisance_init = ntemplate;

	for (size_t i = 0; i < ntemplate.size(); i++){

		rjMcMCNuisance& n = nuisance_init[i];

		std::string s = ninitial[i];

		if (strcasecmp(s, "DataFile") == 0){
			n.value = DBL_MIN;
		}

		if (strncasecmp(s, "Column", 6) == 0){
			size_t col;
			sscanf(ninitial[i].c_str(), "Column %lu", &col);
			col = col - 1;//referenced from 1 not 0
			double v;
			sscanf(CurrentRecordFields[col].c_str(), "%lf", &v);
			n.value = v;
		}
		else if (strncasecmp(s, "-Column", 7) == 0){
			size_t col;
			sscanf(ninitial[i].c_str(), "-Column %lu", &col);
			col = col - 1;//referenced from 1 not 0
			double v;
			sscanf(CurrentRecordFields[col].c_str(), "%lf", &v);
			n.value = -v;
		}
		else{
			//Stay at initial constant
			n.value = n.value;
		}


		if (n.value == DBL_MIN){
			switch (n.type){
			case TX_HEIGHT:
				n.value = IG.tx_height; break;
			case TX_ROLL:
				n.value = IG.tx_roll; break;
			case TX_PITCH:
				n.value = IG.tx_pitch; break;
			case TX_YAW:
				n.value = IG.tx_yaw;	break;
			case TXRX_DX:
				n.value = IG.txrx_dx; break;
			case TXRX_DY:
				n.value = IG.txrx_dy; break;
			case TXRX_DZ:
				n.value = IG.txrx_dz; break;
			case RX_ROLL:
				n.value = IG.rx_roll; break;
			case RX_PITCH:
				n.value = IG.rx_pitch; break;
			case RX_YAW:
				n.value = IG.rx_yaw; break;
			case TXRX_DISTANCE:
				n.value = sqrt(IG.txrx_dx * IG.txrx_dx + IG.txrx_dz * IG.txrx_dz);
				break;
			case TXRX_ANGLE:
				n.value = R2D*atan2(IG.txrx_dz, IG.txrx_dx);
				break;
			default:break;
			}
		}
	}
}
void rjmcmc1dTDEmInverter::sample()
{
	std::string desc = strprint("%lu %lu %lu %lf %lf %lf %lu", CurrentRecord, flightnumber, linenumber, fidnumber, xord, yord, ndata);

	rjMcMC1DSampler::reset();
	description = desc;
	rjMcMC1DSampler::sample();
	
	std::string dstr = results_string();

	sFilePathParts fpp = getfilepathparts(OutputDataFile);
	std::string hdrfile     = fpp.directory + fpp.prefix + ".hdr";
	std::string aseggdffile = fpp.directory + fpp.prefix + ".dfn";

	//Output header file
	if (exists(hdrfile) == false){				
		OI.write_simple_header(hdrfile);
		OI.write_aseggdf_header(aseggdffile);		
	}

	//Output data record	
	FILE* fp = fileopen(OutputDataFile, "a");
	fprintf(fp, dstr.c_str());
	fclose(fp);

	writemapstofile();
	writechainstofile();
}
std::string rjmcmc1dTDEmInverter::results_string()
{
	size_t ndepthcells = pmap.npbins();
	std::vector<cHistogramStats<double>> hs = pmap.hstats();
	std::vector<double> hlike = pmap.modelmap(mHighestLikelihood);
	std::vector<double> lmfit = pmap.modelmap(mLowestMisfit);
	if (param_value == LOG10){
		for (size_t j = 0; j < ndepthcells; j++){
			hs[j].p10 = pow10(hs[j].p10);
			hs[j].p50 = pow10(hs[j].p50);
			hs[j].p90 = pow10(hs[j].p90);
			hs[j].mode = pow10(hs[j].mode);
			hs[j].mean = pow10(hs[j].mean);
			hlike[j] = pow10(hlike[j]);
			lmfit[j] = pow10(lmfit[j]);
		}
	}
	double misfit_lowest = mLowestMisfit.misfit() / double(ndata);

	size_t n = 0;
	double sum = 0.0;
	for (size_t ci = 0; ci < nchains; ci++){
		for (size_t mi = 0; mi < mChainInfo[ci].modelchain.size(); mi++){
			sum += mChainInfo[ci].modelchain[mi].misfit();
			n++;
		}
	}
	double misfit_average = sum / double(n) / double(ndata);
		
	std::string buf;
	//Id		
	OI.addfield("uniqueid", 'I', 12, 0);
	OI.setcomment("Inversion sequence number");
	buf += strprint("%12lu", CurrentRecord);

	OI.addfield("survey", 'I', 12, 0);
	OI.setcomment("Survey number");
	buf += strprint("%12lu", surveynumber);

	OI.addfield("date", 'I', 12, 0);
	OI.setcomment("Date number");
	buf += strprint("%12lu", datenumber);

	OI.addfield("flight", 'I', 12, 0);
	OI.setcomment("Flight number, IntrepidFlightNumber");
	buf += strprint("%12lu", flightnumber);

	OI.addfield("line", 'I', 12, 0);
	OI.setcomment("Line number, IntrepidLineNumber");
	buf += strprint("%12lu", linenumber);

	OI.addfield("fiducial", 'F', 12, 2);
	OI.setcomment("Fiducial number, IntrepidFiducial");
	buf += strprint("%12.2lf", fidnumber);

	//Location
	OI.addfield("easting", 'F', 9, 1);
	OI.setunits("m"); OI.setcomment("IntrepidX");
	buf += strprint("%9.1lf", xord);

	OI.addfield("northing", 'F', 10, 1);
	OI.setunits("m"); OI.setcomment("IntrepidY");
	buf += strprint("%10.1lf", yord);

	OI.addfield("elevation", 'F', 10, 2);
	OI.setunits("m"); OI.setcomment("Ground elevation relative to sea-level");
	buf += strprint("%10.2lf", elevation);

	OI.addfield("altimeter", 'F', 8, 2);
	OI.setunits("m"); OI.setcomment("Height of altimeter above ground-level");
	buf += strprint("%8.2lf", altimeter);

	///////////////////
	
	OI.addfield("nchains", 'I', 9, 0);
	OI.setcomment("Number of chains");	
	buf += strprint("%9lu", nchains);

	OI.addfield("nsamples", 'I', 9, 0);
	OI.setcomment("Number of samples per chain");
	buf += strprint("%9lu", nsamples);

	OI.addfield("nburnin", 'I', 9, 0);
	OI.setcomment("Number of samples in burn-in");
	buf += strprint("%9lu", nburnin);

	OI.addfield("sampletime", 'F', 8, 2);
	OI.setunits("s"); OI.setcomment("Sampling wall time in seconds");
	buf += strprint("%8.2lf", samplingtime);

	OI.addfield("misfit_lowest", 'E', 10, 6);
	OI.setcomment("Lowest misfit on any chain");
	buf += strprint("%15.6le", misfit_lowest);
	
	OI.addfield("misfit_average", 'E', 10, 6);
	OI.setcomment("Average misfit over all chains");
	buf += strprint("%15.6le", misfit_average);
	
	OI.addfield("ndepthcells", 'I', 4, 0);
	OI.setcomment("Number of depth cells or bins in histograms");
	buf += strprint("%4lu", ndepthcells);
	
	OI.addfield("depthcelltop", 'F', 8, 2, ndepthcells);
	OI.setunits("m"); OI.setcomment("Depths of cell tops");
	for (size_t j = 0; j < ndepthcells; j++){
		double dcc = pmap.toppbin(j);
		if (param_position == LOG10){
			dcc = pow10(dcc);
		}
		buf += strprint("%8.2lf", dcc);
	}
	

	size_t width = 15;
	size_t decimals = 6;
	std::string fmt = "%15.6le";

	OI.addfield("conductivity_mean", 'E', width, decimals, ndepthcells);
	OI.setunits("S/m"); OI.setcomment("Conductivity of mean model");
	for (size_t j = 0; j < ndepthcells; j++){
		buf += strprint(fmt.c_str(), hs[j].mean);
	}

	OI.addfield("conductivity_mode", 'E', width, decimals, ndepthcells);
	OI.setunits("S/m"); OI.setcomment("Conductivity of mode model");
	for (size_t j = 0; j < ndepthcells; j++){
		buf += strprint(fmt.c_str(), hs[j].mode);
	}

	OI.addfield("conductivity_p50", 'E', width, decimals, ndepthcells);
	OI.setunits("S/m"); OI.setcomment("Conductivity at 50th percentile or median");	
	for (size_t j = 0; j < ndepthcells; j++){
		buf += strprint(fmt.c_str(), hs[j].p50);
	}
	
	OI.addfield("conductivity_p10", 'E', width, decimals, ndepthcells);
	OI.setunits("S/m"); OI.setcomment("Conductivity at 10th percentile");	
	for (size_t j = 0; j < ndepthcells; j++){
		buf += strprint(fmt.c_str(), hs[j].p10);
	}

	OI.addfield("conductivity_p90", 'E', width, decimals, ndepthcells);
	OI.setunits("S/m"); OI.setcomment("Conductivity at 90th percentile");
	for (size_t j = 0; j < ndepthcells; j++){
		buf += strprint(fmt.c_str(), hs[j].p90);
	}
		
	OI.addfield("conductivity_highestlikelihood", 'E', width, decimals, ndepthcells);
	OI.setunits("S/m"); OI.setcomment("Conductivity of highest likelihood model");
	for (size_t j = 0; j < ndepthcells; j++){
		buf += strprint(fmt.c_str(), hlike[j]);
	}

	OI.addfield("conductivity_lowestmisfit", 'E', width, decimals, ndepthcells);
	OI.setunits("S/m"); OI.setcomment("Conductivity of howest misfit model");
	for (size_t j = 0; j < ndepthcells; j++){
		buf += strprint(fmt.c_str(), lmfit[j]);
	}
	
	OI.addfield("changepoint", 'I', 11, 0, ndepthcells);
	OI.setunits("counts"); OI.setcomment("Changepoint histogram");
	for (size_t j = 0; j < ndepthcells; j++){
		buf += strprint("%11lu", pmap.changepoint()[j]);
	}

	size_t nn = nnuisances();
	for (size_t j = 0; j < nn; j++){
		std::string nstr = mChainInfo[0].modelchain[0].nuisances[j].typestring();
		cStats<double> s(nmap.nuisance[j]);

		std::string hs;
		hs = nstr + "_inv_mean";		
		OI.addfield(hs.c_str(), 'F', 8, 2);		
		buf += strprint("%8.2lf", s.mean);

		hs = nstr + "_inv_std";		
		OI.addfield(hs.c_str(), 'F', 8, 2);
		buf += strprint("%8.2lf", s.std);
	}
	buf += strprint("\n");
	OI.lockfields();
	return buf;
}

std::string rjmcmc1dTDEmInverter::prefixstring()
{
	std::string s = "seq";
	s += strprint(".%08lu", CurrentRecord);
	s += strprint(".%lu", linenumber);
	s += strprint(".%lf", fidnumber);
	return s;
}
void rjmcmc1dTDEmInverter::writemapstofile()
{
	if (SaveMaps == false)return;
	if ((CurrentRecord - HeaderLines - FirstRecord) / SubSample % SaveMapsRate != 0)return;

	std::string fileprefix = prefixstring();
	std::string fname = MapsDirectory + fileprefix + ".pmap";
	FILE* fp = fileopen(fname.c_str(), "wb");

	writeheader(fp);
	writemodel(fp, mHighestLikelihood);
	writemodel(fp, mLowestMisfit);
	nmap.writedata(fp);
	pmap.writedata(fp);
	for (size_t ci = 0; ci < nchains; ci++){
		mChainInfo[ci].writeconvergencedata(ci, fp);
	}

	fclose(fp);
}
void rjmcmc1dTDEmInverter::writechainstofile()
{
	if (SaveChains == false)return;
	if ((CurrentRecord - HeaderLines - FirstRecord) / SubSample % SaveChainsRate != 0)return;

	std::string fileprefix = prefixstring();
	std::string fname = ChainsDirectory + fileprefix + ".chain";
	FILE* fp = fileopen(fname.c_str(), "w");
	bwrite(fp, (uint32_t)nchains);
	for (size_t ci = 0; ci < nchains; ci++){
		mChainInfo[ci].writemodelchain_binary(ci, fp);
	}
	fclose(fp);
}

sTDEmGeometry  rjmcmc1dTDEmInverter::getgeometry(const rjMcMC1DModel& m)
{
	sTDEmGeometry  OG = IG;

	bool angledistance = false;
	double angle = 0.0;
	double distance = 0.0;

	for (size_t i = 0; i < m.nuisances.size(); i++){
		const rjMcMCNuisance& n = m.nuisances[i];
		switch (n.type){
		case TX_HEIGHT:
			OG.tx_height = n.value; break;
		case TX_ROLL:
			OG.tx_roll = n.value; break;
		case TX_PITCH:
			OG.tx_pitch = n.value; break;
		case TX_YAW:
			OG.tx_yaw = n.value; break;
		case TXRX_DX:
			OG.txrx_dx = n.value; break;
		case TXRX_DY:
			OG.txrx_dy = n.value; break;
		case TXRX_DZ:
			OG.txrx_dz = n.value; break;
		case RX_ROLL:
			OG.rx_roll = n.value; break;
		case RX_PITCH:
			OG.rx_pitch = n.value; break;
		case RX_YAW:
			OG.rx_yaw = n.value; break;
		case TXRX_DISTANCE:
			angledistance = true;
			distance = n.value; break;
		case TXRX_ANGLE:
			angledistance = true;
			angle = n.value; break;
		default:
			break;
		}
	}

	if (angledistance == true){
		OG.txrx_dx = distance*cos(D2R*angle);
		OG.txrx_dz = distance*sin(D2R*angle);
	}
	return OG;
}
std::vector<double> rjmcmc1dTDEmInverter::collect(const sTDEmSystemInfo& S, const cTDEmSystem& T)
{
	std::vector<double> v(S.nchans);
	std::vector<double> x, y, z;
	if (S.useTotal){
		if (S.useX) x = T.X + T.PrimaryX;
		if (S.useY) y = T.Y + T.PrimaryY;
		if (S.useZ) z = T.Z + T.PrimaryZ;
	}
	else{
		if (S.useX) x = T.X;
		if (S.useY) y = T.Y;
		if (S.useZ) z = T.Z;
	}

	size_t nx = x.size();
	size_t ny = y.size();
	for (size_t i = 0; i < x.size(); i++) v[i] = x[i];
	for (size_t i = 0; i < y.size(); i++) v[i + nx] = y[i];
	for (size_t i = 0; i < z.size(); i++) v[i + nx + ny] = z[i];

	return v;
}
std::vector<double> rjmcmc1dTDEmInverter::forwardmodel(const rjMcMC1DModel& m)
{
	size_t nl = m.nlayers();
	std::vector<double> c = m.getvalues();
	if (param_value == LOG10){
		for (size_t i = 0; i < nl; i++){
			c[i] = pow10(c[i]);
		}
	}

	std::vector<double> t = m.getthicknesses();
	sTDEmGeometry  G = getgeometry(m);
	std::vector<double> pred(ndata);

	size_t di = 0;
	for (size_t i = 0; i < nsystems; i++){
		sTDEmSystemInfo& S = SV[i];
		cTDEmSystem& T = S.T;
		T.setconductivitythickness(c, t);
		T.setgeometry(G);
		T.setupcomputations();
		T.Earth.calculation_type = CT_FORWARDMODEL;
		T.Earth.derivative_layer = INT_MAX;
		T.setprimaryfields();
		T.setsecondaryfields();
		std::vector<double> v = collect(S, T);
		for (size_t j = 0; j < v.size(); j++){
			pred[di] = v[j];
			di++;
		}
	}
	return pred;
}


