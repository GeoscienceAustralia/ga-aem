/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _rjmcmc1dtdeminverter_H
#define _rjmcmc1dtdeminverter_H

#include "gaaem_version.h"
#include "blocklanguage.h"
#include "fielddefinition.h"
#include "tdemsystem.h"
#include "file_formats.h"
#include "rjmcmc1d.h"

class cTDEmSystemInfo{

public:
	cTDEmSystem T;	
	size_t ncomps;
	size_t nwindows;
	size_t nchans;
	
	bool useX;
	bool useY;
	bool useZ;
	bool useTotal;
	bool reconstructPrimary;
	bool estimateNoiseFromModel;

	double x_multiplicativenoise;
	double y_multiplicativenoise;
	double z_multiplicativenoise;
	std::vector<double> x_additivenoise;
	std::vector<double> y_additivenoise;
	std::vector<double> z_additivenoise;
	double  oPX,oPY,oPZ;
	std::vector<double> oSX,oSY,oSZ;	
	std::vector<double> oEX,oEY,oEZ;	
	cFieldDefinition fd_oPX,fd_oPY,fd_oPZ;
	cFieldDefinition fd_oSX,fd_oSY,fd_oSZ;
	cFieldDefinition fd_oEX,fd_oEY,fd_oEZ;
};

class rjmcmc1dTDEmInverter : public rjMcMC1DSampler{
	
	public:
		
	size_t nsystems;
	std::vector<cTDEmSystemInfo> SV;		
	std::vector<rjMcMCNuisance> ntemplate;
	std::vector<std::string> ninitial;
	cTDEmGeometry  IG;	
	cOutputFileInfo OI;
		
	bool SaveMaps;
	int  SaveMapsRate;
	bool SaveChains;
	int  SaveChainsRate;
		
	cBlock Control;
	std::string LogFile;	
	double memoryusedatstart;

	std::string InputDataFile;	
	FILE*  fp_indata;		

	size_t Outputcolumn; //output column number
	std::string OutputDirectory;
	std::string OutputDataFile;	
	std::string MapsDirectory;
	std::string ChainsDirectory;
		
	size_t HeaderLines;		
	size_t SubSample;	
	size_t FirstRecord;
	size_t LastRecord;
	size_t CurrentRecord;	
	std::string CurrentRecordString;
	std::vector<std::string> CurrentRecordFields;	
								
	size_t surveynumber;
	size_t datenumber;
	size_t flightnumber;
	size_t linenumber;	
	double fidnumber;
	double timenumber;
	double xord;
	double yord;
	double elevation;	

	cFieldDefinition fd_surveynumber;
	cFieldDefinition fd_datenumber;	
	cFieldDefinition fd_flightnumber;	
	cFieldDefinition fd_linenumber;
	cFieldDefinition fd_fidnumber;	
	cFieldDefinition fd_xord;
	cFieldDefinition fd_yord;
	cFieldDefinition fd_elevation;	
	cFieldDefinition fd_geometry[10];

	rjmcmc1dTDEmInverter(const std::string& executable, const std::string& controlfile, size_t size, size_t rank)
	{
		mpiSize = size;
		mpiRank = rank;
		initialise(executable, controlfile);
		initialise_sampler();
	}
	
	~rjmcmc1dTDEmInverter()
	{
		fclose(fp_indata);
		glog.close();
	}

	void initialise(const std::string& executable, const std::string& controlfile)
	{
		glog.logmsg(0, "Loading control file %s\n", controlfile.c_str());
		Control = cBlock(controlfile);

		cBlock OB = Control.findblock("Output");
		std::string OutputDirectory = OB.getstringvalue("Directory");
		addtrailingseparator(OutputDirectory);
		if (exists(OutputDirectory) == false) {
			glog.logmsg(0, "Creating OutputDirectory: %s\n", OutputDirectory.c_str());
			makedirectory(OutputDirectory);
		}

		std::string s = OutputDirectory + OB.getstringvalue("LogFile");
		std::string rankstr = stringvalue(mpiRank, ".%04lu");
		LogFile = insert_after_filename(s, rankstr);

		glog.logmsg(0, "Opening log file %s\n", LogFile.c_str());
		glog.open(LogFile);
		glog.logmsg(0, "Logfile opened on %s\n", timestamp().c_str());
		glog.logmsg(0, "Executing %s\n", executable.c_str());
		glog.logmsg(0, "Control file %s\n", controlfile.c_str());
		glog.logmsg(0, "Version %s Compiled at %s on %s\n", GAAEM_VERSION, __TIME__, __DATE__);
		glog.logmsg(0, "Working directory %s\n", getcurrentdirectory().c_str());
		glog.logmsg(0, "Processes=%lu\tRank=%lu\n", mpiSize, mpiRank);
		glog.log(Control.get_as_string());

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

		if (SaveMaps) {
			MapsDirectory = OB.getstringvalue("MapsDirectory");
			addtrailingseparator(MapsDirectory);
			if (exists(MapsDirectory) == false) {
				glog.logmsg(0, "Creating MapsDirectory: %s\n", MapsDirectory.c_str());
				makedirectory(MapsDirectory);
			}
		}

		if (SaveChains) {
			//ChainsDirectory = OB.getstringvalue("ChainsDirectory");
			//addtrailingseparator(ChainsDirectory);
			//if (exists(ChainsDirectory) == false) {
			//	glog.logmsg(0, "Creating ChainsDirectory: %s\n", ChainsDirectory.c_str());
			//	makedirectory(ChainsDirectory);
			//}
		}

		//Load stm file		
		initialise_systems();
		getcolumnnumbers();
	}
	
	void initialise_systems()
	{
		nsystems = (size_t)Control.getintvalue("NumberOfSystems");
		SV.resize(nsystems);
		ndata = 0;
		for (size_t i = 0; i < nsystems; i++) {
			cTDEmSystemInfo& S = SV[i];
			std::string str = strprint("EMSystem%lu", i + 1);
			cBlock b = Control.findblock(str);

			cTDEmSystem& T = S.T;
			std::string stmfile = b.getstringvalue("SystemFile");
			glog.logmsg(0, "Reading system file %s\n", stmfile.c_str());
			T.readsystemdescriptorfile(stmfile);

			glog.log("==============System file %s\n", stmfile.c_str());
			glog.log(T.STM.get_as_string());
			glog.log("==========================================================================\n");

			S.nwindows = T.NumberOfWindows;
			S.useX = b.getboolvalue("UseXComponent");
			S.useY = b.getboolvalue("UseYComponent");
			S.useZ = b.getboolvalue("UseZComponent");
			S.useTotal = b.getboolvalue("InvertTotalField");
			S.reconstructPrimary = b.getboolvalue("ReconstructPrimaryFieldFromInputGeometry");
			S.estimateNoiseFromModel = b.getboolvalue("EstimateNoiseFromModel");
			S.ncomps = 0;
			if (S.useX) {
				S.ncomps++;
				ndata += S.nwindows;
				if (S.estimateNoiseFromModel) {
					S.x_multiplicativenoise = b.getdoublevalue("XMultiplicativeNoise");
					S.x_additivenoise = b.getdoublevector("XAdditiveNoise");
				}
				S.oSX.resize(S.nwindows);
				S.oEX.resize(S.nwindows);
			}
			if (S.useY) {
				S.ncomps++;
				ndata += S.nwindows;
				if (S.estimateNoiseFromModel) {
					S.y_multiplicativenoise = b.getdoublevalue("YMultiplicativeNoise");
					S.y_additivenoise = b.getdoublevector("YAdditiveNoise");
				}
				S.oSY.resize(S.nwindows);
				S.oEY.resize(S.nwindows);
			}
			if (S.useZ) {
				S.ncomps++;
				ndata += S.nwindows;
				if (S.estimateNoiseFromModel) {
					S.z_multiplicativenoise = b.getdoublevalue("ZMultiplicativeNoise");
					S.z_additivenoise = b.getdoublevector("ZAdditiveNoise");
				}
				S.oSZ.resize(S.nwindows);
				S.oEZ.resize(S.nwindows);
			}
			S.nchans = S.nwindows * S.ncomps;
		}
		obs.resize(ndata);
		err.resize(ndata);
	}
	
	void getcolumnnumbers()
	{
		cBlock b = Control.findblock("Input.Columns");

		fd_surveynumber.initialise(b, "SurveyNumber");
		fd_datenumber.initialise(b, "DateNumber");
		fd_flightnumber.initialise(b, "FlightNumber");
		fd_linenumber.initialise(b, "LineNumber");
		fd_fidnumber.initialise(b, "FidNumber");
		fd_xord.initialise(b, "Easting");
		fd_yord.initialise(b, "Northing");
		fd_elevation.initialise(b, "GroundElevation");

		cTDEmGeometry g;
		for (size_t gi = 0; gi < g.size(); gi++) {
			fd_geometry[gi].initialise(b, g.fname(gi));
		}

		for (size_t i = 0; i < nsystems; i++) {
			cTDEmSystemInfo& S = SV[i];
			std::string str = strprint("EMSystem%zu", i + 1);

			cBlock c = Control.findblock(str);
			S.fd_oPX.initialise(c, "XComponentPrimary");
			S.fd_oPY.initialise(c, "YComponentPrimary");
			S.fd_oPZ.initialise(c, "ZComponentPrimary");

			S.fd_oSX.initialise(c, "XComponentSecondary");
			S.fd_oSY.initialise(c, "YComponentSecondary");
			S.fd_oSZ.initialise(c, "ZComponentSecondary");

			S.fd_oEX.initialise(c, "XComponentNoise");
			S.fd_oEY.initialise(c, "YComponentNoise");
			S.fd_oEZ.initialise(c, "ZComponentNoise");
		}
	}
	
	void initialise_sampler()
	{
		cBlock b = Control.findblock("Sampler");
		nl_min = (size_t)b.getintvalue("NLayersMin");
		nl_max = (size_t)b.getintvalue("NLayersMax");

		pmax = b.getdoublevalue("DepthMax");
		std::string str = b.getstringvalue("DepthScaling");
		param_position = cParameterization(str);
		if (param_position.islog10()) {
			pmax = std::log10(pmax);
		}
		size_t ndepth = (size_t)b.getintvalue("NDepthCells");

		vmin = b.getdoublevalue("ConductivityMin");
		vmax = b.getdoublevalue("ConductivityMax");
		str = b.getstringvalue("ConductivityScaling");
		param_value = cParameterization(str);
		if (param_value.islog10()) {
			vmin = std::log10(vmin);
			vmax = std::log10(vmax);
		}

		size_t ncond = (size_t)b.getintvalue("NConductivityCells");
		pmap.initialise(nl_min, nl_max, pmax, ndepth, vmin, vmax, ncond);

		size_t nc = b.getsizetvalue("NChains");
		chains.resize(nc);
		nsamples = b.getsizetvalue("NSamples");
		nburnin = b.getsizetvalue("NBurnIn");
		thinrate = b.getsizetvalue("ThinRate");
						
		if(b.getvalue("HighTemperature", temperature_high) == false){
			temperature_high = 2.5;
		}

		for (size_t i = 0; i < rjMcMCNuisance::number_of_types(); i++) {
			std::string str = strprint("Nuisance%lu", i + 1);
			cBlock c = b.findblock(str);
			if (c.Entries.size() == 0)break;

			rjMcMCNuisance n;
			std::vector<std::string> ninit;
			std::string typestr = c.getstringvalue("Type");
			n.settype(typestr);

			std::string sinit = c.getstringvalue("Initial");
			if (strcasecmp(sinit, "DataFile") == 0) {
				n.value = DBL_MIN;
			}
			else if (strncasecmp(sinit, "Column", 6) == 0) {
				n.value = DBL_MIN;
			}
			else if (strncasecmp(sinit, "-Column", 7) == 0) {
				n.value = DBL_MIN;
			}
			else {
				n.value = c.getdoublevalue("Initial");
			}

			n.min = c.getdoublevalue("Min");
			n.max = c.getdoublevalue("Max");
			if (n.min >= n.max) {
				printf("Error: Min >= Max for nuisance");
				exit(1);
			}
			if (n.max == DBL_MIN) {
				n.max = DBL_MAX;
			}
			n.sd_valuechange = c.getdoublevalue("STDValueChange");
			ntemplate.push_back(n);
			ninitial.push_back(sinit);
		}

	}
	
	bool readnextrecord()
	{
		if (filegetline(fp_indata, CurrentRecordString) == false) {
			return false;
		}
		else {
			CurrentRecord++;
			if (CurrentRecord > LastRecord + HeaderLines) {
				return false;
			}
			else {
				return true;
			}
		}
	}
	
	bool readnextrecord_thisprocess()
	{
		//Skip to next record for this process	
		while (readnextrecord()) {
			int  r = ((int)CurrentRecord - (int)HeaderLines - (int)FirstRecord);
			if (r < 0)continue;
			size_t n = (size_t)r % SubSample;
			size_t b = (size_t)floor((double)r / (double)SubSample);
			size_t p = b % mpiSize;
			if (n == 0 && p == mpiRank)return true;
		}
		return false;
	}
	
	void parsecurrentrecord()
	{
		CurrentRecordFields = fieldparsestring(CurrentRecordString.c_str(), " ,\t\r\n");
		std::vector<std::string>& f = CurrentRecordFields;

		fd_surveynumber.getvalue(f, surveynumber);
		fd_datenumber.getvalue(f, datenumber);
		fd_flightnumber.getvalue(f, flightnumber);
		fd_linenumber.getvalue(f, linenumber);
		fd_fidnumber.getvalue(f, fidnumber);
		fd_xord.getvalue(f, xord);
		fd_yord.getvalue(f, yord);
		fd_elevation.getvalue(f, elevation);

		for (size_t gi = 0; gi < IG.size(); gi++) {
			fd_geometry[gi].getvalue(f, IG[gi]);
		}

		for (size_t i = 0; i < nsystems; i++) {
			cTDEmSystemInfo& S = SV[i];

			if (S.useX) {
				S.fd_oPX.getvalue(f, S.oPX);
				S.fd_oSX.getvalue(f, S.oSX, S.nwindows);
				if (S.estimateNoiseFromModel) {
					for (size_t wi = 0; wi < S.nwindows; wi++) {
						double an = S.x_additivenoise[wi];
						double mn = 0.01 * S.x_multiplicativenoise * S.oSX[wi];
						S.oEX[wi] = sqrt(an * an + mn * mn);
					}
				}
				else {
					S.fd_oEX.getvalue(f, S.oEX, S.nwindows);
				}
			}

			if (S.useY) {
				S.fd_oPY.getvalue(f, S.oPY);
				S.fd_oSY.getvalue(f, S.oSY, S.nwindows);
				if (S.estimateNoiseFromModel) {
					for (size_t wi = 0; wi < S.nwindows; wi++) {
						double an = S.y_additivenoise[wi];
						double mn = 0.01 * S.y_multiplicativenoise * S.oSY[wi];
						S.oEY[wi] = sqrt(an * an + mn * mn);
					}
				}
				else {
					S.fd_oEY.getvalue(f, S.oEY, S.nwindows);
				}
			}

			if (S.useZ) {
				S.fd_oPZ.getvalue(f, S.oPZ);
				S.fd_oSZ.getvalue(f, S.oSZ, S.nwindows);
				if (S.estimateNoiseFromModel) {
					for (size_t wi = 0; wi < S.nwindows; wi++) {
						double an = S.z_additivenoise[wi];
						double mn = 0.01 * S.z_multiplicativenoise * S.oSZ[wi];
						S.oEZ[wi] = sqrt(an * an + mn * mn);
					}
				}
				else {
					S.fd_oEZ.getvalue(f, S.oEZ, S.nwindows);
				}
			}
		}
		set_data();
		set_nuisance();
	}
	
	void set_data()
	{
		size_t di = 0;
		for (size_t i = 0; i < nsystems; i++) {
			cTDEmSystemInfo& S = SV[i];
			cTDEmSystem& T = S.T;

			if (S.reconstructPrimary) {
				T.setgeometry(IG);
				T.LEM.calculation_type = CT_FORWARDMODEL;
				T.LEM.derivative_layer = INT_MAX;
				T.setprimaryfields();

				S.oPX = T.PrimaryX;
				S.oPY = T.PrimaryY;
				S.oPZ = T.PrimaryZ;
			}

			if (S.useX) {
				for (size_t wi = 0; wi < S.nwindows; wi++) {
					err[di] = S.oEX[wi];
					obs[di] = S.oSX[wi];
					if (S.useTotal)obs[di] += S.oPX;
					di++;
				}
			}
			if (S.useY) {
				for (size_t wi = 0; wi < S.nwindows; wi++) {
					err[di] = S.oEY[wi];
					obs[di] = S.oSY[wi];
					if (S.useTotal)obs[di] += S.oPY;
					di++;
				}
			}
			if (S.useZ) {
				for (size_t wi = 0; wi < S.nwindows; wi++) {
					err[di] = S.oEZ[wi];
					obs[di] = S.oSZ[wi];
					if (S.useTotal)obs[di] += S.oPZ;
					di++;
				}
			}
		}
	}
	
	void set_nuisance()
	{
		nuisance_init = ntemplate;

		for (size_t i = 0; i < ntemplate.size(); i++) {

			rjMcMCNuisance& n = nuisance_init[i];

			std::string s = ninitial[i];

			if (strcasecmp(s, "DataFile") == 0) {
				n.value = DBL_MIN;
			}

			if (strncasecmp(s, "Column", 6) == 0) {
				int col;
				sscanf(ninitial[i].c_str(), "Column %d", &col);
				col = col - 1;//referenced from 1 not 0
				double v;
				sscanf(CurrentRecordFields[col].c_str(), "%lf", &v);
				n.value = v;
			}
			else if (strncasecmp(s, "-Column", 7) == 0) {
				int col;
				sscanf(ninitial[i].c_str(), "-Column %d", &col);
				col = col - 1;//referenced from 1 not 0
				double v;
				sscanf(CurrentRecordFields[col].c_str(), "%lf", &v);
				n.value = -v;
			}
			else {
				//Stay at initial constant
				n.value = n.value;
			}


			if (n.value == DBL_MIN) {
				switch (n.type) {
				case rjMcMCNuisance::Type::TX_HEIGHT:
					n.value = IG.tx_height; break;
				case rjMcMCNuisance::Type::TX_ROLL:
					n.value = IG.tx_roll; break;
				case rjMcMCNuisance::Type::TX_PITCH:
					n.value = IG.tx_pitch; break;
				case rjMcMCNuisance::Type::TX_YAW:
					n.value = IG.tx_yaw;	break;
				case rjMcMCNuisance::Type::TXRX_DX:
					n.value = IG.txrx_dx; break;
				case rjMcMCNuisance::Type::TXRX_DY:
					n.value = IG.txrx_dy; break;
				case rjMcMCNuisance::Type::TXRX_DZ:
					n.value = IG.txrx_dz; break;
				case rjMcMCNuisance::Type::RX_ROLL:
					n.value = IG.rx_roll; break;
				case rjMcMCNuisance::Type::RX_PITCH:
					n.value = IG.rx_pitch; break;
				case rjMcMCNuisance::Type::RX_YAW:
					n.value = IG.rx_yaw; break;
				case rjMcMCNuisance::Type::TXRX_DISTANCE:
					n.value = sqrt(IG.txrx_dx * IG.txrx_dx + IG.txrx_dz * IG.txrx_dz);
					break;
				case rjMcMCNuisance::Type::TXRX_ANGLE:
					n.value = R2D * atan2(IG.txrx_dz, IG.txrx_dx);
					break;
				default:break;
				}
			}
		}
	}
	
	void sample()
	{		
		rjMcMC1DSampler::reset();		
		rjMcMC1DSampler::sample();
		std::string dstr = results_string();

		sFilePathParts fpp = getfilepathparts(OutputDataFile);
		std::string hdrfile = fpp.directory + fpp.prefix + ".hdr";
		std::string aseggdffile = fpp.directory + fpp.prefix + ".dfn";

		//Output header file
		if (exists(hdrfile) == false) {
			OI.write_simple_header(hdrfile);
			OI.write_aseggdf_header(aseggdffile);
		}

		//Output data record	
		FILE* fp = fileopen(OutputDataFile, "a");
		fprintf(fp, dstr.c_str());
		fclose(fp);
		
		write_maps_to_file_netcdf();
		
	}
	
	std::string results_string()
	{
		size_t ndepthcells = pmap.npbins();
		std::vector<cHistogramStats<double>> hs = pmap.hstats();
		std::vector<double> hlike = pmap.modelmap(HighestLikelihood);
		std::vector<double> lmfit = pmap.modelmap(LowestMisfit);
		if (param_value.islog10()) {
			for (size_t j = 0; j < ndepthcells; j++) {
				hs[j].p10 = pow10(hs[j].p10);
				hs[j].p50 = pow10(hs[j].p50);
				hs[j].p90 = pow10(hs[j].p90);
				hs[j].mode = pow10(hs[j].mode);
				hs[j].mean = pow10(hs[j].mean);
				hlike[j] = pow10(hlike[j]);
				lmfit[j] = pow10(lmfit[j]);
			}
		}
		double misfit_lowest = LowestMisfit.get_misfit() / double(ndata);

		size_t n = 0;
		double sum = 0.0;		
		for (size_t mi = 0; mi < ensemble.size(); mi++) {
			sum += ensemble[mi].get_misfit();
			n++;
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

		///////////////////

		OI.addfield("nchains", 'I', 9, 0);
		OI.setcomment("Number of chains");
		buf += strprint("%9lu", nchains());

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
		for (size_t j = 0; j < ndepthcells; j++) {
			double dcc = pmap.toppbin(j);
			if (param_position.islog10()) {
				dcc = pow10(dcc);
			}
			buf += strprint("%8.2lf", dcc);
		}


		size_t width = 15;
		size_t decimals = 6;
		std::string fmt = "%15.6le";

		OI.addfield("conductivity_mean", 'E', width, decimals, ndepthcells);
		OI.setunits("S/m"); OI.setcomment("Conductivity of mean model");
		for (size_t j = 0; j < ndepthcells; j++) {
			buf += strprint(fmt.c_str(), hs[j].mean);
		}

		OI.addfield("conductivity_mode", 'E', width, decimals, ndepthcells);
		OI.setunits("S/m"); OI.setcomment("Conductivity of mode model");
		for (size_t j = 0; j < ndepthcells; j++) {
			buf += strprint(fmt.c_str(), hs[j].mode);
		}

		OI.addfield("conductivity_p50", 'E', width, decimals, ndepthcells);
		OI.setunits("S/m"); OI.setcomment("Conductivity at 50th percentile or median");
		for (size_t j = 0; j < ndepthcells; j++) {
			buf += strprint(fmt.c_str(), hs[j].p50);
		}

		OI.addfield("conductivity_p10", 'E', width, decimals, ndepthcells);
		OI.setunits("S/m"); OI.setcomment("Conductivity at 10th percentile");
		for (size_t j = 0; j < ndepthcells; j++) {
			buf += strprint(fmt.c_str(), hs[j].p10);
		}

		OI.addfield("conductivity_p90", 'E', width, decimals, ndepthcells);
		OI.setunits("S/m"); OI.setcomment("Conductivity at 90th percentile");
		for (size_t j = 0; j < ndepthcells; j++) {
			buf += strprint(fmt.c_str(), hs[j].p90);
		}

		OI.addfield("conductivity_highestlikelihood", 'E', width, decimals, ndepthcells);
		OI.setunits("S/m"); OI.setcomment("Conductivity of highest likelihood model");
		for (size_t j = 0; j < ndepthcells; j++) {
			buf += strprint(fmt.c_str(), hlike[j]);
		}

		OI.addfield("conductivity_lowestmisfit", 'E', width, decimals, ndepthcells);
		OI.setunits("S/m"); OI.setcomment("Conductivity of howest misfit model");
		for (size_t j = 0; j < ndepthcells; j++) {
			buf += strprint(fmt.c_str(), lmfit[j]);
		}

		OI.addfield("changepoint", 'I', 11, 0, ndepthcells);
		OI.setunits("counts"); OI.setcomment("Changepoint histogram");
		for (size_t j = 0; j < ndepthcells; j++) {
			buf += strprint("%11lu", pmap.changepoint()[j]);
		}

		size_t nn = nnuisances();
		for (size_t j = 0; j < nn; j++) {
			std::string nstr = ensemble[0].nuisances[j].typestring();
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

	std::string prefixstring()
	{
		std::string s = "seq";
		s += strprint(".%08lu", CurrentRecord);
		s += strprint(".%lu", linenumber);
		s += strprint(".%lf", fidnumber);
		return s;
	}
	
	void write_maps_to_file_netcdf()
	{
		if (SaveMaps == false)return;
		if ((CurrentRecord - HeaderLines - FirstRecord) / SubSample % SaveMapsRate != 0)return;

		std::string fileprefix = prefixstring();
		std::string fname = MapsDirectory + fileprefix + ".pmap";		
		std::string ncfilepath = MapsDirectory + fileprefix + ".nc";

		NcFile nc(ncfilepath, NcFile::FileMode::replace);
		NcGroupAtt a;
		a = nc.putAtt("survey", NcType::nc_DOUBLE, (double) surveynumber);
		a = nc.putAtt("date", NcType::nc_DOUBLE, (double) datenumber);
		a = nc.putAtt("flight", NcType::nc_DOUBLE, (double) flightnumber);
		a = nc.putAtt("line", NcType::nc_DOUBLE, (double) linenumber);
		a = nc.putAtt("fiducial", NcType::nc_DOUBLE, (double) fidnumber);
		a = nc.putAtt("time", NcType::nc_DOUBLE, (double) timenumber);
		a = nc.putAtt("x", NcType::nc_DOUBLE, xord);
		a = nc.putAtt("y", NcType::nc_DOUBLE, yord);
		a = nc.putAtt("elevation", NcType::nc_DOUBLE, elevation);
		rjMcMC1DSampler::writemapstofile_netcdf(nc);
	}
	
	cTDEmGeometry getgeometry(const rjMcMC1DModel& m)
	{
		cTDEmGeometry  OG = IG;

		bool angledistance = false;
		double angle = 0.0;
		double distance = 0.0;

		for (size_t i = 0; i < m.nuisances.size(); i++) {
			const rjMcMCNuisance& n = m.nuisances[i];
			switch (n.type) {
			case rjMcMCNuisance::Type::TX_HEIGHT:
				OG.tx_height = n.value; break;
			case rjMcMCNuisance::Type::TX_ROLL:
				OG.tx_roll = n.value; break;
			case rjMcMCNuisance::Type::TX_PITCH:
				OG.tx_pitch = n.value; break;
			case rjMcMCNuisance::Type::TX_YAW:
				OG.tx_yaw = n.value; break;
			case rjMcMCNuisance::Type::TXRX_DX:
				OG.txrx_dx = n.value; break;
			case rjMcMCNuisance::Type::TXRX_DY:
				OG.txrx_dy = n.value; break;
			case rjMcMCNuisance::Type::TXRX_DZ:
				OG.txrx_dz = n.value; break;
			case rjMcMCNuisance::Type::RX_ROLL:
				OG.rx_roll = n.value; break;
			case rjMcMCNuisance::Type::RX_PITCH:
				OG.rx_pitch = n.value; break;
			case rjMcMCNuisance::Type::RX_YAW:
				OG.rx_yaw = n.value; break;
			case rjMcMCNuisance::Type::TXRX_DISTANCE:
				angledistance = true;
				distance = n.value; break;
			case rjMcMCNuisance::Type::TXRX_ANGLE:
				angledistance = true;
				angle = n.value; break;
			default:
				break;
			}
		}

		if (angledistance == true) {
			OG.txrx_dx = distance * cos(D2R * angle);
			OG.txrx_dz = distance * sin(D2R * angle);
		}
		return OG;
	}
	
	std::vector<double> collect(const cTDEmSystemInfo& S, const cTDEmSystem& T)
	{
		std::vector<double> v(S.nchans);
		std::vector<double> x, y, z;
		if (S.useTotal) {
			if (S.useX) x = T.X + T.PrimaryX;
			if (S.useY) y = T.Y + T.PrimaryY;
			if (S.useZ) z = T.Z + T.PrimaryZ;
		}
		else {
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
	
	std::vector<double> forwardmodel(const rjMcMC1DModel& m)
	{
		std::vector<double> c = m.getvalues();
		if (param_value.islog10()){
			pow10_apply(c);
		}

		std::vector<double> t = m.getthicknesses();
		cTDEmGeometry  G = getgeometry(m);
		std::vector<double> pred(ndata);

		size_t di = 0;
		for (size_t i = 0; i < nsystems; i++) {
			cTDEmSystemInfo& S = SV[i];
			cTDEmSystem& T = S.T;
			T.setconductivitythickness(c, t);
			T.setgeometry(G);
			T.setupcomputations();
			T.LEM.calculation_type = CT_FORWARDMODEL;
			T.LEM.derivative_layer = INT_MAX;
			T.setprimaryfields();
			T.setsecondaryfields();
			std::vector<double> v = collect(S, T);
			for (size_t j = 0; j < v.size(); j++) {
				pred[di] = v[j];
				di++;
			}
		}
		return pred;
	}
};

#endif
