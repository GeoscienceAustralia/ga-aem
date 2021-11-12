/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _csbsinverter_H
#define _csbsinverter_H

#include <stdio.h>
#include <sstream>
#include <vector>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <functional>
#include <variant>

#include "general_types.h"
#include "cinverter.h"
#include "airborne_types.h"
#include "tdemsystem.h"
#include "string_utils.h"
#include "vector_utils.h"

class cGeomStruct {

public:
	cTDEmGeometry input;
	cTDEmGeometry ref;
	cTDEmGeometry std;
	cTDEmGeometry min;
	cTDEmGeometry max;
	cTDEmGeometry tfr;
	cTDEmGeometry invmodel;
};

class cEarthStruct {

public:	
	cEarth1D ref;
	cEarth1D std;
	cEarth1D min;
	cEarth1D max;
	cEarth1D invmodel;

	void sanity_check() {

		size_t nc = ref.conductivity.size();
		size_t nt = ref.thickness.size();

		std::ostringstream oss;
		if (nc != nt + 1) {
			oss << "The conductivity and/or thickness do not have the correct number of layers\n";
		}

		if (ref.conductivity.size() > 0) {
			if (::min(ref.conductivity) <= 0) oss << "The conductivity ref is <= 0 in at least one layer\n";
		}	

		if (std.conductivity.size() > 0) {
			if (::min(std.conductivity) <= 0) oss << "The conductivity std is <= 0\n";
		}

		if (min.conductivity.size() > 0) {		
			if (min.conductivity.size() != nc) oss << "The conductivity min does not have the correct number of layer\n";
			if (max.conductivity.size() != nc) oss << "The conductivity max does not have the correct number of layer\n";
			if (::min(min.conductivity) <= 0) oss << "The conductivity min is <= 0 in at least one layer in at least one layer\n";
			if (::min(max.conductivity) <= 0) oss << "The conductivity max is <= 0 in at least one layer in at least one layer\n";
			if (::min(max.conductivity - min.conductivity) <= 0) oss << "The conductivity max <= min in at least one layer\n";
			if (::min(ref.conductivity - min.conductivity) <= 0) oss << "The conductivity ref <= min in at least one layer\n";
			if (::min(max.conductivity - ref.conductivity) <= 0) oss << "The conductivity ref >= max in at least one layer\n";
		}

		if (ref.thickness.size() > 0) {
			if (::min(ref.thickness) <= 0) oss << "The thickness ref is <= 0 in at least one layer\n";
		}

		if (std.thickness.size() > 0) {
			if (::min(std.thickness) <= 0) oss << "The thickness std is <= 0 in at least one layer\n";
		}

		if (min.thickness.size() > 0) {
			if (min.thickness.size() != nc) oss << "The thickness min does not have the correct number of layer\n";
			if (max.thickness.size() != nc) oss << "The thickness max does not have the correct number of layer\n";
			if (::min(min.thickness) <= 0) oss << "The thickness min is <= 0 in at least one layer\n";
			if (::min(max.thickness) <= 0) oss << "The thickness max is <= 0 in at least one layer\n";
			if (::min(max.thickness - min.thickness) <= 0) oss << "The thickness max <= min in at least one layer\n";
			if (::min(ref.thickness - min.thickness) <= 0) oss << "The thickness ref <= min in at least one layer\n";
			if (::min(max.thickness - ref.thickness) <= 0) oss << "The thickness ref >= max in at least one layer\n";
		}

		if (oss.str().size() > 0) {
			glog.errormsg(oss.str());
		}

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

/*
class cComponentInfo1 {

public:

	std::string Name;
	bool Use = false;
	double  oP = 0.0;
	std::vector<double>  oS;
	std::vector<double>  oE;
	cFieldDefinition fd_oP;
	cFieldDefinition fd_oS;
	cFieldDefinition fd_oE;
	bool EstimateNoiseFromModel = false;
	std::vector<double> mn;
	std::vector<double> an;
	int dataindex = -1;

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
				glog.errormsg(_SRC_, "Must have exactly 1 or nwindows AdditiveNoise values\n");
			};

			if (mn.size() == 1) {
				mn = std::vector<double>(nwindows, mn[0]);
			}
			if (mn.size() != nwindows) {
				glog.errormsg(_SRC_, "Must have exactly 1 or nwindows MultiplicativeNoise values\n");
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
				const double v = 0.01 * mn[w] * oS[w];
				oE[w] = std::hypot(an[w], v);
			}
		}
		else {
			IM->read(fd_oE, oE, nw());
		}
	}
};
*/

class cComponentInfo {

public:

	std::string Name;
	bool Use = false;
	double  oP = 0.0;
	std::vector<double>  oS;
	std::vector<double>  oE;
	cFieldDefinition fd_oP;
	cFieldDefinition fd_oS;
	cFieldDefinition fd_oE;
	bool EstimateNoiseFromModel=false;
	std::vector<double> mn;
	std::vector<double> an;
	int dataindex = -1;

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
			else if (an.size() != nwindows) {				
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
	size_t nwindows=0;
	size_t ncomps=0;
	size_t nchans=0;
	cComponentInfo CompInfo[3];
	int xzIndex=-1;//Start index of XZ in data array

	bool invertXPlusZ=false;
	bool invertPrimaryPlusSecondary=false;
	bool reconstructPrimary=false;		
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

class cSBSInverter : public cInverter {

private:	
	
	Vector cull(const Vector& vall) const {
		assert(ActiveData.size() == nData);
		assert(vall.size() == nAllData);
		Vector vcull(nData);
		for (size_t i = 0; i < nData; i++) {
			vcull[i] = vall[ActiveData[i]];
		}
		return vcull;
	}

	Vector cull(const std::vector<double>& vall) const {
		assert(ActiveData.size() == nData);
		assert(vall.size() == nAllData);
		Vector vcull(nAllData);
		for (size_t i = 0; i < nData; i++) {
			vcull[i] = vall[ActiveData[i]];
		}
		return vcull;
	}

	Matrix cull(const Matrix& mall) const {
		assert(ActiveData.size() == nData);
		assert(mall.rows() == nAllData);
		assert(mall.cols() == nParam);
		Matrix mcull(nData, nParam);
		for (size_t i = 0; i < nData; i++) {
			mcull.row(i) = mall.row(ActiveData[i]);
		}
		return mcull;
	}

public:			
			
	using cIFDMap = cKeyVec<std::string, cInvertibleFieldDefinition, caseinsensetiveequal<std::string>>;

	int    BeginGeometrySolveIteration=0;
	bool   FreeGeometry = false;
	size_t nAllData = 0;
	
	Matrix Wc;
	Matrix Wt;
	Matrix Wg;
	Matrix Wr;
	Matrix Ws;
	Matrix Wq;

	double AlphaC=0.0;
	double AlphaT=0.0;
	double AlphaG=0.0;
	double AlphaS=0.0;
	double AlphaQ=0.0;

	cSBSInverter(const std::string& controlfile, const int& size, const int& rank, const bool& usingopenmp, const std::string commandline) 
					: cInverter(controlfile, size, rank, usingopenmp, commandline)
	{
		std::cout << "Constructing cSBSInverter\n";
		_GSTPUSH_
		try {			
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

	~cSBSInverter() {
		std::cout << "Destroying cSBSInverter\n";
	}
	
	cOutputOptions OO;
	std::vector<cTDEmSystemInfo> SV;		
	size_t nsystems = 0;			
	size_t pointsoutput = 0;
	size_t nlayers = 0;
	size_t ngeomparam = 0;
	int _gindex_ = -1;
			
	//Column definitions		
	cInvertibleFieldDefinition fdC;
	cInvertibleFieldDefinition fdT;					
	cIFDMap fdG;
	
	//Sample instances
	struct SampleId {
		int uniqueid = -1;
		int survey = -1;
		int date = -1;
		int flight = -1;
		int line = -1;
		double fiducial = -1.0;
		double x = -1.0;
		double y = -1.0;
		double elevation = 0.0;
	} Id;
			
	cKeyVec<std::string,cFdVrnt,caseinsensetiveequal<std::string>> AncFld;
	cGeomStruct G;		
	cEarthStruct E;
		
	void loadcontrolfile(const std::string& filename)
	{		
		glog.logmsg(0, "Loading control file %s\n", filename.c_str());
		Control = cBlock(filename);
		cBlock ob = Control.findblock("Output");
		cBlock ib = Control.findblock("Input");

		OO = cOutputOptions(ob);
		Verbose = ob.getboolvalue("verbose");

		std::string suffix = stringvalue(Rank, ".%04d");
		OO.LogFile = insert_after_filename(OO.LogFile, suffix);
		openlogfile(); //load this first to get outputlogfile opened

		//Load control file
		parseoptions();
		initialise_systems();


		if (cInputManager::isnetcdf(ib)) {
#if !defined HAVE_NETCDF
			glog.errormsg(_SRC_, "Sorry NETCDF I/O is not available in this executable\n");
#endif			
			IM = std::make_unique<cNetCDFInputManager>(ib);
			std::string s = IM->datafilename();
		}
		else {
			IM = std::make_unique<cASCIIInputManager>(ib);
		}

		if (cOutputManager::isnetcdf(ob)) {
#if !defined HAVE_NETCDF
			glog.errormsg(_SRC_, "Sorry NETCDF I/O is not available in this executable\n");
#endif			
			OM = std::make_unique<cNetCDFOutputManager>(ob, Size, Rank);
		}
		else {
			OM = std::make_unique<cASCIIOutputManager>(ob, Size, Rank);
		}
		OM->opendatafile(IM->datafilename(), IM->subsamplerate());
	}

	bool solve_thickness() const
	{
		return fdT.solve;
	};
	
	bool solve_conductivity() const {
		return fdC.solve;
	};
	
	bool solve_geometry(const std::string& e) {						
		return fdG.cref(e).solve;
	};
	
	bool solve_geometry() const {
		if (_gindex_ == -1) return false;
		return true;
	};
	
	std::string dumppath() const
	{
		std::string s = OO.DumpPath(IM->record(), CIS.iteration);
		return s;
	};

	void dump_record_number() {
		std::ofstream of(dumppath() + "record.dat");
		of << "Record\t" << IM->record() << std::endl;;
	}

	const size_t cindex(const size_t& li) {		
		if (solve_conductivity() == false) {
			glog.errormsg("Out of boundes in cindex()\n");
		}
		return fdC.index + li;
	}
	
	const size_t tindex(const size_t& li) {
		if (solve_thickness() == false) {
			glog.errormsg("Out of boundes in tindex()\n");
		}
		return fdT.index + li;
	}
	
	const size_t gindex(const size_t& pi) {
		if (solve_geometry() == false) {
			glog.errormsg("Out of boundes in gindex\n");
		}
		return pi + _gindex_;
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
	
	void parseoptions()
	{
		cBlock b = Control.findblock("Options");		
		AlphaC = b.getdoublevalue("AlphaConductivity");
		AlphaT = b.getdoublevalue("AlphaThickness");
		AlphaG = b.getdoublevalue("AlphaGeometry");		
		AlphaS = b.getdoublevalue("AlphaSmoothness");
		AlphaQ = b.getdoublevalue("AlphaHomogeneous");
		
		BeginGeometrySolveIteration = b.getintvalue("BeginGeometrySolveIteration");						
		if (!isdefined(BeginGeometrySolveIteration)) {
			BeginGeometrySolveIteration = 0;
		}

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
		cBlock b = Control.findblock("Input.AncillaryFields");
		set_field_definitions_ancillary(b);
		if (AncFld.keyindex("line") < 0) {
			glog.errormsg("Must specify a linenumber field\n");
		}

		b = Control.findblock("Input.Geometry");
		fdG = set_field_definitions_geometry(b);

		b = Control.findblock("Input.Earth");
		fdC = cInvertibleFieldDefinition(b, "Conductivity");
		fdT = cInvertibleFieldDefinition(b, "Thickness");
	}

	void set_field_definitions_ancillary(const cBlock& parent) {
		const cBlock& b = parent;
		for (size_t i = 0; i < b.Entries.size(); i++) {
			std::string key   = b.key(i);
			std::string value = b.value(i);
			cFieldDefinition fd(parent, key);						
			cFdVrnt fdvrnt(fd,cVrnt());						
			IM->set_variant_type(fd.varname, fdvrnt.vnt);
			AncFld.add(key, fdvrnt);
		}		
	}
	
	cIFDMap set_field_definitions_geometry(const cBlock& parent)
	{				
		cIFDMap g;
		for (size_t i = 0; i < cTDEmGeometry::size(); i++) {			
			std::string key = cTDEmGeometry::element_name(i);
			cInvertibleFieldDefinition f(parent, key);
			bool a = g.add(key, f);
			if (a == false) {
				std::string msg = strprint("Parameter %s has already been already added\n", key.c_str());
				glog.errormsg(msg);
			}
		}
		return g;
	}
	
	void initialise_systems()
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
		nAllData = 0;
		for (size_t si = 0; si < nsystems; si++) {
			cTDEmSystemInfo& S = SV[si];
			if (S.invertXPlusZ) {
				S.xzIndex = (int) nAllData;
				S.CompInfo[0].dataindex = -1;
				S.CompInfo[2].dataindex = -1;
				nAllData += S.nwindows;

				if (S.CompInfo[1].Use) {
					S.CompInfo[1].dataindex = (int) nAllData;
					nAllData += S.nwindows;
				}
			}
			else {
				for (size_t i = 0; i < 3; i++) {
					if (S.CompInfo[i].Use) {
						S.CompInfo[i].dataindex = (int) nAllData;
						nAllData += S.nwindows;
					}
				}
			}
		}				
	}

	void setup_parameters()
	{
		bool status = Control.getvalue("Input.Earth.Conductivity.NumberOfLayers",nlayers);
		if (status == false) {
			std::stringstream msg;
			msg << "The NumberOfLayers must be specified in Input.Columns.Conductivity\n";
			glog.errormsg(msg.str());
		}

		nParam = 0;
		ngeomparam = 0;
		if (solve_conductivity()) {			
			fdC.index = (int)nParam;			
			nParam += nlayers;
		}

		if (solve_thickness()) {
			fdT.index = (int)nParam;
			nParam += nlayers - 1;
		}


		//Geometry params	
		_gindex_ = (int)nParam;
		for (size_t i = 0; i < cTDEmGeometry::size(); i++) {
			std::string gname = cTDEmGeometry::element_name(i);
			cInvertibleFieldDefinition& a = fdG.cref(gname);
			if (a.solve) {
				a.index = (int)nParam;
				nParam++;
				ngeomparam++;
			}
			else {
				a.index = -1;
			}
		}		
		if (ngeomparam == 0) _gindex_ = -1;//reset if none;
		
		RefParam.resize(nParam);
		RefParamStd.resize(nParam);
	}
	
	bool read_record()
	{		
		if (IM->parserecord() == false) return false;		
		bool readstatus = true;
		bool status;				
		Id.uniqueid = (int) IM->record();
		
		status = read_ancillary_fields();
		
		status = readgeometry(fdG);
		
		status = IM->read(fdC.input, E.ref.conductivity, nlayers); if (status == false) readstatus = false;				
		if (solve_conductivity()) {
			status = IM->read(fdC.ref, E.ref.conductivity, nlayers); if (status == false) readstatus = false;
			status = IM->read(fdC.std, E.std.conductivity, nlayers); if (status == false) readstatus = false;
			status = IM->read(fdC.min, E.min.conductivity, nlayers); if (status == false) readstatus = false;
			status = IM->read(fdC.max, E.max.conductivity, nlayers); if (status == false) readstatus = false;
		}
		
		status = IM->read(fdT.input, E.ref.thickness, nlayers - 1); if (status == false) readstatus = false;
		if (solve_thickness()) {
			status = IM->read(fdT.ref, E.ref.thickness, nlayers - 1); if (status == false) readstatus = false;
			status = IM->read(fdT.std, E.std.thickness, nlayers - 1); if (status == false) readstatus = false;
			status = IM->read(fdT.min, E.min.thickness, nlayers - 1); if (status == false) readstatus = false;
			status = IM->read(fdT.max, E.max.thickness, nlayers - 1); if (status == false) readstatus = false;
		}
		E.sanity_check();

		for (size_t si = 0; si < nsystems; si++) {
			read_system_data(si);
		}
		return readstatus;
	}

	template<typename T>
	bool set_ancillary_id(const std::string key, T& value) {
		int ki = AncFld.keyindex(key);
		if (ki >= 0) {
			value = std::get<T>(AncFld[ki].second.vnt);
			return true;
		}
		return false;
	}

	bool read_ancillary_fields() {
		for (size_t i = 0; i < AncFld.size(); i++) {
			IM->readfdvnt(AncFld[i].second);
		}
		set_ancillary_id("Survey", Id.survey);
		set_ancillary_id("Date", Id.date);
		set_ancillary_id("Flight", Id.flight);
		set_ancillary_id("Line", Id.line);
		set_ancillary_id("Fiducial", Id.fiducial);
		set_ancillary_id("X", Id.x);
		set_ancillary_id("Y", Id.y);
		set_ancillary_id("GroundElevation", Id.elevation);
		return true;
	}

	void read_system_data(size_t sysindex)
	{
		cTDEmSystemInfo& S = SV[sysindex];
		S.CompInfo[0].readdata(IM);
		S.CompInfo[1].readdata(IM);
		S.CompInfo[2].readdata(IM);
	}
	
	bool initialise_data(){
		std::vector<double> obs(nAllData);
		std::vector<double> err(nAllData);
		std::vector<double> pred(nAllData);
		for (size_t si = 0; si < nsystems; si++) {
			cTDEmSystemInfo& S = SV[si];
			cTDEmSystem& T = S.T;
			if (S.reconstructPrimary) {
				T.setgeometry(G.tfr);
				T.LEM.calculation_type = cLEM::CalculationType::FORWARDMODEL;
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
		ActiveData.clear();
		for (size_t i = 0; i < nAllData; i++){
			if(!isnull(obs[i])  && !isnull(err[i])) ActiveData.push_back(i);			
		}
		nData = ActiveData.size();

		if (nData != nAllData) {
			size_t ncull = nAllData - nData;
			OutputMessage += strprint(", %d null data/noise were culled",(int)ncull);
		}
		Err = cull(err);
		Obs = cull(obs);

		//Check for zero Error values		
		int nzeroerr = 0;
		for (size_t i = 0; i < nData; i++) {
			if (Err[i] == 0.0) nzeroerr++;
		}
		if (nzeroerr > 0) {			
			OutputMessage += strprint(", Skipped %d noise values were 0.0", nzeroerr);			
			return false;
		}				
		return true;
	}
		
	void initialise_parameters(){
		if (solve_conductivity()) {
			for (size_t i = 0; i < nlayers; i++) {
				RefParam[cindex(i)] = log10(E.ref.conductivity[i]);
				RefParamStd[cindex(i)] = E.std.conductivity[i];
			}
		}

		if (solve_thickness()) {
			for (size_t i = 0; i < nlayers - 1; i++) {
				RefParam[tindex(i)] = log10(E.ref.thickness[i]);
				RefParamStd[tindex(i)] = E.std.thickness[i];
			}
		}

		for (int i = 0; i < cTDEmGeometry::size(); i++) {
			std::string gname = cTDEmGeometry::element_name(i);
			auto a = fdG.cref(gname);
			if (a.index > 0) {
				RefParam[a.index] = G.ref[gname];
				RefParamStd[a.index] = G.std[gname];
			}
		}		
	}

	void initialise_Wc(){
		Wc = Matrix::Zero(nParam, nParam);		
		if (solve_conductivity() == false)return;

		std::vector<double> t(nlayers);
		if (nlayers == 1) {
			t[0] = 1;
		}
		else if (nlayers == 2) {
			t[0] = E.ref.thickness[0];
			t[1] = E.ref.thickness[0];
		}
		else {
			for (size_t i = 0; i < (nlayers - 1); i++) {
				t[i] = E.ref.thickness[i];
			}
			t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3])*t[nlayers - 2];
		}


		double tsum = 0.0;
		for (size_t i = 0; i < nlayers; i++)tsum += t[i];
		double tavg = tsum / (double)nlayers;

		double s = AlphaC / (double)(nlayers);
		for (size_t i = 0; i < nlayers; i++) {
			size_t p = cindex(i);
			Wc(p,p) = s * (t[i] / tavg) / (RefParamStd[p] * RefParamStd[p]);
		}
	}

	void initialise_Wt(){
		Wt = Matrix::Zero(nParam, nParam);
		if (solve_thickness() == false)return;

		double s = AlphaT / (double)(nlayers - 1);
		for (size_t i = 0; i < nlayers - 1; i++) {
			size_t p = tindex(i);
			Wt(p,p) = s / (RefParamStd[p] * RefParamStd[p]);
		}

	}

	void initialise_Wg(){		
		Wg = Matrix::Zero(nParam, nParam);
		if (ngeomparam <= 0)return;

		double s = AlphaG / (double)ngeomparam;
		for (size_t i = 0; i < ngeomparam; i++) {
			size_t p = gindex(i);
			Wg(p,p) = s / (RefParamStd[p] * RefParamStd[p]);
		}
	}
		
	void initialise_L_Ws_1st_derivative()
	{
		Ws = Matrix::Zero(nParam, nParam);
		if (AlphaS == 0 || nlayers < 3) return;
		if (solve_conductivity() == false) return;

		std::vector<double> t(nlayers);
		for (size_t i = 0; i < (nlayers - 1); i++) {
			t[i] = E.ref.thickness[i];
		}
		t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3])*t[nlayers - 2];


		double tsum = 0.0;
		for (size_t i = 0; i < nlayers; i++)tsum += t[i];
		double tavg = tsum / (double)nlayers;

		Matrix L = Matrix::Zero(nlayers - 1, nParam);
		size_t neqn = 0;
		for (size_t li = 1; li < nlayers; li++) {
			const size_t& pindex = cindex(li);
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
		Ws = Matrix::Zero(nParam, nParam);
		if (AlphaS == 0 || nlayers < 3) return;
		if (solve_conductivity() == false) return;

		std::vector<double> t(nlayers);
		if (nlayers == 1) {
			t[0] = 1.0;
		}
		else {
			for (size_t i = 0; i < (nlayers - 1); i++) {
				t[i] = E.ref.thickness[i];
			}
			t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3])*t[nlayers - 2];
		}

		double tsum = 0.0;
		for (size_t i = 0; i < nlayers; i++)tsum += t[i];
		double tavg = tsum / (double)nlayers;

		Matrix L = Matrix::Zero(nlayers - 2, nParam);
		size_t neqn = 0;
		for (size_t li = 1; li < nlayers - 1; li++) {
			const size_t& pindex = cindex(li);
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
	
	void initialise_Ws() {
		if (SmoothnessMethod == eSmoothnessMethod::DERIVATIVE_1ST) {
			initialise_L_Ws_1st_derivative();
		}
		else if (SmoothnessMethod == eSmoothnessMethod::DERIVATIVE_2ND) {
			initialise_L_Ws_2nd_derivative();
		}
	}

	void initialise_Wq()
	{
		Wq = Matrix::Zero(nParam, nParam);
		if (AlphaQ == 0) return;
		if (solve_conductivity() == false) return;

		std::vector<double> t(nlayers);
		for (size_t i = 0; i < (nlayers - 1); i++) {
			t[i] = E.ref.thickness[i];
		}
		t[nlayers - 1] = (t[nlayers - 2] / t[nlayers - 3]) * t[nlayers - 2];


		double tsum = 0.0;
		for (size_t i = 0; i < nlayers; i++)tsum += t[i];
		double tavg = tsum / (double)nlayers;

		Matrix L = Matrix::Zero(nlayers, nParam);
		size_t neqn = 0;
		for (size_t li = 0; li < nlayers; li++) {
			const size_t& lpindex = cindex(li);
			for (size_t ki = 0; ki < nlayers; ki++) {
				const size_t& kpindex = cindex(ki);				
				double s = std::sqrt(t[li] / tavg);//sqrt because it gets squared in L'L				
				if (li == ki) {
					L(lpindex, kpindex) = 1.0;
				}
				else {
					L(lpindex, kpindex) = -1.0 / ((double)nlayers - 1);
				}
				neqn++;
			}
		}
		Wq = L.transpose() * L;
		Wq *= (AlphaQ / (double)(nlayers));
		std::cerr << Wq;
	}

	void initialise_Wr() {
		initialise_Wc();
		initialise_Wt();
		initialise_Wg();

		Wr = Matrix::Zero(nParam, nParam);
		if (AlphaC > 0.0) Wr += Wc;
		if (AlphaT > 0.0) Wr += Wt;
		if (AlphaG > 0.0) Wr += Wg;
	}

	void initialise_Wm(){
		initialise_Wq();
		initialise_Ws();
		initialise_Wr();		
		Wm = Wr + Ws + Wq;
	}

	void dump_W_matrices(){		
		if (OO.Dump) {
			writetofile(Wc, dumppath() + "Wc.dat");
			writetofile(Wt, dumppath() + "Wt.dat");
			writetofile(Wg, dumppath() + "Wg.dat");
			writetofile(Wr, dumppath() + "Wr.dat");
			writetofile(Ws, dumppath() + "Ws.dat");
			writetofile(Wm, dumppath() + "Wm.dat");
			writetofile(Wd, dumppath() + "Wd.dat");
		}
	}
	
	Vector parameter_change(const double& lambda, const Vector& m_old, const Vector& pred)
	{
		Vector m_new = solve_linear_system(lambda, m_old, pred);
		Vector dm = m_new - m_old;

		if (fdC.bound()) {			
			for (size_t li = 0; li < nlayers; li++) {
				const size_t pindex = cindex(li);				
				const double lmin = std::log10(E.min.conductivity[li]);
				const double lmax = std::log10(E.max.conductivity[li]);
				if(m_new[pindex] < lmin) {
					if (Verbose) {
						std::cerr << rec_it_str() << std::endl;
						std::cerr << "Lower conductivity bound reached" << std::endl;
						std::cerr << "\t li=" << li << "\tdm=" << dm[pindex] << "\tm=" << m_old[pindex] << "\tm+dm=" << m_new[pindex] << std::endl;
						std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(m_new[pindex]) << std::endl;
					}
					dm[pindex] = lmin - m_old[pindex];
					if (Verbose) std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(dm[pindex] + m_old[pindex]) << std::endl;
				}
				else if (m_new[pindex] > lmax) {
					if (Verbose) {
						std::cerr << rec_it_str() << std::endl;
						std::cerr << "Upper conductivity bound reached" << std::endl;
						std::cerr << "\t li=" << li << "\tdm=" << dm[pindex] << "\tm=" << m_old[pindex] << "\tm+dm=" << m_new[pindex] << std::endl;
						std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(m_new[pindex]) << std::endl;
					}
					dm[pindex] = lmax - m_old[pindex];
					if (Verbose) std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(dm[pindex] + m_old[pindex]) << std::endl;
				}
			}
		}
		
		if (fdT.bound()) {
			for (size_t li = 0; li < nlayers - 1; li++) {
				const size_t pindex = tindex(li);
				const double lmin = std::log10(E.min.thickness[li]);
				const double lmax = std::log10(E.max.thickness[li]);				
				if (m_new[pindex] < lmin) {
					if (Verbose) {
						std::cerr << rec_it_str() << std::endl;
						std::cerr << "Lower thickness bound reached" << std::endl;
						std::cerr << "\t li=" << li << "\tdm=" << dm[pindex] << "\tm=" << m_old[pindex] << "\tm+dm=" << m_new[pindex] << std::endl;
						std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(m_new[pindex]) << std::endl;
					}
					dm[pindex] = lmin - m_old[pindex];
					if (Verbose) std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(dm[pindex] + m_old[pindex]) << std::endl;
				}
				else if (m_new[pindex] > lmax) {
					if (Verbose) {
						std::cerr << rec_it_str() << std::endl;
						std::cerr << "Upper thickness bound reached" << std::endl;
						std::cerr << "\t li=" << li << "\tdm=" << dm[pindex] << "\tm=" << m_old[pindex] << "\tm+dm=" << m_new[pindex] << std::endl;
						std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(m_new[pindex]) << std::endl;
					}
					dm[pindex] = lmax - m_old[pindex];
					if (Verbose) std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(dm[pindex] + m_old[pindex]) << std::endl;
				}				
			}
		}

		for (size_t i = 0; i < cTDEmGeometry::size(); i++) {
			const std::string ename = cTDEmGeometry::element_name(i);			
			const cInvertibleFieldDefinition& e = fdG.cref(ename);
			if (e.bound()) {
				const size_t pindex = e.index;
				const double emin = G.min[ename];
				const double emax = G.max[ename];
				if (m_new[pindex] < emin) {										
					if (Verbose) {
						std::cerr << rec_it_str() << std::endl;
						std::cerr << "Lower " << ename << " bound reached" << std::endl;
						std::cerr << "\tdm=" << dm[pindex] << "\tm=" << m_old[pindex] << "\tm+dm=" << m_new[pindex] << std::endl;
					}
					dm[pindex] = emin - m_old[pindex];
				}
				else if (m_new[pindex] > emax) {
					if (Verbose) {
						std::cerr << rec_it_str() << std::endl;
						std::cerr << "Upper " << ename << " bound reached" << std::endl;
						std::cerr << "\tdm=" << dm[pindex] << "\tm=" << m_old[pindex] << "\tm+dm=" << m_new[pindex] << std::endl;
					}
					dm[pindex] = emax - m_old[pindex];
				}
			}
		}
		
		return dm;
	}

	cEarth1D get_earth(const Vector& parameters)
	{
		cEarth1D e = E.ref;
		if (solve_conductivity()) {
			for (size_t li = 0; li < nlayers; li++) {
				e.conductivity[li] = pow10(parameters[cindex(li)]);
			}
		}

		if (solve_thickness()) {
			for (size_t li = 0; li < nlayers - 1; li++) {
				e.thickness[li] = pow10(parameters[tindex(li)]);
			}
		}
		return e;
	}

	cTDEmGeometry get_geometry(const Vector& parameters)
	{
		cTDEmGeometry g = G.input;
		for (int i = 0; i < cTDEmGeometry::size(); i++) {
			std::string gname = cTDEmGeometry::element_name(i);
			auto p = fdG.cref(gname);
			if (p.solve) {
				g[gname] = parameters[p.index];
			}
		}		
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

	void forwardmodel(const Vector& parameters, Vector& predicted) {
		Matrix dummy;	
		_forwardmodel_(parameters, predicted, dummy, false);		
	}

	void forwardmodel_and_jacobian(const Vector& parameters, Vector& predicted, Matrix& jacobian) {
		_forwardmodel_(parameters, predicted, jacobian, true);
	}

	void _forwardmodel_(const Vector& parameters, Vector& predicted, Matrix& jacobian, bool computederivatives)
	{
		Vector pred_all(nAllData);
		Matrix J_all;

		cEarth1D      e = get_earth(parameters);
		cTDEmGeometry g = get_geometry(parameters);
		for (size_t si = 0; si < nsystems; si++) {
			cTDEmSystemInfo& S = SV[si];
			cTDEmSystem& T = S.T;
			const size_t nw = T.NumberOfWindows;
			T.setconductivitythickness(e.conductivity, e.thickness);
			T.setgeometry(g);

			//Forwardmodel
			T.LEM.calculation_type = cLEM::CalculationType::FORWARDMODEL;
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
					pred_all[wi + S.xzIndex] = xzfm[wi];
					if (S.CompInfo[1].Use) pred_all[wi + S.CompInfo[1].dataindex] = yfm[wi];
				}
			}
			else {
				for (size_t wi = 0; wi < nw; wi++) {
					if (S.CompInfo[0].Use) pred_all[wi + S.CompInfo[0].dataindex] = xfm[wi];
					if (S.CompInfo[1].Use) pred_all[wi + S.CompInfo[1].dataindex] = yfm[wi];
					if (S.CompInfo[2].Use) pred_all[wi + S.CompInfo[2].dataindex] = zfm[wi];
				}				
			}
						
			if (computederivatives) {					
				if (computederivatives) {
					J_all.resize(nAllData, nParam);
					J_all.setZero();
				}
				std::vector<double> xdrv(nw);
				std::vector<double> ydrv(nw);
				std::vector<double> zdrv(nw);

				if (solve_conductivity()) {
					for (size_t li = 0; li < nlayers; li++) {
						const size_t& pindex = cindex(li);
						T.LEM.calculation_type = cLEM::CalculationType::CONDUCTIVITYDERIVATIVE;
						T.LEM.derivative_layer = li;
						T.setprimaryfields();
						T.setsecondaryfields();

						fillDerivativeVectors(S, xdrv, ydrv, zdrv);
						//multiply by natural log(10) as parameters are in logbase10 units
						double sf = log(10.0)*e.conductivity[li];
						xdrv *= sf; ydrv *= sf; zdrv *= sf;
						fillMatrixColumn(J_all, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
					}
				}

				if (solve_thickness()) {
					for (size_t li = 0; li < nlayers - 1; li++) {
						size_t pindex = tindex(li);
						T.LEM.calculation_type = cLEM::CalculationType::THICKNESSDERIVATIVE;
						T.LEM.derivative_layer = li;
						T.setprimaryfields();
						T.setsecondaryfields();
						fillDerivativeVectors(S, xdrv, ydrv, zdrv);
						//multiply by natural log(10) as parameters are in logbase10 units
						double sf = log(10.0)*e.thickness[li];
						xdrv *= sf; ydrv *= sf; zdrv *= sf;
						fillMatrixColumn(J_all, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
					}
				}
				
				if (FreeGeometry) {					
					
					if (solve_geometry("tx_height")) {
						const size_t pindex = fdG.cref("tx_height").index;						
						T.LEM.calculation_type = cLEM::CalculationType::HDERIVATIVE;
						T.LEM.derivative_layer = INT_MAX;
						T.setprimaryfields();
						T.setsecondaryfields();
						fillDerivativeVectors(S, xdrv, ydrv, zdrv);
						fillMatrixColumn(J_all, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
					}

					if (solve_geometry("txrx_dx")) {
						const size_t pindex = fdG.cref("txrx_dx").index;
						T.LEM.calculation_type = cLEM::CalculationType::XDERIVATIVE;
						T.LEM.derivative_layer = INT_MAX;
						T.setprimaryfields();
						T.setsecondaryfields();
						fillDerivativeVectors(S, xdrv, ydrv, zdrv);
						fillMatrixColumn(J_all, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
					}

					if (solve_geometry("txrx_dy")) {
						const size_t pindex = fdG.cref("txrx_dy").index;
						T.LEM.calculation_type = cLEM::CalculationType::YDERIVATIVE;
						T.LEM.derivative_layer = INT_MAX;
						T.setprimaryfields();
						T.setsecondaryfields();
						fillDerivativeVectors(S, xdrv, ydrv, zdrv);
						fillMatrixColumn(J_all, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
					}

					if (solve_geometry("txrx_dz")) {
						const size_t pindex = fdG.cref("txrx_dz").index;
						T.LEM.calculation_type = cLEM::CalculationType::ZDERIVATIVE;
						T.LEM.derivative_layer = INT_MAX;
						T.setprimaryfields();
						T.setsecondaryfields();
						fillDerivativeVectors(S, xdrv, ydrv, zdrv);
						fillMatrixColumn(J_all, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
					}

					if (solve_geometry("rx_pitch")) {
						const size_t pindex = fdG.cref("rx_pitch").index;
						T.drx_pitch(xfm, zfm, g.rx_pitch, xdrv, zdrv);
						ydrv *= 0.0;
						fillMatrixColumn(J_all, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
					}

					if (solve_geometry("rx_roll")) {
						const size_t pindex = fdG.cref("rx_roll").index;
						T.drx_roll(yfm, zfm, g.rx_roll, ydrv, zdrv);
						xdrv *= 0.0;
						fillMatrixColumn(J_all, S, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
					}
				}				
			}
		}
		predicted = cull(pred_all);
		if(computederivatives) jacobian = cull(J_all);

		if (Verbose) {
			std::cerr << "\n-----------------\n";
			std::cerr << "It " << CIS.iteration + 1 << std::endl;
			std::cerr << J_all;
			std::cerr << "\n-----------------\n";
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

	void fillMatrixColumn(Matrix& M, cTDEmSystemInfo& S, const size_t& pindex, const std::vector<double>& xfm, const std::vector<double>& yfm, const std::vector<double>& zfm, const std::vector<double>& xzfm, const std::vector<double>& xdrv, const std::vector<double>& ydrv, const std::vector<double>& zdrv)
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

	void save_iteration_file(const cIterationState& S) {
		std::ofstream ofs(dumppath() + "iteration.dat");
		ofs << "Iteration "  << S.iteration << std::endl;
		ofs << "Lambda " << S.lambda << std::endl;
		ofs << "TargetPhiD " << S.targetphid << std::endl;
		ofs << "PhiD " << S.phid << std::endl;		
		ofs << "PhiM " << S.phim << std::endl;
		ofs << "PhiC " << S.phic << std::endl;
		ofs << "PhiT " << S.phit << std::endl;
		ofs << "PhiG " << S.phig << std::endl;
		ofs << "PhiS " << S.phis << std::endl;
		ofs << "PhiQ " << S.phiq << std::endl;
	};
	
	void writeresult(const int& pointindex, const cIterationState& S)
	{		
		const int& pi = pointindex;
		OM->begin_point_output();
		
		//Ancillary	
		OM->writefield(pi, Id.uniqueid, "uniqueid", "Inversion sequence number", UNITLESS, 1, NC_UINT, DN_NONE, 'I', 12, 0);
		for (size_t i = 0; i<AncFld.size(); i++) {
			cFdVrnt& fdv = AncFld[i].second;
			cAsciiColumnField c;
			std::string fname = fdv.fd.varname;
			IM->get_acsiicolumnfield(fname, c);			
			OM->writevrnt(pi, fdv.vnt, c);
		}

		//Geometry Input
		bool invertedfieldsonly = false;
		for (size_t i = 0; i < G.input.size(); i++) {
			if (invertedfieldsonly && solvegeometryindex(i) == false)continue;
			OM->writefield(pi, G.input[i], "input_" + G.input.element_name(i), "Input " + G.input.description(i), G.input.units(i), 1, NC_FLOAT, DN_NONE, 'F', 9, 2);
		}

		//Geometry Modelled		
		invertedfieldsonly = true;
		for (size_t i = 0; i < G.invmodel.size(); i++) {
			if (invertedfieldsonly && solvegeometryindex(i) == false)continue;
			OM->writefield(pi, G.invmodel[i], "inverted_" + G.invmodel.element_name(i), "Inverted " + G.invmodel.description(i), G.invmodel.units(i), 1, NC_FLOAT, DN_NONE, 'F', 9, 2);
		}
				
		//ndata
		OM->writefield(pi,
			nData, "ndata", "Number of data in inversion", UNITLESS,
			1, NC_UINT, DN_NONE, 'I', 4, 0);

		//Earth	
		OM->writefield(pi,
			nlayers,"nlayers","Number of layers ", UNITLESS,
			1, NC_UINT, DN_NONE, 'I', 4, 0);
		
		OM->writefield(pi,
			E.invmodel.conductivity, "conductivity", "Layer conductivity", "S/m",
			E.invmodel.conductivity.size(), NC_FLOAT, DN_LAYER, 'E', 15, 6);
		
		double bottomlayerthickness = 100.0;
		if (solve_thickness() == false && nlayers > 1) {
			bottomlayerthickness = E.invmodel.thickness[nlayers - 2];
		}
		std::vector<double> thickness = E.invmodel.thickness;
		thickness.push_back(bottomlayerthickness);

		OM->writefield(pi,
			thickness, "thickness", "Layer thickness", "m",
			thickness.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
					
				
		if (OO.PositiveLayerTopDepths) {			
			std::vector<double> dtop = E.invmodel.layer_top_depth();
			OM->writefield(pi,
				dtop, "depth_top", "Depth to top of layer", "m",
				dtop.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
		}

		if (OO.NegativeLayerTopDepths) {
			std::vector<double> ndtop = -1.0*E.invmodel.layer_top_depth();
			OM->writefield(pi,
				ndtop, "depth_top_negative", "Negative of depth to top of layer", "m",
				ndtop.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
		}
		
		if (OO.PositiveLayerBottomDepths) {
			std::vector<double> dbot = E.invmodel.layer_bottom_depth();
			OM->writefield(pi,
				dbot, "depth_bottom", "Depth to bottom of layer", "m",
				dbot.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
		}

		if (OO.NegativeLayerBottomDepths) {
			std::vector<double> ndbot = -1.0 * E.invmodel.layer_bottom_depth();
			OM->writefield(pi,
				ndbot, "depth_bottom_negative", "Negative of depth to bottom of layer", "m",
				ndbot.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
		}

		if (OO.InterfaceElevations) {			
			std::vector<double> etop = E.invmodel.layer_top_depth();
			etop += Id.elevation;
			OM->writefield(pi,
				etop, "elevation_interface", "Elevation of interface", "m",
				etop.size(), NC_FLOAT, DN_LAYER, 'F', 9, 2);
		}
				
		if (OO.ParameterSensitivity) {
			std::vector<double> ps = copy(ParameterSensitivity);
			if (solve_conductivity()) {
				std::vector<double> v(ps.begin() + cindex(0), ps.begin() + cindex(0) + nlayers);
				OM->writefield(pi,
					v, "conductivity_sensitivity", "Conductivity parameter sensitivity", UNITLESS,
					v.size(), NC_FLOAT, DN_LAYER, 'E', 15, 6);
			}
			
			if (solve_thickness()) {
				std::vector<double> v(ps.begin() + tindex(0), ps.begin() + tindex(0) + nlayers-1);
				v.push_back(0.0);//halfspace layer not a parameter
				OM->writefield(pi,
					v, "thickness_sensitivity", "Thickness parameter sensitivity", UNITLESS,
					v.size(), NC_FLOAT, DN_LAYER, 'E', 15, 6);
			}

			size_t k = 0;
			for (size_t gi = 0; gi < G.input.size(); gi++) {				
				if (solvegeometryindex(gi) == true) {
					std::string name = "inverted_" + G.input.element_name(gi) + "_sensitivity";
					std::string desc = G.input.description(gi) + " parameter sensitivity";
					OM->writefield(pi,
						ps[gindex(k)], name, desc, UNITLESS,
						1, NC_FLOAT, DN_NONE, 'E', 15, 6);
					k++;
				}
			}
		}

		if (OO.ParameterUncertainty) {
			std::vector<double> pu = copy(ParameterUncertainty);
			if (solve_conductivity()) {
				std::vector<double> v(pu.begin() + cindex(0), pu.begin() + cindex(0) + nlayers);
				OM->writefield(pi,
					v, "conductivity_uncertainty", "Conductivity parameter uncertainty", "log10(S/m)",
					v.size(), NC_FLOAT, DN_LAYER, 'E', 15, 6);
			}

			if (solve_thickness()) {
				std::vector<double> v(pu.begin() + tindex(0), pu.begin() + tindex(0) + nlayers - 1);
				v.push_back(0.0);//halfspace layer not a parameter
				OM->writefield(pi,
					v, "thickness_uncertainty", "Thickness parameter uncertainty", "log10(m)",
					v.size(), NC_FLOAT, DN_LAYER, 'E', 15, 6);
			}

			size_t k = 0;
			for (size_t gi = 0; gi < G.input.size(); gi++) {
				if (solvegeometryindex(gi) == false) continue;
				std::string name = "inverted_" + G.input.element_name(gi) + "_uncertainty";
				std::string desc = G.input.description(gi) + " parameter uncertainty";
				OM->writefield(pi,
					pu[gindex(k)], name, desc, G.input.units(gi),
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
		OM->writefield(pi, AlphaC, "AlphaC", "AlphaConductivity inversion parameter", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, AlphaT, "AlphaT", "AlphaThickness inversion parameter", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, AlphaG, "AlphaG", "AlphaGeometry inversion parameter", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, AlphaS, "AlphaS", "AlphaSmoothness inversion parameter", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);		
		OM->writefield(pi, AlphaQ, "AlphaQ", "AlphaHomogeneous inversion parameter", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, S.phid, "PhiD", "Normalised data misfit", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, S.phim, "PhiM", "Combined model norm", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, S.phic, "PhiC", "Conductivity reference model norm", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, S.phit, "PhiT", "Thickness reference model norm", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, S.phig, "PhiG", "Geometry reference model norm", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, S.phis, "PhiS", "Smoothness model norm", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, S.phiq, "PhiQ", "Homogeneity model norm", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, S.lambda, "Lambda", "Lambda regularization parameter", UNITLESS, 1, NC_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, S.iteration, "Iterations", "Number of iterations", UNITLESS, 1, NC_UINT, DN_NONE, 'I', 4, 0);
				
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

	bool solvegeometryindex(const size_t index) {		
		//eGeometryElementType getype = cTDEmGeometry::elementtype(index);		
		return fdG.cref(cTDEmGeometry::element_name(index)).solve;
	}
	
	bool readgeometry(cIFDMap& map)
	{
		bool status = true;
		for (size_t i = 0; i < cTDEmGeometry::size(); i++) {
			std::string ename = cTDEmGeometry::element_name(i);
			const cInvertibleFieldDefinition e = map.cref(ename);
			bool inpstatus = IM->read(e.input, G.input[i]);
			bool refstatus = IM->read(e.ref, G.ref[i]);

			if (refstatus == false && inpstatus == true) {
				G.ref[i] = G.input[i];
				refstatus = true;
			}
			else if (inpstatus == false && refstatus == true) {
				G.input[i] = G.ref[i];
				inpstatus = true;
			}

			if (inpstatus == false) {
				std::ostringstream msg;
				msg << "Error: no 'Input or Ref' defined for " << ename << std::endl;
				glog.errormsg(msg.str());
			}

			if (refstatus == false) {
				std::ostringstream msg;
				msg << "Error: no 'Ref or Input' defined for " << ename << std::endl;
				glog.errormsg(msg.str());
			}

			bool tfrstatus = IM->read(e.tfr, G.tfr[i]);
			if (tfrstatus == false) {
				G.tfr[i] = G.input[i];
			}

			if (e.solve) {
				bool stdstatus = IM->read(e.std, G.std[i]);
				if (stdstatus == false) {			
					std::ostringstream msg;
					msg << "Error: no 'Std' defined for "<< ename << std::endl;					
					glog.errormsg(msg.str());
				}

				bool minstatus = IM->read(e.min, G.min[i]);
				bool maxstatus = IM->read(e.max, G.max[i]);
			}
		}
		return status;
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
		
	void dump_first_iteration() {		
			const std::string dp = dumppath();
			makedirectorydeep(dumppath());

			write(Obs, dp + "observed.dat");
			write(Err, dp + "observed_std.dat");

			G.ref.write(dp + "geometry_start.dat");
			E.ref.write(dp + "earth_start.dat");

			G.ref.write(dp + "geometry_ref.dat");
			E.ref.write(dp + "earth_ref.dat");

			G.std.write(dp + "geometry_std.dat");
			E.std.write(dp + "earth_std.dat");
									
			std::ofstream ofs(dp+"Id.dat");
			char sep = '\n';			
			
			ofs << Id.uniqueid << sep;
			ofs << Id.survey << sep;
			ofs << Id.date << sep;
			ofs << Id.flight << sep;
			ofs << Id.line << sep;
			ofs << Id.fiducial << sep;
			ofs << Id.x << sep;
			ofs << Id.y << sep;
			ofs << Id.elevation << sep;						
	}

	void dump_iteration(const cIterationState& state) {
		const std::string dp = dumppath();
		makedirectorydeep(dp);
		writetofile(Obs, dp + "d.dat");
		writetofile(Err, dp + "e.dat");
		writetofile(state.pred, dp + "g.dat");
		cEarth1D e = get_earth(state.param);
		cTDEmGeometry g = get_geometry(state.param);
		e.write(dumppath() + "earth_inv.dat");
		g.write(dumppath() + "geometry_inv.dat");
		save_iteration_file(state);
	}

	bool initialise_sample() {	
		CIS = cIterationState();
		bool status = initialise_data();
		if (status == false) return false;
		initialise_parameters();
		initialise_Wd();
		initialise_Wm();		
		dump_W_matrices();
		return true;
	}

	void iterate() {
		_GSTITEM_				
		CIS.iteration = 0;
		CIS.lambda = 1e8;
		CIS.param = RefParam;		
		forwardmodel(CIS.param, CIS.pred);
		CIS.phid = phiData(CIS.pred);
		CIS.targetphid = CIS.phid;
		CIS.phim = phiModel(CIS.param, CIS.phic, CIS.phit, CIS.phig, CIS.phis, CIS.phim);

		TerminationReason = "Has not terminated";

		if (OO.Dump) {
			dump_first_iteration();
			dump_iteration(CIS);
		}
					
		double percentchange = 100.0;
		bool   keepiterating = true;
		while (keepiterating == true) {
			if (CIS.iteration >= MaxIterations) {
				keepiterating = false;
				TerminationReason = "Too many iterations";
			}
			else if (CIS.phid <= MinimumPhiD) {
				keepiterating = false;
				TerminationReason = "Reached minimum";
			}
			else if (percentchange < MinimumImprovement) {
				keepiterating = false;
				TerminationReason = "Small % improvement";
			}
			else {				
				if (CIS.iteration+1 >= BeginGeometrySolveIteration) FreeGeometry = true;
				else FreeGeometry = false;
				//if ((CIS.iteration+1)%2) FreeGeometry = false;
				//else FreeGeometry = true;
				
				Vector g;
				forwardmodel_and_jacobian(CIS.param, g, J);
				
				double targetphid = std::max(CIS.phid*0.7, MinimumPhiD);
				cTrial t  = targetsearch(CIS.lambda, targetphid);
				Vector dm = parameter_change(t.lambda, CIS.param, CIS.pred);
				Vector m = CIS.param + (t.stepfactor * dm);
				
				forwardmodel(m,g);
				double phid = phiData(g);

				percentchange = 100.0 * (CIS.phid - phid) / (CIS.phid);
				if (phid < CIS.phid) {				
					CIS.iteration++;
					CIS.param = m;
					CIS.pred = g;
					CIS.targetphid = targetphid;										
					CIS.phid   = phid;
					CIS.lambda = t.lambda;
					CIS.phim = phiModel(CIS.param, CIS.phic, CIS.phit, CIS.phig, CIS.phis, CIS.phiq);
					if (OO.Dump) dump_iteration(CIS);
				}						
			}			
		} 
		
		E.invmodel = get_earth(CIS.param);
		G.invmodel = get_geometry(CIS.param);
		forwardmodel_and_jacobian(CIS.param, CIS.pred, J);
		set_predicted();		
		ParameterSensitivity = compute_parameter_sensitivity();
		ParameterUncertainty = compute_parameter_uncertainty();
	}

	bool invert() {
		_GSTITEM_
		OutputMessage = "";
		if (read_record() == false) {
			OutputMessage += ", Skipping - could not parse record";
			return false;
		}

		if (initialise_sample() == false) {
			return false;
		}

		iterate();
		return true;
	}

	std::string record_id() {
		std::ostringstream s;			
		s << "Rec " << ixd(6) << 1 + IM->record();
		s << " Fl " << ixd(3) << Id.flight;
		s << " Ln " << ixd(7) << Id.line;
		s << " Fd " << fxd(10,2) << Id.fiducial;
		return s.str();
	}

	std::string record_result(const double& etime) {
		std::ostringstream s;						
		s << " Its="  << ixd(3) << CIS.iteration;
		s << " Phid=" << fxd(6,2) << CIS.phid;		
		s << " Time=" << fxd(4,1) << etime;
		s << " " << TerminationReason;
		s << " " << OutputMessage;		
		return s.str();
	}

	int execute() {
		_GSTITEM_		
		bool readstatus = true;
		int paralleljob = 0;		
		do{	
			int record = paralleljob*(int)IM->subsamplerate();			
			if ((paralleljob % Size) == Rank) {								
				readstatus = IM->read_record(record);				
				if (readstatus) {
					bool valid = IM->is_record_valid();
					if (valid == true) {
						if (OO.Dump) {
							dump_record_number();
						}

						double t1 = gettime();
						bool invstatus = invert();
						double t2 = gettime();
						double etime = t2 - t1;

						std::ostringstream s;
						if (invstatus) {							
							writeresult(record, CIS);							
							s << record_id();
							s << record_result(etime);
							s << std::endl;
							glog.logmsg(s.str());
							if (OutputMessage.size() > 0) {
								std::cerr << s.str();
							}
						}
						else {
							s << record_id();
							s << "Skipping: " << OutputMessage;
							s << std::endl;
							glog.logmsg(s.str());
							std::cerr << s.str();
						}
						
					}
				}
			}	
			paralleljob++;
		} while (readstatus == true);
		glog.close();
		return 0;
	}	

	double phiModel(const Vector& p)
	{
		double phic, phit, phig, phis, phiq;
		return phiModel(p, phic, phit, phig, phis, phiq);
	}

	double phiModel(const Vector& p, double& phic, double& phit, double& phig, double& phis, double& phiq)
	{
		phic = phiC(p);
		phit = phiT(p);
		phig = phiG(p);
		phis = phiS(p);
		phiq = phiQ(p);

		double v = phic + phit + phig + phis + phiq;
		return v;
	}

	double phiC(const Vector& p)
	{
		if (AlphaC == 0.0) return 0.0;
		Vector v = p - RefParam;
		return mtDm(v, Wc);
	}

	double phiT(const Vector& p)
	{
		if (AlphaT == 0.0)return 0.0;
		if (solve_thickness() == false)return 0.0;
		Vector v = p - RefParam;
		return mtDm(v, Wt);
	}

	double phiG(const Vector& p)
	{
		if (AlphaG == 0.0)return 0.0;
		Vector v = p - RefParam;
		return mtDm(v, Wg);
	}

	double phiS(const Vector& p)
	{
		if (AlphaS == 0)return 0.0;
		else return mtAm(p, Ws);
	}

	double phiQ(const Vector& p)
	{
		if (AlphaQ == 0)return 0.0;
		else return mtAm(p, Wq);
	}

	Vector solve_linear_system(const double& lambda, const Vector& param, const Vector& pred)
	{
		// Phi = (d-g(m)+Jm) Wd (d-g(m)+Jm) + lambda ( (m-m0)' Wr (m-m0) + m' Ws m) )
		//Ax = b
		//A = [J'WdJ + lambda (Wr + Ws)]
		//x = m(n+1)
		//b = J'Wd(d - g(m) + Jm) + lambda*Wr*m0
		//dm = m(n+1) - m = x - m

		
		const Vector& m = param;
		const Vector& g = pred;
		const Vector& d = Obs;		
		const Vector& e = Err;
		const Vector& m0 = RefParam;

		Matrix V = Wd;
		if (NormType == eNormType::L1) {
			for (size_t i = 0; i < nData; i++) {
				const double r = (d[i] - g[i]) / e[i];
				V(i, i) *= 1.0 / std::abs(r);
			}
		}

		Matrix JtV = J.transpose() * V;
		Matrix JtVJ = JtV * J;

		Vector b = JtV * (d - g + J * m) + lambda * (Wr * m0);
		Matrix A = JtVJ + lambda * Wm;
		Vector x = pseudoInverse(A) * b;
		return x;
	}
};

#endif
