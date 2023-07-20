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
#include "string_utils.h"
#include "vector_utils.h"

#include "airborne_types.h"
#include "cinverter.h"
#include "tdemsystem.h"
#include "tdemsysteminfo.h"
#include "samplebunch.h"
#include <Eigen/Cholesky>
#include <Eigen/LU>

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
		std::string p = DumpBasePath;
		p += strprint("rec_%07d", (int)datafilerecord + 1) + pathseparatorstring();
		p += strprint("it_%03d", (int)iteration) + pathseparatorstring();
		return p;
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

class cConstraint {

private:	
	bool refmodeldiff = false;
	std::vector<std::string> allowed_methods;

public:
	bool alreadyparsed = false;
	Matrix W;			
	double alpha = 0.0;
	std::string key;
	std::string method;
	std::string initials;
	std::string description;
	
	std::string alpha_field_name() const {
		return "Alpha_" + initials;
	};

	std::string phi_field_name() const {
		return "Phi_" + initials;
	};

	std::string alpha_field_description() const {
		return description + " constraint alpha parameter";
	};

	std::string phi_field_description() const {
		return description + " constraint model norm";
	};

	std::string matrix_name() const {
		return "W_" + initials;
	};

		
	cConstraint() {};

	cConstraint(
		const std::string& _key,		
		const std::vector<std::string>& _allowed_methods,
		const std::string& _initials,
		const std::string& _field_description		
	) {
		key = _key;		
		allowed_methods = _allowed_methods;
		initials = _initials;
		description = _field_description;		
	};

	static std::string get_key(const std::string& controlstring) {		
		std::istringstream is(controlstring);
		std::string k;
		is >> k;		
		return k;
	}

	bool verify_method() {		
		if (allowed_methods.size() == 0 && method.size() == 0) return true;

		for (size_t i = 0; i < allowed_methods.size(); i++) {
			if (strcasecmp(method, allowed_methods[i]) == 0) {
				return true;
			}
		}
		glog.errormsg(_SRC_, "Unknown linear constraint method %s for %s\n", method.c_str(), key.c_str());
		return false;
	}

	const bool& operates_on_difference_from_reference_model() const {
		return refmodeldiff;
	}

	void set_operates_on_difference_from_reference_model() {
		refmodeldiff = true;
	}

	double phi(const Vector& m, const Vector& m0) const {
		if (alpha == 0) return 0.0;
		if (operates_on_difference_from_reference_model()) {
			Vector p = m - m0;
			return mtAm(p, W);
		}
		else return mtAm(m, W);
	}

	void write_W_matrix(const std::string& dp) {
		writetofile(W, dp + matrix_name() + ".dat");
	}
};

class cLinearConstraint : public cConstraint {

private:
	

public:
	
	cLinearConstraint() {};
	
	cLinearConstraint(const std::string& _key, const std::vector<std::string>& _allowed_methods, const std::string& _initials, const std::string& _description)
		: cConstraint(_key, _allowed_methods, _initials, _description) {};

	bool parse(const std::string& controlstring) {
		
		std::istringstream is(controlstring);
		std::string inputkey;
		is >> inputkey;
		if (key != inputkey) return false;

		if (alreadyparsed == true) {
			std::ostringstream msg;
			msg << "Constraint has already been set: " << std::endl << controlstring << std::endl;
			glog.errormsg(msg.str());
		};

		is >> alpha;
		is >> method;
		verify_method();
		alreadyparsed = true;
		return true;
	}
};


class cNonLinearConstraint : public cConstraint {

private:
		

public:
	double _sd_;//Holding value for err std from control file

	Matrix J;//Jacobian of non-linear constraint
	Vector data;//Data of non-linear constraint
	Vector err;//Data error std of non-linear constraint
					
	cNonLinearConstraint() {};

	cNonLinearConstraint(const std::string& _key, const std::vector<std::string>& _allowed_methods, const std::string& _subscript, const std::string& _description)
		: cConstraint(_key, _allowed_methods, _subscript, _description) {};

	bool parse(const std::string& controlstring) {
		
		std::istringstream is(controlstring);
		std::string inputkey;
		is >> inputkey;
		if (key != inputkey) return false;
		
		if (alreadyparsed == true) {
			std::ostringstream msg;
			msg << "Constraint has already been set: " << std::endl << controlstring << std::endl;
			glog.errormsg(msg.str());
		};

		is >> alpha;
		is >> method;		
		is >> _sd_;		
		alreadyparsed = true;
		return true;
	}


	double phi(const Vector& predicted) const {
		if (alpha == 0) return 0.0;
		Vector delta = data - predicted;
		return mtAm(delta, W);
	}
	
};

class cSBSInverter : public cInverter {

	double ErrorAddition = 0.0;
	using cIFDMap = cKeyVec<std::string, cInvertibleFieldDefinition, caseinsensetiveequal<std::string>>;
	const size_t XCOMP = 0;
	const size_t YCOMP = 1;
	const size_t ZCOMP = 2;
	const size_t XZAMP = 3;
	std::vector<std::vector<std::vector<std::vector<int>>>> _dindex_;

	int    BeginGeometrySolveIteration = 0;
	bool   FreeGeometry = false;	
	Matrix Wr;//Composite reference model matrix

	size_t StartRecord = 0; // 1-based first record to be inverted
	size_t EndRecord = INT_MAX;// 1-based last record to be inverted
	size_t Subsample = 0;
	size_t nSoundings = 0;
	size_t nBunchSubsample = 0;	
	size_t nDataPerSounding = 0;
	size_t nAllData = 0;			
	size_t nLayers = 0;
	size_t nParamPerSounding = 0;
	size_t nGeomParamPerSounding = 0;
	size_t nScalingParam = 0;
	size_t cOffset = 0;//Offset within sample of conductivity parameters
	size_t tOffset = 0;//Offset within sample of thickness parameters
	//size_t gOffset = 0;//Offset within sample of geometry parameters
	//size_t sOffset = 0;//Offset within whole inversion

	size_t nSystems = 0;
	size_t pointsoutput = 0;
	std::vector<cGeomStruct> G;
	std::vector<cEarthStruct> E;		
	cOutputOptions OO;
	std::vector<cTDEmSystemInfo> SV;

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
	};
	std::vector<SampleId> Id;
	std::vector<cKeyVec<std::string, cFdVrnt, caseinsensetiveequal<std::string>>> AncFld;
	
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
		Vector vcull(nData);
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
		
	const size_t& nsoundings() {
		return nSoundings;
	}

	cLinearConstraint LCrefc;//Conductivity reference model constraint
	cLinearConstraint LCreft;//Thickness reference model constraint
	cLinearConstraint LCrefg;//Geometry reference model constraint
	cLinearConstraint LCrefs;//Scaling factor reference model constraint

	cLinearConstraint LCvcsmth;//Vertical conductivity smoothness constraint
	cLinearConstraint LCvcsim;//Vertical Conductivity Similarity

	cLinearConstraint LClatc;//Lateral conductivity constraint	
	cLinearConstraint LClatg;//Lateral geometry constraint

	cNonLinearConstraint NLCcablen;//Cable length constraint

	cSBSInverter(const std::string& controlfile, const int& size, const int& rank, const bool& usingopenmp, const std::string commandline) 
					: cInverter(controlfile, size, rank, usingopenmp, commandline)
	{		
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
		//std::cout << "Destroying cSBSInverter\n";
	}
			
	void loadcontrolfile(const std::string& filename)
	{		
		glog.logmsg(0, "Loading control file %s\n", filename.c_str());
		Control = cBlock(filename);
		cBlock ob = Control.findblock("Output");
		cBlock ib = Control.findblock("Input");

		OO = cOutputOptions(ob);
		Verbose = ob.getboolvalue("verbose");

		
		if (Rank == 0) {
			std::string od = extractfiledirectory(OO.LogFile);
			makedirectorydeep(od);			
		}
		#if defined _MPI_ENABLED
			cMpiEnv::world_barrier();
		#endif
		
		std::string suffix = stringvalue(Rank, ".%04d");
		OO.LogFile = insert_after_filename(OO.LogFile, suffix);
		openlogfile(); //load this first to get outputlogfile opened

		//Load control file
		parse_options();
		
		//Setup InputManager
		if (cInputManager::isnetcdf(ib)) {
			#if defined HAVE_NETCDF
				IM = std::make_unique<cNetCDFInputManager>(ib);
			#else
				glog.errormsg(_SRC_, "Sorry NETCDF I/O is not available in this executable\n");
			#endif			
			//std::string s = IM->datafilename();
		}
		else {
			IM = std::make_unique<cASCIIInputManager>(ib);
		}
		IM->set_subsample_rate(Subsample);
		
		initialise_systems();		

		//Setup OutputManager
		if (cOutputManager::isnetcdf(ob)) {
			#if defined HAVE_NETCDF
				OM = std::make_unique<cNetCDFOutputManager>(ob, Size, Rank);
			#else
				glog.errormsg(_SRC_, "Sorry NETCDF I/O is not available in this executable\n");
			#endif						
		}
		else {
			OM = std::make_unique<cASCIIOutputManager>(ob, Size, Rank);
		}

		if (Rank == 0) {
			std::string od = extractfiledirectory(OM->datafilename());
			makedirectorydeep(od);
		}
		#if defined _MPI_ENABLED
			cMpiEnv::world_barrier();
		#endif

		OM->opendatafile(IM->datafilename(), IM->subsamplerate());
	}

	bool solve_thickness() const
	{
		return fdT.solve;
	};
	
	bool solve_conductivity() const {
		return fdC.solve;
	};
	
	bool solve_geometry_element(const std::string& gname) const {								
		return fdG.cref(gname).solve;		
	};

	bool solve_geometry_index(const size_t index) const {		
		return fdG.cref(cTDEmGeometry::element_name(index)).solve;
	}
	
	bool solve_geometry() const {		
		if (nGeomParamPerSounding > 0) return true;
		return true;
	};

	bool solve_scalingfactors(const size_t sysi, const size_t ci) {
		return SV[sysi].CompInfo[ci].fdSF.solve;
	}

	bool solve_scalingfactors() {
		if (nScalingParam > 0) return true;
		else return false;
	}

	std::string bunch_id() {
		const size_t si = Bunch.master_index();
		const size_t& record = Bunch.master_record();
		std::ostringstream s;
		s << "Rec " << ixd(6) << 1 + record;
		s << " Flt " << ixd(3) << Id[si].flight;
		s << " Line " << ixd(7) << Id[si].line;
		s << " Fid " << fxd(10, 2) << Id[si].fiducial;
		return s.str();
	}

	std::string bunch_result(const double& etime) {
		std::ostringstream s;
		s << " Its=" << ixd(3) << CIS.iteration;
		s << " Phid=" << fxd(6, 2) << CIS.phid;
		s << " Time=" << fxd(4, 1) << etime;
		s << " " << TerminationReason;
		s << " " << OutputMessage;
		//s << " nF= " << nForwards / CIS.iteration;
		//s << " nJ= " << nJacobians;
		return s.str();
	}

	std::string dumppath() const
	{
		const size_t& record  = Bunch.master_record();
		std::string s = OO.DumpPath(record, CIS.iteration);
		return s;
	};

	void dump_record_number() {
		const size_t& record = Bunch.master_record();
		std::ofstream of(dumppath() + "record.dat");
		of << "Record\t" << record << std::endl;
	}

	int cindex(const size_t& si, const size_t& li) {
		if (solve_conductivity() == false) {
			glog.errormsg("Out of boundes in cindex()\n");
		}
		return (int) (si*nParamPerSounding + cOffset + li);
	}

	int tindex(const size_t& si, const size_t& li) {
		if (solve_thickness() == false) {
			glog.errormsg("Out of boundes in tindex()\n");
		}		
		return (int) (si*nParamPerSounding + tOffset + li);
	}
	
	int gindex(const size_t& si, const std::string& gname) const {				
		cInvertibleFieldDefinition val;
		bool status = fdG.get(gname,val);
		if (status) {
			if (val.offset >= 0) return (int)(si * nParamPerSounding + val.offset);
		} 
		return -1;
	}

    int gindex(const size_t& si, const size_t& gi) const {		
		int goff = fdG[gi].second.offset;
		if (goff < 0) return -1;
		return (int)(si*nParamPerSounding + goff);
	}

	int sfindex(const size_t& sysi, const size_t& ci) const {		
		int pi = SV[sysi].CompInfo[ci].fdSF.offset;
		if (pi < 0) return -1;		
		return pi;
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
	
	void parse_options()
	{
		cBlock b = Control.findblock("Options");				
		if (b.getvalue("StartRecord", StartRecord) == false) {
			StartRecord = 1;
		}

		if (b.getvalue("EndRecord", EndRecord) == false) {
			EndRecord = INT_MAX;
		}

		if (b.getvalue("Subsample", Subsample) == false) {
			Subsample = 1;
		}
		
		if (b.getvalue("SoundingsPerBunch", nSoundings)==false) {
			nSoundings = 1;
		}

		if (b.getvalue("BunchSubsample", nBunchSubsample) == false) {
			nBunchSubsample = 1;
		}

		if (b.getvalue("ErrorAddition", ErrorAddition) == false) {
			ErrorAddition = 0.0;
		}
		
		cBlock cb = b.findblock("Constraints");
		parse_constraints(cb);
				
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
		
		MaxIterations = b.getsizetvalue("MaximumIterations");
		MinimumPhiD = b.getdoublevalue("MinimumPhiD");
		MinimumImprovement = b.getdoublevalue("MinimumPercentageImprovement");
	}
	
	void parse_constraints(const cBlock& b) {
		
		LCrefc = cLinearConstraint("ConductivityReferenceModel", { }, "RefCon", "Conductivity reference model");
		LCreft = cLinearConstraint("ThicknessReferenceModel", { }, "RefThk", "Thickness reference model");
		LCrefg = cLinearConstraint("GeometryReferenceModel", { }, "RefGeom", "Geometry reference model");
		LCrefs = cLinearConstraint("ScalingFactorsReferenceModel", { }, "RefScalingFactors", "ScalingFactors reference model");

		LCvcsmth   = cLinearConstraint("VerticalConductivity", {"Minimise1stDerivatives","Minimise2ndDerivatives" }, "VCsmth", "Vertical conductivity smoothness");
		LCvcsim = cLinearConstraint("VerticalConductivitySimilarity", { }, "VCsim", "Vertical conductivity similarity");
		LClatc   = cLinearConstraint("LateralConductivity", {"Minimise1stDerivatives", "Minimise2ndDerivatives","Similarity" }, "LatCon", "Lateral conductivity smoothness");		
		LClatg   = cLinearConstraint("LateralGeometry", { "Similarity","MinimiseAccelerations","MinimiseAccelerationDerivatives","Minimise2ndDerivativesOfDifferenceFromReferenceModel" }, "LatGeom", "Lateral geometry smoothness");		
		NLCcablen  = cNonLinearConstraint("CableLength", { "Input","InputBunchMean","BunchSimilarity" }, "CabLength", "Cable length");
		
		for (size_t i = 0; i < b.Entries.size(); i++) {			
			const std::string cstr = b.Entries[i];
			std::string key = cConstraint::get_key(cstr);
			if (key[0] == '/')continue;
			if (LCrefc.parse(cstr)){
				LCrefc.set_operates_on_difference_from_reference_model();
			}
			else if (LCreft.parse(cstr)) {
				LCreft.set_operates_on_difference_from_reference_model();
			}
			else if (LCrefg.parse(cstr)) {
				LCrefg.set_operates_on_difference_from_reference_model();
			}
			else if (LCrefs.parse(cstr)) {
				LCrefs.set_operates_on_difference_from_reference_model();
			}
			else if (LCvcsmth.parse(cstr)) {
				//nothing to do
			}
			else if (LCvcsim.parse(cstr)) {
				//nothing to do
			}
			else if (LClatc.parse(cstr)) {
				//nothing to do
			}
			else if (LClatg.parse(cstr)) {
				if (LClatg.method == "Similarity") {
					LClatg.set_operates_on_difference_from_reference_model();
				}
				if (LClatg.method == "Minimise2ndDerivativesOfDifferenceFromReferenceModel") {
					LClatg.set_operates_on_difference_from_reference_model();
				}
			}
			else if(NLCcablen.parse(cstr)) {				
				//nothing to do
			}
			else {
				std::stringstream msg;
				msg << "Unknown constraint " << key << std::endl;
				glog.errormsg(msg.str());
			}			
		}				
	}

	void set_field_definitions()
	{
		cBlock b = Control.findblock("Input.AncillaryFields");
		set_field_definitions_ancillary(b);

		if (nSoundings > 1) {
			if (AncFld[0].keyindex("line") < 0) {
				glog.errormsg("Must specify a linenumber field\n");
			}
		}

		b = Control.findblock("Input.Geometry");
		fdG = set_field_definitions_geometry(b);

		b = Control.findblock("Input.Earth");
		bool status = b.getvalue("NumberOfLayers", nLayers);
		if (status == false) {
			std::stringstream msg;
			msg << "The NumberOfLayers must be specified in Input.Earth\n";
			glog.errormsg(msg.str());
		}

		fdC = cInvertibleFieldDefinition(b, "Conductivity");
		if (nLayers > 1) {
			fdT = cInvertibleFieldDefinition(b, "Thickness");
		}
	}

	void set_field_definitions_ancillary(const cBlock& parent) {
		AncFld.resize(nSoundings);
		const cBlock& b = parent;
		for (size_t i = 0; i < b.Entries.size(); i++) {
			std::string key   = b.key(i);
			std::string value = b.value(i);
			if (key[0] == '/') continue;//Skip commented out fields
			cFieldDefinition fd(parent, key);						
			cFdVrnt fdvrnt(fd,cVrnt());			
			IM->set_variant_type(fd, fdvrnt.vnt);
			for (size_t si = 0; si < nSoundings; si++) {
				AncFld[si].add(key, fdvrnt);
			}
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
	
	void setup_parameters()
	{
		Id.resize(nSoundings);
		E.resize(nSoundings);
		G.resize(nSoundings);

		nParamPerSounding = 0;
		nGeomParamPerSounding = 0;
		cOffset = 0;
		tOffset = 0;
		//gOffset = 0;		
		//sOffset = 0;

		if (solve_conductivity()) {
			fdC.offset = 0;
			tOffset += nLayers;
			//gOffset += nLayers;
			nParamPerSounding += nLayers;
		}

		if (solve_thickness()) {
			fdT.offset = (int)tOffset;
			//gOffset += nLayers-1;
			nParamPerSounding += nLayers - 1;
		}

		//Geometry params			
		for (size_t gi = 0; gi < cTDEmGeometry::size(); gi++) {
			std::string gname = cTDEmGeometry::element_name(gi);
			cInvertibleFieldDefinition& g = fdG.ref(gname);
			if (g.solve) {
				g.offset = (int)nParamPerSounding;
				nGeomParamPerSounding++;
				nParamPerSounding++;
			}
			else {
				g.offset = -1;
			}
		}

		//Scaling params
		for (size_t sysi = 0; sysi < SV.size(); sysi++){
			for (size_t ci = 0; ci < 3; ci++) {
				if (solve_scalingfactors(sysi, ci)){
					SV[sysi].CompInfo[ci].fdSF.offset = (int)(nParamPerSounding * nSoundings + nScalingParam);
					nScalingParam++;
				}
			}
		}

		nParam = nParamPerSounding*nSoundings + nScalingParam;				
		RefParam.resize(nParam);
		RefParamStd.resize(nParam);

		if (nGeomParamPerSounding == 0) {
			LCrefg.alpha = 0.0;
			NLCcablen.alpha = 0.0;
		}

		if (nScalingParam == 0) {
			LCrefs.alpha = 0.0;			
		}
		
		if (nSoundings == 1) {
			LClatc.alpha = 0.0;
			LClatg.alpha = 0.0;
		}
	}
	
	void initialise_Wc() {
		cLinearConstraint& C = LCrefc;
		C.W = Matrix::Zero(nParam, nParam);
		if (solve_conductivity() == false)return;

		for (size_t si = 0; si < nSoundings; si++) {
			const cEarthStruct& e = E[si];
			std::vector<double> t(nLayers);
			if (nLayers == 1) {
				t[0] = 1;
			}
			else if (nLayers == 2) {
				t[0] = e.ref.thickness[0];
				t[1] = e.ref.thickness[0];
			}
			else {
				for (size_t i = 0; i < (nLayers - 1); i++) {
					t[i] = e.ref.thickness[i];
				}
				t[nLayers - 1] = (t[nLayers - 2] / t[nLayers - 3]) * t[nLayers - 2];
			}

			double tsum = 0.0;
			for (size_t li = 0; li < nLayers; li++)tsum += t[li];
			double tavg = tsum / (double)nLayers;

			double s = C.alpha / (double)(nLayers * nSoundings);
			for (size_t li = 0; li < nLayers; li++) {
				int p = cindex(si, li);
				C.W(p, p) = s * (t[li] / tavg) / (RefParamStd[p] * RefParamStd[p]);
			}
		}
	}

	void initialise_Wt() {
		cLinearConstraint& C = LCreft;
		C.W = Matrix::Zero(nParam, nParam);
		if (solve_thickness() == false)return;

		const double s = C.alpha / (double)((nLayers - 1) * nSoundings);
		for (size_t si = 0; si < nSoundings; si++) {
			for (size_t li = 0; li < nLayers - 1; li++) {
				const int pi = tindex(si, li);
				C.W(pi, pi) = s / (RefParamStd[pi] * RefParamStd[pi]);
			}
		}
	}

	void initialise_Wg() {
		cLinearConstraint& C = LCrefg;
		C.W = Matrix::Zero(nParam, nParam);
		if (nGeomParamPerSounding <= 0)return;

		double s = C.alpha / (double)(nGeomParamPerSounding * nSoundings);
		for (size_t si = 0; si < nSoundings; si++) {
			for (size_t gi = 0; gi < cTDEmGeometry::size(); gi++) {
				const int pi = gindex(si, gi);
				if (pi >= 0) {
					C.W(pi, pi) = s / (RefParamStd[pi] * RefParamStd[pi]);
				}
			}
		}
	}

	void initialise_Ws() {
		cLinearConstraint& C = LCrefs;
		C.W = Matrix::Zero(nParam, nParam);
		if (solve_scalingfactors() == false)return;

		double s = C.alpha / (double)(nScalingParam);
		for (size_t sysi = 0; sysi < SV.size(); sysi++) {
			for (size_t ci = 0; ci < 3; ci++) {
				const int pi = sfindex(sysi, ci);
				if (pi >= 0) {
					C.W(pi, pi) = s / (RefParamStd[pi] * RefParamStd[pi]);
				}
			}
		}
	}

	void initialise_VC() {
		cLinearConstraint& C = LCvcsmth;
		C.W = Matrix::Zero(nParam, nParam);
		if (solve_conductivity() == false) return;
		if (C.alpha == 0 ) return;
		if (C.method == "Minimise1stDerivatives") {			
			initialise_VC_1st_derivative(C);
		}
		else if (C.method == "Minimise2ndDerivatives") {
			initialise_VC_2nd_derivative(C);			
		}		
	}

	void initialise_VC_1st_derivative(cLinearConstraint& C)
	{		
		if (nLayers < 3) return;		
		Matrix L = Matrix::Zero(nSoundings * (nLayers - 1), nParam);
		size_t nrows = 0;
		for (size_t si = 0; si < nSoundings; si++) {
			const cEarthStruct& e = E[si];
			std::vector<double> t = e.ref.dummy_thickness();
			double tavg = mean(t);
			for (size_t li = 1; li < nLayers; li++) {
				const int pi0 = cindex(si, li - 1);
				const int pi1 = cindex(si, li);
				double t1 = t[li - 1];
				double t2 = t[li];
				double d12 = (t1 + t2) / 2.0;
				double s = std::sqrt(t2 / tavg);//sqrt because it gets squared in L'L		
				L(nrows, pi0) = -s / d12;
				L(nrows, pi1) = s / d12;
				nrows++;
			}
		}
		C.W = L.transpose() * L;
		C.W *= (C.alpha / (double)(nrows));
	}

	void initialise_VC_2nd_derivative(cLinearConstraint& C)
	{		
		if (nLayers < 3) return;		
		Matrix L = Matrix::Zero(nSoundings * nLayers, nParam);
		size_t nrows = 0;
		for (size_t si = 0; si < nSoundings; si++) {
			const cEarthStruct& e = E[si];
			std::vector<double> t = e.ref.dummy_thickness();
			double tavg = mean(t);
			for (size_t li = 1; li < nLayers - 1; li++) {
				const int pi0 = cindex(si, li - 1);
				const int pi1 = cindex(si, li);
				const int pi2 = cindex(si, li + 1);
				double t1 = t[li - 1];
				double t2 = t[li];
				double t3 = t[li + 1];
				double d12 = (t1 + t2) / 2.0;
				double d23 = (t2 + t3) / 2.0;
				double s = std::sqrt(t2 / tavg);//sqrt because it gets squared in L'L		
				L(nrows, pi0) = s / d12;
				L(nrows, pi1) = -s / d12 - s / d23;
				L(nrows, pi2) = s / d23;
				nrows++;
			}

			//Minimise 1st deriv on first interface
			int pi0 = cindex(si, 0);
			int pi1 = cindex(si, 1);
			double s = std::sqrt(t[0] / tavg);
			L(nrows, pi0) = -1.0 * s / t[0];
			L(nrows, pi1) =  1.0 * s / t[0];
			nrows++;

			//Minimise 1st deriv on last interface
			pi0 = cindex(si, nLayers-2);
			pi1 = cindex(si, nLayers-1);
			s = std::sqrt(t[nLayers-2] / tavg);
			L(nrows, pi0) = -1.0 * s / t[nLayers-1];
			L(nrows, pi1) =  1.0 * s / t[nLayers-1];
			nrows++;			
		}
		C.W = L.transpose() * L;
		C.W *= (C.alpha / (double)(nrows));
	}
	
	void initialise_LC() {
		cLinearConstraint& C = LClatc;
		C.W = Matrix::Zero(nParam, nParam);
		if (solve_conductivity() == false) return;
		if (C.alpha == 0) return;

		if (C.method == "Minimise1stDerivatives") {
			initialise_LC_1st_derivative(C);
		}
		else if (C.method == "Minimise2ndDerivatives") {
			initialise_LC_2nd_derivative(C);
		}
		else if (C.method == "Similarity") {
			initialise_LC_similarity(LClatc);
		}
	}

	void initialise_LC_1st_derivative(cLinearConstraint& C)
	{
		if (nSoundings < 2) return;
		Matrix L = Matrix::Zero((nSoundings - 1) * nLayers, nParam);
		size_t nrows = 0;
		for (size_t si = 1; si < nSoundings; si++) {
			double d = std::hypot(Id[si].x - Id[si - 1].x, Id[si].y - Id[si - 1].y);
			for (size_t li = 0; li < nLayers; li++) {
				const int pi0 = cindex(si - 1, li);
				const int pi1 = cindex(si, li);
				L(nrows, pi0) = 1.0 / d;
				L(nrows, pi1) = -1.0 / d;
				nrows++;
			}
		}
		C.W = L.transpose() * L;
		C.W *= (C.alpha / (double)(nrows));
	}

	void initialise_LC_2nd_derivative(cLinearConstraint& C)
	{
		if (nSoundings < 3) return;
		Matrix L = Matrix::Zero((nSoundings - 2) * nLayers, nParam);
		size_t nrows = 0;
		for (size_t si = 1; si < nSoundings - 1; si++) {
			double d01 = std::hypot(Id[si].x - Id[si - 1].x, Id[si].y - Id[si - 1].y);
			double d12 = std::hypot(Id[si].x - Id[si + 1].x, Id[si].y - Id[si + 1].y);
			for (size_t li = 0; li < nLayers; li++) {
				const int pi0 = cindex(si - 1, li);
				const int pi1 = cindex(si, li);
				const int pi2 = cindex(si + 1, li);
				L(nrows, pi0) = 1.0 / d01;
				L(nrows, pi1) = -1.0 / d01 - 1.0 / d12;
				L(nrows, pi2) = 1.0 / d12;
				nrows++;
			}
		}
		C.W = L.transpose() * L;
		C.W *= (C.alpha / (double)(nrows));
	}

	void initialise_LC_similarity(cLinearConstraint& C)
	{
		if (nSoundings < 2) return;
		Matrix L = Matrix::Zero(nSoundings * nLayers, nParam);
		size_t nrows = 0;
		for (size_t li = 0; li < nLayers; li++) {
			for (size_t si = 0; si < nSoundings; si++) {
				const int psi = cindex(si, li);
				double std = RefParamStd[psi];
				L(nrows, psi) = 1.0/std;
				for (size_t ri = 0; ri < nSoundings; ri++) {
					if (ri != si) {
						const int pri = cindex(ri, li);
						L(nrows, pri) = -1.0 / (std*(double)(nSoundings - 1));
					}
				}
				nrows++;
			}
		}
		C.W = L.transpose() * L;
		C.W *= C.alpha / (double)(nrows);
	}

	void initialise_QC(cLinearConstraint& C)
	{
		C.W = Matrix::Zero(nParam, nParam);
		if (C.alpha == 0) return;
		if (solve_conductivity() == false) return;
		Matrix L = Matrix::Zero(nLayers * nSoundings, nParam);

		size_t nrows = 0;
		for (size_t si = 0; si < nSoundings; si++) {
			const cEarthStruct& e = E[si];
			std::vector<double> t = e.ref.dummy_thickness();
			double tavg = mean(t);
			//Loop over constraints equations
			for (size_t li = 0; li < nLayers; li++) {
				//double s = std::sqrt(t[li] / tavg);//sqrt because it gets squared in L'L									
				const int lpindex = cindex(si, li);
				L(nrows, lpindex) = 1.0;
				//Loop over layers for this equation
				for (size_t ki = 0; ki < nLayers; ki++) {
					const int kpindex = cindex(si, ki);					
					if(li!=ki){
						L(nrows, kpindex) = -1.0 / ((double)nLayers - 1.0);
					}
				}
				nrows++;
			}
		}
		C.W = L.transpose() * L;
		C.W *= (C.alpha / (double)(nrows));
	}
	
	void initialise_LG() {
		cLinearConstraint& C = LClatg;
		C.W = Matrix::Zero(nParam, nParam);
		if (solve_geometry() == false) return;
		if (C.alpha == 0) return;

		if (C.method == "MinimiseAccelerations") {
			initialise_LG_accelerations(C);
		}
		else if (C.method == "MinimiseAccelerationDerivatives") {
			initialise_LG_acceleration_derivatives(C);
		}
		else if (C.method == "Similarity") {
			initialise_LG_similarity(C);
		}
		else if (C.method == "Minimise2ndDerivativesOfDifferenceFromReferenceModel") {
			initialise_LG_accelerations(C);
		}
		
	}

	void initialise_LG_accelerations(cLinearConstraint& C)
	{				
		if (nSoundings < 3) return;		
		Matrix L = Matrix::Zero((nSoundings - 2) * nGeomParamPerSounding, nParam);
		size_t nrows = 0;
		for (size_t gi = 0; gi < cTDEmGeometry::size(); gi++) {
			if (solve_geometry_index(gi) == false)continue;
			for (size_t si = 1; si < nSoundings - 1; si++) {				
				double d01 = std::hypot(Id[si].x - Id[si - 1].x, Id[si].y - Id[si - 1].y);
				double d12 = std::hypot(Id[si].x - Id[si + 1].x, Id[si].y - Id[si + 1].y);

				const int pi0 = gindex(si - 1, gi);
				const int pi1 = gindex(si, gi);
				const int pi2 = gindex(si + 1, gi);
				const double std = RefParamStd[pi1];
				L(nrows, pi0) = 1.0 / (d01*std);
				L(nrows, pi1) = -1.0 / (d01*std) - 1.0 / (d12*std);
				L(nrows, pi2) = 1.0 / (d12*std);
				nrows++;
			}
		}
		C.W = L.transpose() * L;
		C.W *= (C.alpha / (double)(nrows));
	}

	void initialise_LG_acceleration_derivatives(cLinearConstraint& C)
	{		
		if (nSoundings < 5) return;		
		Matrix L = Matrix::Zero(nSoundings * nGeomParamPerSounding, nParam);
		size_t nrows = 0;
		
		double d = 0.0;
		for (size_t j = 0; j < nSoundings-1; j++) {
			d += std::hypot(Id[j].x - Id[j + 1].x, Id[j].y - Id[j + 1].y);
		}
		d = d / (double)(nSoundings-1);//average sample distance

		for (size_t gi = 0; gi < cTDEmGeometry::size(); gi++) {			
			if (solve_geometry_index(gi) == false)continue;
			for (size_t si = 2; si < nSoundings - 2; si++) {																
				const int pi0 = gindex(si - 2, gi);
				const int pi1 = gindex(si - 1, gi);
				const int pi2 = gindex(si, gi);
				const int pi3 = gindex(si + 1, gi);
				const int pi4 = gindex(si + 2, gi);	
				const double std = RefParamStd[pi2];

				L(nrows, pi0) = -1.0 / (d * std);
				L(nrows, pi1) = 4.0 / (d * std);
				L(nrows, pi2) = -6.0 / (d * std);
				L(nrows, pi3) = 4.0 / (d * std);
				L(nrows, pi4) = -1.0 / (d * std);
				nrows++;
			}			
		}
		C.W = L.transpose() * L;
		C.W *= (C.alpha / (double)(nrows));
	}

	void initialise_LG_similarity(cLinearConstraint& C)
	{				
		if (nSoundings < 2) return;		
		Matrix L = Matrix::Zero(nSoundings * nGeomParamPerSounding, nParam);
		size_t nrows = 0;		
		for (size_t gi = 0; gi < cTDEmGeometry::size(); gi++) {
			if (solve_geometry_index(gi) == false)continue;
			for (size_t si = 0; si < nSoundings; si++) {
				const int spi = gindex(si, gi);				
				const double std = RefParamStd[spi];
				L(nrows, spi) = 1.0 / std;
				for (size_t ri = 0; ri < nSoundings; ri++) {
					const int rpi = gindex(ri, gi);					
					if (si != ri) {
						L(nrows, rpi) = -1.0 / (std * (double)(nSoundings-1));
					}
				}
				nrows++;
			}
		}		
		C.W = L.transpose() * L;
		C.W *= C.alpha / (double)(nrows);
	}

	void initialise_CableLengthConstraint() {
		cNonLinearConstraint& C = NLCcablen;
		C.W = Matrix::Zero(nSoundings, nSoundings);
		C.J = Matrix::Zero(nSoundings, nParam);
		C.err = Vector::Zero(nSoundings);
		C.data = Vector::Zero(nSoundings);
		if (C.alpha == 0.0) return;
		if (nGeomParamPerSounding <= 0) return;

		double s = C.alpha / (double)(nSoundings);
		C.data.resize(nSoundings);
		for (size_t si = 0; si < nSoundings; si++) {
			C.err[si] = C._sd_;

			if (NLCcablen.method == "Input") {
				C.data[si] = G[si].input.txrx_dr();
			}
			else if (NLCcablen.method == "InputBunchMean") {
				C.data[si] = G[si].input.txrx_dr();
			}
			else if (NLCcablen.method == "BunchSimilarity") {
				C.data[si] = 0.0;
			}
			else {
				glog.errormsg("");
			}
			C.W(si, si) = s / (C.err[si] * C.err[si]);
		}

		if (NLCcablen.method == "InputBunchMean") {
			double mn = C.data.mean();
			for (size_t si = 0; si < nSoundings; si++) {
				C.data[si] = mn;
			}
		}

		//std::cout << C.data.transpose() << std::endl;
		//std::cout << C.err.transpose() << std::endl;
	}

	Vector CableLengths(const Vector& m) {
		Vector cablelength(nSoundings);
		const std::vector<cTDEmGeometry> gv = get_geometry(m);
		for (size_t si = 0; si < nSoundings; si++) {
			const cTDEmGeometry& g = gv[si];
			const double dr = g.txrx_dr();
			cablelength[si] = dr;
		}
		return cablelength;
	}

	Vector CableLengthConstraint_forward(const Vector& m) {
		Vector predicted = CableLengths(m);
		if (NLCcablen.method == "BunchSimilarity") {
			double clmean = predicted.mean();
			for (size_t si = 0; si < nSoundings; si++) {
				predicted[si] = predicted[si] - clmean;
			}
		}
		return predicted;
	}

	void CableLengthConstraint_jacobian(const Vector& m) {
		cNonLinearConstraint& C = NLCcablen;
		if (C.alpha == 0.0) return;
		if (nGeomParamPerSounding <= 0) return;

		const std::vector<cTDEmGeometry> gv = get_geometry(m);
		double s = C.alpha / (double)(nSoundings);
		for (size_t si = 0; si < nSoundings; si++) {
			const cTDEmGeometry& g = gv[si];
			const double dr = g.txrx_dr();
			const int pix = gindex(si, "txrx_dx");
			const int piy = gindex(si, "txrx_dy");
			const int piz = gindex(si, "txrx_dz");

			double f = 1.0;
			if (NLCcablen.method == "BunchSimilarity") {
				f = (double)(nSoundings - 1) / (double)nSoundings;
			}

			if (pix >= 0) {
				C.J(si, pix) = f * g.txrx_dx / dr;
			}
			if (piy >= 0) {
				C.J(si, piy) = f * g.txrx_dy / dr;
			}
			if (piz >= 0) {
				C.J(si, piz) = f * g.txrx_dz / dr;
			}
		}
	}

	void initialise_Wr() {
		initialise_Wc();
		initialise_Wt();
		initialise_Wg();
		initialise_Ws();

		Wr = Matrix::Zero(nParam, nParam);
		if (LCrefc.alpha > 0.0) Wr += LCrefc.W;
		if (LCreft.alpha > 0.0) Wr += LCreft.W;
		if (LCrefg.alpha > 0.0) Wr += LCrefg.W;
		if (LCrefs.alpha > 0.0) Wr += LCrefs.W;
	}

	void initialise_Wm() {
		initialise_Wr();		
		initialise_VC();
		initialise_LC();				
		initialise_LG();
		initialise_QC(LCvcsim);
		initialise_CableLengthConstraint();
		Wm = Wr + LCvcsmth.W + LCvcsim.W + LClatc.W + LClatg.W;
	}

	void dump_W_matrices() {
		if (OO.Dump) {						
			const std::string dp = dumppath();			
			makedirectorydeep(dp);
			writetofile(Wd, dp + "Wd.dat");			
			writetofile(Wr, dp + "Wr.dat");
			
			LCrefc.write_W_matrix(dp);			
			LCreft.write_W_matrix(dp);
			LCrefg.write_W_matrix(dp);
			LCrefs.write_W_matrix(dp);
			LCvcsmth.write_W_matrix(dp);
			LCvcsim.write_W_matrix(dp);
			LClatc.write_W_matrix(dp);
			LClatg.write_W_matrix(dp);
			NLCcablen.write_W_matrix(dp);
			writetofile(Wm, dp + "Wm.dat");			
		}
	}

	const int& dindex(const size_t& sampleindex, const size_t& systemindex, const size_t& componentindex, const size_t& windowindex) {
		const size_t& si = sampleindex;
		const size_t& sysi = systemindex;		
		const size_t& ci = componentindex;
		const size_t& wi = windowindex;	
		return _dindex_[si][sysi][ci][wi];
	};

	void initialise_systems()
	{
		set_fftw_lock();
		std::vector<cBlock> B = Control.findblocks("EMSystem");
		nSystems = B.size();
		SV.resize(nSystems);		
		for (size_t sysi = 0; sysi < nSystems; sysi++) {
			SV[sysi].initialise(B[sysi], nSoundings);				
			SV[sysi].set_units(IM.get());
		}
		unset_fftw_lock();		
	}
	
	void setup_data()
	{		
		nAllData = 0;
		_dindex_.resize(nSoundings);
		for (size_t si = 0; si < nSoundings; si++) {
			_dindex_[si].resize(nSystems);
			for (size_t sysi = 0; sysi < nSystems; sysi++) {
				_dindex_[si][sysi].resize(4);//4 because of xzinversion				
				for (size_t ci = 0; ci < 4; ci++) {
					_dindex_[si][sysi][ci].resize(SV[sysi].nwindows);
					for (size_t wi = 0; wi < SV[sysi].nwindows; wi++) {
						_dindex_[si][sysi][ci][wi] = -1;
					}
				}
			}
		}		

		int di = 0;
		for (size_t si = 0; si < nSoundings; si++) {
			for (size_t sysi = 0; sysi < nSystems; sysi++) {
				cTDEmSystemInfo& S = SV[sysi];				
				if (S.invertXPlusZ) {
					nAllData += S.nwindows;
					for (size_t wi = 0; wi < S.nwindows; wi++) {
						_dindex_[si][sysi][XZAMP][wi] = di;
						di++;
					}

					if (S.CompInfo[YCOMP].Use) {
						nAllData += S.nwindows;
						for (size_t wi = 0; wi < S.nwindows; wi++) {
							_dindex_[si][sysi][YCOMP][wi] = di;
							di++;
						}
					}
				}
				else {
					for (size_t ci = 0; ci < 3; ci++) {
						cTDEmComponentInfo& c = S.CompInfo[ci];
						if (c.Use) {
							nAllData += S.nwindows;
							for (size_t wi = 0; wi < S.nwindows; wi++) {
								_dindex_[si][sysi][ci][wi] = di;
								di++;
							}
						}
					}
				}
			}
		}
	}

	bool initialise_bunch_data(){
		std::vector<double> obs(nAllData);
		std::vector<double> err(nAllData);
		std::vector<double> pred(nAllData);

		for (size_t si = 0; si < nSoundings; si++) {
			for (size_t sysi = 0; sysi < nSystems; sysi++) {
				cTDEmSystemInfo& S = SV[sysi];
				cTDEmSystem& T = S.T;
				if (S.reconstructPrimary) {
					T.setgeometry(G[si].tfr);
					T.LEM.calculation_type = cLEM::CalculationType::FORWARDMODEL;
					T.LEM.derivative_layer = INT_MAX;
					T.setprimaryfields();
					
					if (S.CompInfo[XCOMP].Use) S.CompInfo[XCOMP].data[si].P = T.PrimaryX;
					if (S.CompInfo[YCOMP].Use) S.CompInfo[YCOMP].data[si].P = T.PrimaryY;
					if (S.CompInfo[ZCOMP].Use) S.CompInfo[ZCOMP].data[si].P = T.PrimaryZ;
				}

				if (S.invertXPlusZ) {
					for (size_t wi = 0; wi < S.nwindows; wi++) {
						//XZ Amplitude						
						int di = dindex(si, sysi, XZAMP, wi);
						double X = S.CompInfo[XCOMP].data[si].S[wi];
						double Z = S.CompInfo[ZCOMP].data[si].S[wi];						
						if (S.invertPrimaryPlusSecondary) {
							X += S.CompInfo[XCOMP].data[si].P;
							Z += S.CompInfo[ZCOMP].data[si].P;														
						}
						obs[di] = std::hypot(X, Z);

						const double& Xerr = S.CompInfo[XCOMP].data[si].E[wi];
						const double& Zerr = S.CompInfo[ZCOMP].data[si].E[wi];
						if (obs[di] == 0.0) {
							err[di] = std::hypot(Xerr, Zerr);							
						}
						else {
							err[di] = std::hypot(X * Xerr, Z * Zerr) / obs[di];
						}

						//Y Comp
						if (S.CompInfo[YCOMP].Use) {
							int di = dindex(si, sysi, YCOMP, wi);							
							obs[di] = S.CompInfo[YCOMP].data[si].S[wi];
							if (S.invertPrimaryPlusSecondary) {
								obs[di] += S.CompInfo[YCOMP].data[si].P;
							}
							err[di] = S.CompInfo[YCOMP].data[si].E[wi];
						}
					}
				}
				else {
					for (size_t ci = 0; ci < 3; ci++) {
						if (S.CompInfo[ci].Use == false) continue;
						for (size_t wi = 0; wi < S.nwindows; wi++) {							
							int di = dindex(si, sysi, ci, wi);
							obs[di] = S.CompInfo[ci].data[si].S[wi];
							if (S.invertPrimaryPlusSecondary) {
								obs[di] += S.CompInfo[ci].data[si].P;
							}
							err[di] = S.CompInfo[ci].data[si].E[wi];
						}
					}
				}
			}
		}

		if (ErrorAddition > 0.0) {
			err = err + ErrorAddition;
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
		//std::cerr << std::endl << Obs << std::endl;

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
		
	void initialise_bunch_parameters() {

		for (size_t si = 0; si < nSoundings; si++) {
			const cEarthStruct& e = E[si];
			const cGeomStruct& g = G[si];
			if (solve_conductivity()) {
				for (size_t li = 0; li < nLayers; li++) {
					RefParam[cindex(si,li)] = log10(e.ref.conductivity[li]);
					RefParamStd[cindex(si,li)] = e.std.conductivity[li];
				}
			}

			if (solve_thickness()) {
				for (size_t li = 0; li < nLayers - 1; li++) {
					RefParam[tindex(si,li)] = log10(e.ref.thickness[li]);
					RefParamStd[tindex(si,li)] = e.std.thickness[li];
				}
			}

			for (int gi = 0; gi < cTDEmGeometry::size(); gi++) {
				std::string gname = cTDEmGeometry::element_name(gi);
				const int pi = gindex(si, gname);								
				if (pi >= 0) {						
					RefParam[pi] = g.ref[gname];
					RefParamStd[pi] = g.std[gname];
				}
			}

			//Scaling params				
			for (size_t sysi = 0; sysi < SV.size(); sysi++) {
				for (size_t ci = 0; ci < 3; ci++) {
					const int pi = sfindex(sysi,ci);
					if (pi >= 0) {		
						double ref = SV[sysi].CompInfo[ci].SF.ref;
						double std = SV[sysi].CompInfo[ci].SF.std;
						RefParam[pi]    = SV[sysi].CompInfo[ci].SF.ref;
						RefParamStd[pi] = SV[sysi].CompInfo[ci].SF.std;
					}					
				}
			}


		}
	}
		
	Vector parameter_change(const double& lambda, const Vector& m_old, const Vector& pred)
	{
		Vector m_new = solve_linear_system(lambda, m_old, pred);
		Vector dm = m_new - m_old;
				
		if (fdC.bound()) {			
			for (size_t si = 0; si < nSoundings; si++) {
				const cEarthStruct& e = E[si];
				for (size_t li = 0; li < nLayers; li++) {
					const int pindex = cindex(si,li);
					const double lmin = std::log10(e.min.conductivity[li]);
					const double lmax = std::log10(e.max.conductivity[li]);
					if (m_new[pindex] < lmin) {
						if (Verbose) {
							//std::cerr << rec_it_str() << std::endl;
							//std::cerr << "Lower conductivity bound reached" << std::endl;
							//std::cerr << "\t li=" << li << "\tdm=" << dm[pindex] << "\tm=" << m_old[pindex] << "\tm+dm=" << m_new[pindex] << std::endl;
							//std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(m_new[pindex]) << std::endl;
						}
						dm[pindex] = lmin - m_old[pindex];
						//if (Verbose) std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(dm[pindex] + m_old[pindex]) << std::endl;
					}
					else if (m_new[pindex] > lmax) {
						if (Verbose) {
							//std::cerr << rec_it_str() << std::endl;
							//std::cerr << "Upper conductivity bound reached" << std::endl;
							//std::cerr << "\t li=" << li << "\tdm=" << dm[pindex] << "\tm=" << m_old[pindex] << "\tm+dm=" << m_new[pindex] << std::endl;
							//std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(m_new[pindex]) << std::endl;
						}
						dm[pindex] = lmax - m_old[pindex];
						//if (Verbose) std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(dm[pindex] + m_old[pindex]) << std::endl;
					}
				}
			}
		}
		
		if (fdT.bound()) {
			for (size_t si = 0; si < nSoundings; si++) {
				const cEarthStruct& e = E[si];
				for (size_t li = 0; li < nLayers - 1; li++) {
					const int pindex = tindex(si,li);
					const double lmin = std::log10(e.min.thickness[li]);
					const double lmax = std::log10(e.max.thickness[li]);
					if (m_new[pindex] < lmin) {
						if (Verbose) {
							//std::cerr << rec_it_str() << std::endl;
							//std::cerr << "Lower thickness bound reached" << std::endl;
							//std::cerr << "\t li=" << li << "\tdm=" << dm[pindex] << "\tm=" << m_old[pindex] << "\tm+dm=" << m_new[pindex] << std::endl;
							//std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(m_new[pindex]) << std::endl;
						}
						dm[pindex] = lmin - m_old[pindex];
						//if (Verbose) std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(dm[pindex] + m_old[pindex]) << std::endl;
					}
					else if (m_new[pindex] > lmax) {
						if (Verbose) {
							//std::cerr << rec_it_str() << std::endl;
							//std::cerr << "Upper thickness bound reached" << std::endl;
							//std::cerr << "\t li=" << li << "\tdm=" << dm[pindex] << "\tm=" << m_old[pindex] << "\tm+dm=" << m_new[pindex] << std::endl;
							//std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(m_new[pindex]) << std::endl;
						}
						dm[pindex] = lmax - m_old[pindex];
						//if (Verbose) std::cerr << "\t li=" << li << "\tdm=" << pow10(dm[pindex]) << "\tm=" << pow10(m_old[pindex]) << "\tm+dm=" << pow10(dm[pindex] + m_old[pindex]) << std::endl;
					}
				}
			}
		}

		for (size_t si = 0; si < nSoundings; si++) {
			cGeomStruct& g = G[si];
			for (size_t i = 0; i < cTDEmGeometry::size(); i++) {
				const std::string ename = cTDEmGeometry::element_name(i);
				const cInvertibleFieldDefinition& e = fdG.cref(ename);
				if (e.bound()) {
					const int pi = gindex(si,ename);
					const double emin = g.min[ename];
					const double emax = g.max[ename];
					if (m_new[pi] < emin) {
						if (Verbose) {
							//std::cerr << rec_it_str() << std::endl;
							//std::cerr << "Lower " << ename << " bound reached" << std::endl;
							//std::cerr << "\tdm=" << dm[pi] << "\tm=" << m_old[pi] << "\tm+dm=" << m_new[pi] << std::endl;
						}
						dm[pi] = emin - m_old[pi];
					}
					else if (m_new[pi] > emax) {
						if (Verbose) {
							//std::cerr << rec_it_str() << std::endl;
							//std::cerr << "Upper " << ename << " bound reached" << std::endl;
							//std::cerr << "\tdm=" << dm[pi] << "\tm=" << m_old[pi] << "\tm+dm=" << m_new[pi] << std::endl;
						}
						dm[pi] = emax - m_old[pi];
					}
				}
			}
		}
		
		return dm;
	}

	std::vector<cEarth1D> get_earth(const Vector& parameters)
	{
		std::vector<cEarth1D> ev(nSoundings);;
		for (size_t si = 0; si < nSoundings; si++) {
			ev[si] = E[si].ref;
			if (solve_conductivity()) {				
				for (size_t li = 0; li < nLayers; li++) {
					ev[si].conductivity[li] = pow10(parameters[cindex(si, li)]);
				}
			}			

			if (solve_thickness()) {				
				for (size_t li = 0; li < nLayers - 1; li++) {
					ev[si].thickness[li] = pow10(parameters[tindex(si, li)]);
				}
			}			
		}
		return ev;
	}

	std::vector<cTDEmGeometry> get_geometry(const Vector& parameters)
	{
		std::vector<cTDEmGeometry> gv(nSoundings);
		for (size_t si = 0; si < nSoundings; si++) {
			gv[si] = G[si].input;
			for (int gi = 0; gi < cTDEmGeometry::size(); gi++) {
				const std::string& gname = cTDEmGeometry::element_name(gi);				
				const int pi = gindex(si, gname);
				if (pi>=0) {
					gv[si][gname] = parameters[pi];
				}
			}
		}
		return gv;
	}

	
	std::vector<double> get_scalefactors(const size_t sysi, const Vector& parameters)
	{
		std::vector<double> sf(3);		
		const cTDEmSystemInfo& S = SV[sysi];								
		//bookmark
		for (int ci = 0; ci < 3; ci++) {
			sf[ci] = 1.0;
			const int pi = sfindex(sysi, ci);
			if (pi >= 0) {
				sf[ci] = parameters[pi];
			}
		}		
		return sf;
	}
	
	
	void set_predicted(const Vector& parameters)
	{				
		std::vector<cEarth1D> ev = get_earth(parameters);
		std::vector<cTDEmGeometry> gv = get_geometry(parameters);
		for (size_t sysi = 0; sysi < nSystems; sysi++) {
			cTDEmSystemInfo& S = SV[sysi];
			S.predicted.resize(nSoundings);

			cTDEmSystem& T = S.T;
			const size_t nw = T.NumberOfWindows;
			for (size_t si = 0; si < nSoundings; si++) {
				const cEarth1D& e = ev[si];
				const cTDEmGeometry& g = gv[si];
				T.setconductivitythickness(e.conductivity, e.thickness);
				T.setgeometry(g);

				//Forwardmodel
				T.LEM.calculation_type = cLEM::CalculationType::FORWARDMODEL;
				T.LEM.derivative_layer = INT_MAX;
				T.setupcomputations();
				T.setprimaryfields();
				T.setsecondaryfields();

				cTDEmData& d = S.predicted[si];
				d.xcomponent().Primary = T.PrimaryX;
				d.ycomponent().Primary = T.PrimaryY;
				d.zcomponent().Primary = T.PrimaryZ;
				d.xcomponent().Secondary = T.X;
				d.ycomponent().Secondary = T.Y;
				d.zcomponent().Secondary = T.Z;				
			}
		}				
	}
	
	void forwardmodel(const Vector& parameters, Vector& predicted) {
		Matrix dummy;	
		nForwards++;
		forwardmodel_impl(parameters, predicted, dummy, false);		
	}

	void forwardmodel_and_jacobian(const Vector& parameters, Vector& predicted, Matrix& jacobian) {
		nForwards++;
		nJacobians++;
		forwardmodel_impl(parameters, predicted, jacobian, true);
	}

	void forwardmodel_impl(const Vector& parameters, Vector& predicted, Matrix& jacobian, bool computederivatives)
	{
		Vector pred_all(nAllData);
		Matrix J_all;
		if (computederivatives) {
			J_all.resize(nAllData, nParam);
			J_all.setZero();
		}

		std::vector<cEarth1D> ev = get_earth(parameters);
		std::vector<cTDEmGeometry> gv = get_geometry(parameters);				
		for (size_t sysi = 0; sysi < nSystems; sysi++) {			
			cTDEmSystemInfo& S = SV[sysi];
			cTDEmSystem& T = S.T;

			std::vector<double> scalefactors = get_scalefactors(sysi, parameters);

			const size_t nw = T.NumberOfWindows;
			for (size_t si = 0; si < nSoundings; si++) {
				const cEarth1D& e = ev[si];
				const cTDEmGeometry& g = gv[si];
				T.setconductivitythickness(e.conductivity, e.thickness);
				T.setgeometry(g);

				//Forwardmodel
				T.LEM.calculation_type = cLEM::CalculationType::FORWARDMODEL;
				T.LEM.derivative_layer = INT_MAX;
				T.setupcomputations();
				T.setprimaryfields();
				T.setsecondaryfields();

				std::vector<double> xfm = T.X * scalefactors[XCOMP];
				std::vector<double> yfm = T.Y * scalefactors[YCOMP];
				std::vector<double> zfm = T.Z * scalefactors[ZCOMP];
				std::vector<double> xzfm;
				if (S.invertPrimaryPlusSecondary) {
					xfm += T.PrimaryX * scalefactors[XCOMP];
					yfm += T.PrimaryY * scalefactors[YCOMP];
					zfm += T.PrimaryZ * scalefactors[ZCOMP];
				}

				if (S.invertXPlusZ) {
					xzfm.resize(T.NumberOfWindows);
					for (size_t wi = 0; wi < T.NumberOfWindows; wi++) {
						xzfm[wi] = std::hypot(xfm[wi], zfm[wi]);
					}
				}

				if (S.invertXPlusZ) {
					for (size_t wi = 0; wi < nw; wi++) {
						const int& di = dindex(si, sysi, XZAMP, wi);
						pred_all[di] = xzfm[wi];
						if (S.CompInfo[1].Use){
							pred_all[dindex(si, sysi, YCOMP, wi)] = yfm[wi];
						}
					}
				}
				else {
					for (size_t wi = 0; wi < nw; wi++) {
						if (S.CompInfo[XCOMP].Use) pred_all[dindex(si, sysi, XCOMP, wi)] = xfm[wi];
						if (S.CompInfo[YCOMP].Use) pred_all[dindex(si, sysi, YCOMP, wi)] = yfm[wi];
						if (S.CompInfo[ZCOMP].Use) pred_all[dindex(si, sysi, ZCOMP, wi)] = zfm[wi];
					}
				}

				if (computederivatives) {					
					std::vector<double> xdrv(nw);
					std::vector<double> ydrv(nw);
					std::vector<double> zdrv(nw);
					
					//bookmark
					for (size_t ci = 0; ci < 3; ci++) {
						if (S.CompInfo[ci].Use) {
							const int pindex = sfindex(sysi, ci);
							if (pindex >= 0) {
								//Here filling with the forward itself as no new computations
								fillDerivativeVectors(S, xdrv, ydrv, zdrv);
								if (ci != XCOMP) {
									xdrv *= 0.0;
								}
								if (ci != YCOMP) {
									ydrv *= 0.0;
								}
								if (ci != ZCOMP) {
									zdrv *= 0.0;
								}
								fillMatrixColumn(J_all, si, sysi, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
							}
						}
					}		

					if (solve_conductivity()) {
						for (size_t li = 0; li < nLayers; li++) {
							const int pindex = cindex(si,li);
							T.LEM.calculation_type = cLEM::CalculationType::CONDUCTIVITYDERIVATIVE;
							T.LEM.derivative_layer = li;
							T.setprimaryfields();
							T.setsecondaryfields();

							fillDerivativeVectors(S, xdrv, ydrv, zdrv);
							//multiply by natural log(10) as parameters are in logbase10 units
							const double f = log(10.0) * e.conductivity[li];
							xdrv *= f; ydrv *= f; zdrv *= f;
							fillMatrixColumn(J_all, si, sysi, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
						}
					}

					if (solve_thickness()) {
						for (size_t li = 0; li < nLayers - 1; li++) {
							const int pindex = tindex(si,li);
							T.LEM.calculation_type = cLEM::CalculationType::THICKNESSDERIVATIVE;
							T.LEM.derivative_layer = li;
							T.setprimaryfields();
							T.setsecondaryfields();
							fillDerivativeVectors(S, xdrv, ydrv, zdrv);
							//multiply by natural log(10) as parameters are in logbase10 units
							double sf = log(10.0) * e.thickness[li];
							xdrv *= sf; ydrv *= sf; zdrv *= sf;
							fillMatrixColumn(J_all, si, sysi, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
						}
					}

					if (FreeGeometry) {
						if (solve_geometry_element("tx_height")) {							
							const size_t pindex = gindex(si, "tx_height");
							T.LEM.calculation_type = cLEM::CalculationType::HDERIVATIVE;
							T.LEM.derivative_layer = INT_MAX;
							T.setprimaryfields();
							T.setsecondaryfields();
							fillDerivativeVectors(S, xdrv, ydrv, zdrv);
							fillMatrixColumn(J_all, si, sysi, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
						}

						if (solve_geometry_element("txrx_dx")) {							
							const size_t pindex = gindex(si,"txrx_dx");
							T.LEM.calculation_type = cLEM::CalculationType::XDERIVATIVE;
							T.LEM.derivative_layer = INT_MAX;
							T.setprimaryfields();
							T.setsecondaryfields();
							fillDerivativeVectors(S, xdrv, ydrv, zdrv);
							fillMatrixColumn(J_all, si, sysi, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
						}

						if (solve_geometry_element("txrx_dy")) {							
							const size_t pindex = gindex(si, "txrx_dy");
							T.LEM.calculation_type = cLEM::CalculationType::YDERIVATIVE;
							T.LEM.derivative_layer = INT_MAX;
							T.setprimaryfields();
							T.setsecondaryfields();
							fillDerivativeVectors(S, xdrv, ydrv, zdrv);
							fillMatrixColumn(J_all, si, sysi, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
						}

						if (solve_geometry_element("txrx_dz")) {
							const size_t pindex = gindex(si, "txrx_dz");
							T.LEM.calculation_type = cLEM::CalculationType::ZDERIVATIVE;
							T.LEM.derivative_layer = INT_MAX;
							T.setprimaryfields();
							T.setsecondaryfields();
							fillDerivativeVectors(S, xdrv, ydrv, zdrv);
							fillMatrixColumn(J_all, si, sysi, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
						}

						if (solve_geometry_element("rx_pitch")) {
							const size_t pindex = gindex(si, "rx_pitch");
							T.drx_pitch(xfm, zfm, g.rx_pitch, xdrv, zdrv);
							ydrv *= 0.0;
							fillMatrixColumn(J_all, si, sysi, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
						}
						
						if (solve_geometry_element("rx_roll")) {
							const size_t pindex = gindex(si, "rx_roll");
							T.drx_roll(yfm, zfm, g.rx_roll, ydrv, zdrv);
							xdrv *= 0.0;
							fillMatrixColumn(J_all, si, sysi, pindex, xfm, yfm, zfm, xzfm, xdrv, ydrv, zdrv);
						}
					}
				}
			}
		}
		predicted = cull(pred_all);
		if(computederivatives) jacobian = cull(J_all);

		if (Verbose && computederivatives) {
			//std::cerr << "\n-----------------\n";
			//std::cerr << "J_all: It " << CIS.iteration + 1 << std::endl;			
			//std::cerr << J_all;			
			//std::cerr << "\n-----------------\n";
		}

		if (OO.Dump && computederivatives) {
			const std::string dp = dumppath();
			writetofile(J_all, dp + "J" + ".dat");
			std::ofstream of(dp + "J1" + ".dat");
			of << J_all;
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

	void fillMatrixColumn(Matrix& M, const size_t& si, const size_t& sysi, const size_t& pindex, const std::vector<double>& xfm, const std::vector<double>& yfm, const std::vector<double>& zfm, const std::vector<double>& xzfm, const std::vector<double>& xdrv, const std::vector<double>& ydrv, const std::vector<double>& zdrv)
	{
		const cTDEmSystemInfo& S = SV[sysi];
		const size_t& nw = S.T.NumberOfWindows;
		if (S.invertXPlusZ) {
			for (size_t wi = 0; wi < nw; wi++) {								
				M(dindex(si, sysi, XZAMP, wi), pindex) = (xfm[wi] * xdrv[wi] + zfm[wi] * zdrv[wi]) / xzfm[wi];
				if (S.CompInfo[1].Use) {
					M(dindex(si, sysi, YCOMP, wi), pindex) = ydrv[wi];
				}
			}
		}
		else {
			for (size_t wi = 0; wi < nw; wi++) {				
				if (S.CompInfo[XCOMP].Use)M(dindex(si, sysi, XCOMP, wi),pindex) = xdrv[wi];
				if (S.CompInfo[YCOMP].Use)M(dindex(si, sysi, YCOMP, wi),pindex) = ydrv[wi];
				if (S.CompInfo[ZCOMP].Use)M(dindex(si, sysi, ZCOMP, wi),pindex) = zdrv[wi];
			}
		}
	}

	void save_iteration_file(const cIterationState& S) {
		std::ofstream ofs(dumppath() + "iteration.dat");
		ofs << S.info_string();
	};
	
	bool read_bunch(const size_t& record) {
		_GSTITEM_

		bool bunchstatus=false;
		int fi = AncFld[0].keyindex("line");
		if (fi > 0) {
			//only if linenumber is specified
			cFieldDefinition& fdline = AncFld[0][fi].second.fd;
			bunchstatus = IM->get_bunch(Bunch, fdline, (int)record, (int)nSoundings, (int)nBunchSubsample);
		}
		else {
			cFieldDefinition fdnone;
			bunchstatus = IM->get_bunch(Bunch, fdnone, (int)record, (int)nSoundings, (int)nBunchSubsample);
		}

		if (bunchstatus == false) {
			return bunchstatus;
		}

		for (size_t si = 0; si < Bunch.size(); si++) {
			const size_t& record = Bunch.record(si);
			bool loadstatus = IM->load_record(record);
			if (loadstatus == false) {
				OutputMessage += ", Skipping - could not load record";
				return false;
			}
			bool valid = IM->is_record_valid();
			if (valid == false) {
				OutputMessage += ", Skipping - record is not valid";
				return false;
			}
			bool readstatus = read_record(si);
			if (valid == false) {
				OutputMessage += ", Skipping - could not read record";
				return false;
			}
		}

		//bookmark
		for (size_t sysi = 0; sysi < SV.size(); sysi++) {
			cTDEmSystemInfo& S = SV[sysi];
			for (size_t ci = 0; ci < 3; ci++) {
				cTDEmComponentInfo& C = S.CompInfo[ci];
				if (C.fdSF.solve) {
					IM->read(C.fdSF.ref, C.SF.ref);
					IM->read(C.fdSF.std, C.SF.std);
					IM->read(C.fdSF.min, C.SF.min);
					IM->read(C.fdSF.max, C.SF.max);
				}
			}
		}

		return true;			
	}

	bool read_record(const size_t& bunchsoundingindex){
		const size_t& si = bunchsoundingindex;
		bool readstatus = true;
		cEarthStruct& e = E[si];
		cGeomStruct& g = G[si];

		if (IM->parse_record() == false) return false;

		bool status;
		Id[si].uniqueid = (int)IM->record();

		status = read_ancillary_fields(si);
		status = read_geometry(si, fdG);
		status = IM->read(fdC.input, e.ref.conductivity, nLayers); if (status == false) readstatus = false;
		if (solve_conductivity()) {
			status = IM->read(fdC.ref, e.ref.conductivity, nLayers); if (status == false) readstatus = false;
			status = IM->read(fdC.std, e.std.conductivity, nLayers); if (status == false) readstatus = false;
			status = IM->read(fdC.min, e.min.conductivity, nLayers); if (status == false) readstatus = false;
			status = IM->read(fdC.max, e.max.conductivity, nLayers); if (status == false) readstatus = false;
		}

		status = IM->read(fdT.input, e.ref.thickness, nLayers - 1); if (status == false) readstatus = false;
		if (solve_thickness()) {
			status = IM->read(fdT.ref, e.ref.thickness, nLayers - 1); if (status == false) readstatus = false;
			status = IM->read(fdT.std, e.std.thickness, nLayers - 1); if (status == false) readstatus = false;
			status = IM->read(fdT.min, e.min.thickness, nLayers - 1); if (status == false) readstatus = false;
			status = IM->read(fdT.max, e.max.thickness, nLayers - 1); if (status == false) readstatus = false;
		}
		e.sanity_check();

		status = IM->read(fdT.max, e.max.thickness, nLayers - 1); if (status == false) readstatus = false;

		for (size_t sysi = 0; sysi < nSystems; sysi++) {
			read_system_data(sysi, si);			
		}
		return readstatus;
	}
	
	bool read_ancillary_fields(const size_t& bunchindex) {
		const size_t& si = bunchindex;
		SampleId& id = Id[si];

		for (size_t fi = 0; fi < AncFld[si].size(); fi++) {
			IM->readfdvnt(AncFld[si][fi].second);
		}
				
		set_ancillary_id(si, "Survey", id.survey);
		set_ancillary_id(si, "Date", id.date);
		set_ancillary_id(si, "Flight", id.flight);
		set_ancillary_id(si, "Line", id.line);
		set_ancillary_id(si, "Fiducial", id.fiducial);		
		set_ancillary_id(si, "X", id.x);
		set_ancillary_id(si, "Y", id.y);
		set_ancillary_id(si, "GroundElevation", id.elevation);
		return true;
	}

	template<typename T>
	bool set_ancillary_id(const size_t& si, const std::string key, T& value) {
		int ki = AncFld[si].keyindex(key);
		if (ki >= 0) {
			//value = std::get<T>(AncFld[si][ki].second.vnt);						
			const cVrnt& v = AncFld[si][ki].second.vnt;						
			if (v.index() == 0) {
				value = (T)std::get<double>(v);
			}
			else if (v.index() == 1) {
				value = (T)std::get<int>(v);
			}
			else if (v.index() == 2) {
				value = (T)std::get<float>(v);
			}
			else {
				glog.errormsg(_SRC_,"Bad variant type\n");
			}
			
			return true;
		}
		return false;
	}

	bool read_geometry(const size_t& bunchindex, cIFDMap& map)
	{
		bool status = true;
		const size_t si = bunchindex;
		cGeomStruct& g = G[si];
		for (size_t gi = 0; gi < cTDEmGeometry::size(); gi++) {
			std::string ename = cTDEmGeometry::element_name(gi);
			const cInvertibleFieldDefinition ge = map.cref(ename);
			bool inpstatus = IM->read(ge.input, g.input[gi]);
			bool refstatus = IM->read(ge.ref, g.ref[gi]);

			if (refstatus == false && inpstatus == true) {
				g.ref[gi] = g.input[gi];
				refstatus = true;
			}
			else if (inpstatus == false && refstatus == true) {
				g.input[gi] = g.ref[gi];
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

			bool tfrstatus = IM->read(ge.tfr, g.tfr[gi]);
			if (tfrstatus == false) {
				g.tfr[gi] = g.input[gi];
			}

			if (ge.solve) {
				bool stdstatus = IM->read(ge.std, g.std[gi]);
				if (stdstatus == false) {
					std::ostringstream msg;
					msg << "Error: no 'Std' defined for " << ename << std::endl;
					glog.errormsg(msg.str());
				}

				bool minstatus = IM->read(ge.min, g.min[gi]);
				bool maxstatus = IM->read(ge.max, g.max[gi]);
			}
		}
		return status;
	}
		
	void read_system_data(size_t& sysindex, const size_t& soundingindex)
	{
		cTDEmSystemInfo& S = SV[sysindex];
		S.CompInfo[XCOMP].readdata(IM, soundingindex);
		S.CompInfo[YCOMP].readdata(IM, soundingindex);
		S.CompInfo[ZCOMP].readdata(IM, soundingindex);
	}

	void dump_first_iteration() {
		
		const std::string dp = dumppath();
		makedirectorydeep(dumppath());

		const size_t si = Bunch.master_index();
		cGeomStruct& g = G[si];
		cEarthStruct& e = E[si];
		SampleId& id = Id[si];
		
		write(Obs, dp + "observed.dat");
		write(Err, dp + "observed_std.dat");

		g.ref.write(dp + "geometry_start.dat");
		e.ref.write(dp + "earth_start.dat");

		g.ref.write(dp + "geometry_ref.dat");
		e.ref.write(dp + "earth_ref.dat");

		g.std.write(dp + "geometry_std.dat");
		e.std.write(dp + "earth_std.dat");

		std::ofstream ofs(dp + "Id.dat");
		char sep = '\n';

		ofs << id.uniqueid << sep;
		ofs << id.survey << sep;
		ofs << id.date << sep;
		ofs << id.flight << sep;
		ofs << id.line << sep;
		ofs << id.fiducial << sep;
		ofs << id.x << sep;
		ofs << id.y << sep;
		ofs << id.elevation << sep;
	}

	void dump_iteration(const cIterationState& state) {
		const std::string dp = dumppath();
		makedirectorydeep(dp);
		writetofile(Obs, dp + "d.dat");
		writetofile(Err, dp + "e.dat");
		writetofile(state.param, dp + "m.dat");
		writetofile(state.pred, dp + "g.dat");
		std::vector<cEarth1D> e = get_earth(state.param);
		std::vector <cTDEmGeometry> g = get_geometry(state.param);
		e[Bunch.master_index()].write(dumppath() + "earth_inv.dat");
		g[Bunch.master_index()].write(dumppath() + "geometry_inv.dat");		
		dump_earth_all(e, dumppath() + "earth_all.dat");
		dump_geometry_all(g, dumppath() + "geometry_all.dat");
		save_iteration_file(state);
	}

	void dump_earth_all(const std::vector <cEarth1D> e, const std::string& path) {
		std::ofstream of(path);
		for (size_t si = 0; si < e.size(); si++) {
			const std::vector<double>& c = e[si].conductivity;
			const std::vector<double>t   = e[si].dummy_thickness();						
			for (size_t li = 0; li < e[si].nlayers(); li++) {
				of << t[li] << " " << c[li] << std::endl;
			}
		}
	}

	void dump_geometry_all(const std::vector <cTDEmGeometry> g, const std::string& path){
		std::ofstream of(path);
		for (size_t si = 0; si < g.size(); si++) {
			for (size_t gi = 0; gi < g[si].size(); gi++) {
				of << g[si][gi] << std::endl;
			}
		}
	}

	bool initialise_bunch() {	
		nForwards = 0;
		nJacobians = 0;
		OutputMessage = "";		
		CIS = cIterationState();		
		bool status = initialise_bunch_data();
		if (status == false) return false;
		initialise_bunch_parameters();
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
		CIS.phim = phiModel(CIS.param);

		TerminationReason = "Has not terminated";

		if (OO.Dump) {
			dump_first_iteration();
			dump_iteration(CIS);
		}
					
		double percentchange = 100.0;
		bool   keepiterating = true;
		while (keepiterating == true) {			
			if (Verbose && nScalingParam >0) {			
				std::vector<double> scalefactors = get_scalefactors(0, CIS.param);
				std::cout << "Scaling Factors ";
				for (size_t ci = 0; ci < 3; ci++) {
					std::cout << scalefactors[ci] << " ";
				}
				std::cout << std::endl;
			}

			if (CIS.iteration >= MaxIterations) {
				keepiterating = false;
				TerminationReason = "Too many iterations";
			}
			else if (CIS.iteration>0 && CIS.phid <= MinimumPhiD) {
				keepiterating = false;
				TerminationReason = "Reached minimum";
			}			
			else if (CIS.phid < 50.0 && CIS.iteration > 10 && percentchange < MinimumImprovement) {
				keepiterating = false;
				TerminationReason = "Small % improvement";
			}
			else {			
				if (Verbose) {
					std::cerr << CIS.info_string();
				}
				if (CIS.iteration+1 >= BeginGeometrySolveIteration) FreeGeometry = true;
				else FreeGeometry = false;
								
				Vector g;
				forwardmodel_and_jacobian(CIS.param, g, J);
				//if (CIS.iteration == 0) {
				//	CIS.lambda = estimate_initial_lambda();
				//	std::cerr << "Initial lambda = " << CIS.lambda << std::endl;
				//}

				double targetphid = std::max(CIS.phid*0.7, MinimumPhiD);

				cTrial t  = lambda_search_target(CIS.lambda, targetphid);
				Vector dm = parameter_change(t.lambda, CIS.param, CIS.pred);
				Vector m = CIS.param + (t.stepfactor * dm);
				
				forwardmodel(m,g);
				double phid = phiData(g);

				percentchange = 100.0 * (CIS.phid - phid) / (CIS.phid);
				if (phid <= CIS.phid) {				
					CIS.iteration++;
					CIS.param = m;
					CIS.pred = g;
					CIS.targetphid = targetphid;										
					CIS.phid   = phid;
					CIS.lambda = t.lambda;
					CIS.phim = phiModel(CIS.param);
					if (OO.Dump) dump_iteration(CIS);
				}						
			}			
		} 
		
		std::vector<cEarth1D> ev = get_earth(CIS.param);
		std::vector<cTDEmGeometry> gv = get_geometry(CIS.param);
		for (size_t si = 0; si < nSoundings; si++) {			
			E[si].invmodel = ev[si];
			G[si].invmodel = gv[si];
		}

		set_predicted(CIS.param);
		forwardmodel_and_jacobian(CIS.param, CIS.pred, J);		
		ParameterSensitivity = compute_parameter_sensitivity();
		ParameterUncertainty = compute_parameter_uncertainty();
	}
	
	int execute() {
		_GSTITEM_				
		bool readstatus = true;
		int paralleljob = 0;			
		do{												
			int record = ((int)StartRecord-1) + paralleljob*(int)IM->subsamplerate();			
			if (record > (EndRecord - 1))break;
			if ((paralleljob % Size) == Rank) {					
				std::ostringstream s;				
				if ((readstatus = read_bunch(record))) {					
					s << bunch_id();
					if (initialise_bunch()) {
						double t1 = gettime();
						iterate();
						double t2 = gettime();
						double etime = t2 - t1;
						write_result(record);						
						s << bunch_result(etime);																		
					}
					else {
						OutputMessage += ", Skipping - could not initialise the bunch";						
					}
					s << std::endl;
					if (OutputMessage.size() > 0) {
						std::cerr << s.str();
					}
					glog.logmsg(s.str());
				}								
			}						
			paralleljob++;
		} while (readstatus == true);
		glog.close();
		return 0;
	}	

	double phiModel(const Vector& m)
	{
		double v = 0.0;
		v += LCrefc.phi(m, RefParam);
		v += LCreft.phi(m, RefParam);
		v += LCrefg.phi(m, RefParam);
		v += LCrefs.phi(m, RefParam);
		v += LCvcsmth.phi(m, RefParam);
		v += LCvcsim.phi(m, RefParam);
		v += LClatc.phi(m, RefParam);
		v += LClatg.phi(m, RefParam);

		Vector clfwd = CableLengthConstraint_forward(m);
		v += NLCcablen.phi(clfwd);
		return v;
	}

	double estimate_initial_lambda()
	{				
		Matrix JtWdJ = J.transpose() * Wd * J;
		
		Eigen::JacobiSVD<Matrix> svd0(JtWdJ);
		Vector s0 = svd0.singularValues();
		//std::cerr << "s0" << std::endl << s0 << std::endl;

		Eigen::JacobiSVD<Matrix> svd1(Wm);
		Vector s1 = svd1.singularValues();
		//std::cerr << "s1" << std::endl << s1 << std::endl;

		//std::cerr << "ratio " << s0[0]/s1[0] << std::endl;		

		double elambda = s0[0] / s1[0] * 1.0e4;
		return elambda;
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
		
		Matrix A = JtVJ + lambda * Wm;
		Vector b = JtV * (d - g + J * m);
		b += lambda * (Wr * m0);

		if (LClatg.operates_on_difference_from_reference_model()) {
			b += lambda * (LClatg.W * m0);
		}

		cNonLinearConstraint& C = NLCcablen;
		if (C.alpha > 0) {		
			CableLengthConstraint_jacobian(m);
			Vector predicted = CableLengthConstraint_forward(m);
			A += C.J.transpose() * C.W.transpose() * C.J;						
			b += C.J.transpose() * C.W.transpose() * (C.data - predicted + C.J * m);
		}
				
		const Eigen::LLT<Matrix> lltOfA(A);										
		if (lltOfA.info() == Eigen::NumericalIssue)
		{
			std::cerr << "The matrix A is possibly non semi-positive definite" << std::endl << A << std::endl;						
		}
		Vector x = lltOfA.solve(b);		
		return x;
	}
	
	void write_result(const int& pointindex)
	{
		const Vector& m = CIS.param;
		const Vector& m0 = RefParam;

		const int& pi = (int)Bunch.master_record();
		const int& si = (int)Bunch.master_index();
		OM->begin_point_output();

		//Ancillary	
		OM->writefield(pi, Id[si].uniqueid, "uniqueid", "Inversion sequence number", UNITLESS, 1, ST_UINT, DN_NONE, 'I', 12, 0);
		for (size_t fi = 0; fi < AncFld[si].size(); fi++) {
			cFdVrnt& fdv = AncFld[si][fi].second;
			cAsciiColumnField c;			
			IM->get_acsiicolumnfield(fdv.fd, c);
			OM->writevrnt(pi, fdv.vnt, c);
		}

		//Geometry Input
		bool invertedfieldsonly = false;
		for (size_t i = 0; i < G[si].input.size(); i++) {
			if (invertedfieldsonly && solve_geometry_index(i) == false)continue;
			OM->writefield(pi, G[si].input[i], "input_" + G[si].input.element_name(i), "Input " + G[si].input.description(i), G[si].input.units(i), 1, ST_FLOAT, DN_NONE, 'F', 9, 2);
		}

		//Geometry Modelled		
		const cTDEmGeometry& g = G[si].invmodel;
		invertedfieldsonly = true;
		for (size_t gi = 0; gi < g.size(); gi++) {
			if (invertedfieldsonly && solve_geometry_index(gi) == false)continue;
			OM->writefield(pi, g[gi], "inverted_" + g.element_name(gi), "Inverted " + g.description(gi), g.units(gi), 1, cOutputField::binarystoragetype::FLOAT, DN_NONE, 'F', 9, 2);
		}

		//ndata
		OM->writefield(pi,
			nData, "ndata", "Number of data in inversion", UNITLESS,
			1, ST_UINT, DN_NONE, 'I', 4, 0);

		//Scaling factors
		if (solve_scalingfactors()) {
			for (size_t sysi = 0; sysi < nSystems; sysi++) {
				cTDEmSystemInfo& S = SV[sysi];
				std::vector<double> sf = get_scalefactors(sysi, m);
				for (size_t ci = 0; ci < 3; ci++) {					
					if (S.CompInfo[ci].Use) {
						std::string comp = S.CompInfo[ci].Name;						
						std::string fname = "scalingfactor" + strprint("_EMSystem_%d_", (int)sysi + 1) + comp;
						std::string fdesc = "Scaling factor" + strprint(" EMSystem %d ", (int)sysi + 1) + comp +"-component";
											
						OM->writefield(pi,
							sf[ci], fname, fdesc, UNITLESS,
							1, ST_FLOAT, DN_NONE, 'F', 6, 3);						
					}
				}
			}
		}

		//Earth	
		const cEarth1D& e = E[si].invmodel;
		OM->writefield(pi,
			nLayers, "nlayers", "Number of layers ", UNITLESS,
			1, ST_UINT, DN_NONE, 'I', 4, 0);

		OM->writefield(pi,
			e.conductivity, "conductivity", "Layer conductivity", "S/m",
			e.conductivity.size(), cOutputField::binarystoragetype::FLOAT, DN_LAYER, 'E', 15, 6);

		if (nLayers > 1) {
			double bottomlayerthickness = 100.0;
			if (solve_thickness() == false && nLayers > 1) {
				bottomlayerthickness = e.thickness[nLayers - 2];
			}
			std::vector<double> thickness = e.thickness;
			thickness.push_back(bottomlayerthickness);

			OM->writefield(pi,
				thickness, "thickness", "Layer thickness", "m",
				thickness.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);


			if (OO.PositiveLayerTopDepths) {
				std::vector<double> dtop = e.layer_top_depth();
				OM->writefield(pi,
					dtop, "depth_top", "Depth to top of layer", "m",
					dtop.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);
			}

			if (OO.NegativeLayerTopDepths) {
				std::vector<double> ndtop = -1.0 * e.layer_top_depth();
				OM->writefield(pi,
					ndtop, "depth_top_negative", "Negative of depth to top of layer", "m",
					ndtop.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);
			}

			if (OO.PositiveLayerBottomDepths) {
				std::vector<double> dbot = e.layer_bottom_depth();
				OM->writefield(pi,
					dbot, "depth_bottom", "Depth to bottom of layer", "m",
					dbot.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);
			}

			if (OO.NegativeLayerBottomDepths) {
				std::vector<double> ndbot = -1.0 * e.layer_bottom_depth();
				OM->writefield(pi,
					ndbot, "depth_bottom_negative", "Negative of depth to bottom of layer", "m",
					ndbot.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);
			}

			if (OO.InterfaceElevations) {
				std::vector<double> etop = e.layer_top_depth();
				etop += Id[si].elevation;
				OM->writefield(pi,
					etop, "elevation_interface", "Elevation of interface", "m",
					etop.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);
			}
		}

		if (OO.ParameterSensitivity) {
			std::vector<double> ps = copy(ParameterSensitivity);
			if (solve_conductivity()) {
				std::vector<double> v(ps.begin() + cindex(si, 0), ps.begin() + cindex(si, 0) + nLayers);
				OM->writefield(pi,
					v, "conductivity_sensitivity", "Conductivity parameter sensitivity", UNITLESS,
					v.size(), ST_FLOAT, DN_LAYER, 'E', 15, 6);
			}

			if (solve_thickness()) {
				std::vector<double> v(ps.begin() + tindex(si, 0), ps.begin() + tindex(si, 0) + nLayers - 1);
				v.push_back(0.0);//halfspace layer not a parameter
				OM->writefield(pi,
					v, "thickness_sensitivity", "Thickness parameter sensitivity", UNITLESS,
					v.size(), ST_FLOAT, DN_LAYER, 'E', 15, 6);
			}

			const cTDEmGeometry& g = G[si].input;
			for (size_t gi = 0; gi < g.size(); gi++) {
				if (solve_geometry_index(gi) == true) {
					const std::string& gname = g.element_name(gi);
					std::string name = "inverted_" + gname + "_sensitivity";
					std::string desc = g.description(gi) + " parameter sensitivity";
					OM->writefield(pi,
						ps[gindex(si, gname)], name, desc, UNITLESS,
						1, ST_FLOAT, DN_NONE, 'E', 15, 6);
				}
			}
		}

		if (OO.ParameterUncertainty) {
			std::vector<double> pu = copy(ParameterUncertainty);
			if (solve_conductivity()) {
				std::vector<double> v(pu.begin() + cindex(si, 0), pu.begin() + cindex(si, 0) + nLayers);
				OM->writefield(pi,
					v, "conductivity_uncertainty", "Conductivity parameter uncertainty", "log10(S/m)",
					v.size(), ST_FLOAT, DN_LAYER, 'E', 15, 6);
			}

			if (solve_thickness()) {
				std::vector<double> v(pu.begin() + tindex(si, 0), pu.begin() + tindex(si, 0) + nLayers - 1);
				v.push_back(0.0);//halfspace layer not a parameter
				OM->writefield(pi,
					v, "thickness_uncertainty", "Thickness parameter uncertainty", "log10(m)",
					v.size(), ST_FLOAT, DN_LAYER, 'E', 15, 6);
			}

			const cTDEmGeometry& g = G[si].input;
			for (size_t gi = 0; gi < g.size(); gi++) {
				if (solve_geometry_index(gi) == false) continue;
				const std::string& gname = g.element_name(gi);
				std::string name = "inverted_" + gname + "_uncertainty";
				std::string desc = g.description(gi) + " parameter uncertainty";
				OM->writefield(pi,
					pu[gindex(si, gname)], name, desc, g.units(gi),
					1, ST_FLOAT, DN_NONE, 'E', 15, 6);
			}
		}


		//ObservedData
		if (OO.ObservedData) {
			for (size_t sysi = 0; sysi < nSystems; sysi++) {
				cTDEmSystemInfo& S = SV[sysi];		
				for (size_t ci = 0; ci < 3; ci++) {					
					if (S.CompInfo[ci].Use) writeresult_emdata(pi,
						sysi, S.CompInfo[ci].Name,
						"observed", "Observed",
						'E', 15, 6, S.CompInfo[ci].data[si].P, S.CompInfo[ci].data[si].S, S.invertPrimaryPlusSecondary,S.units);
				}
			}
		}


		//Noise Estimates
		if (OO.NoiseEstimates) {
			for (size_t sysi = 0; sysi < nSystems; sysi++) {
				cTDEmSystemInfo& S = SV[sysi];
				for (size_t ci = 0; ci < 3; ci++) {
					if (S.CompInfo[ci].Use) writeresult_emdata(pi,
						sysi, S.CompInfo[ci].Name,
						"noise", "Estimated noise",
						'E', 15, 6, 0.0, S.CompInfo[ci].data[si].E, false, S.units);
				}
			}
		}

		//PredictedData
		if (OO.PredictedData) {
			for (size_t sysi = 0; sysi < nSystems; sysi++) {
				cTDEmSystemInfo& S = SV[sysi];
				for (size_t ci = 0; ci < 3; ci++) {
					if (S.CompInfo[ci].Use) writeresult_emdata(pi,
						sysi, S.CompInfo[ci].Name, "predicted", "Predicted", 'E', 15, 6,
						S.predicted[si].component(ci).Primary,
						S.predicted[si].component(ci).Secondary,
						S.invertPrimaryPlusSecondary,S.units);
				}
			}
		}

		//Inversion parameters and norms				
		write_result(pi, LCrefc, m, m0);
		write_result(pi, LCreft, m, m0);
		write_result(pi, LCrefg, m, m0);
		write_result(pi, LCrefs, m, m0);
		write_result(pi, LCvcsmth, m, m0);
		write_result(pi, LCvcsim, m, m0);
		write_result(pi, LClatc, m, m0);
		write_result(pi, LClatg, m, m0);

		Vector clfwd = CableLengthConstraint_forward(m);
		write_result(pi, NLCcablen, clfwd);

		OM->writefield(pi, CIS.phid, "PhiD", "Normalised data misfit", UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, CIS.phim, "PhiM", "Combined model norm", UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, CIS.lambda, "Lambda", "Lambda regularization parameter", UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pi, CIS.iteration, "Iterations", "Number of iterations", UNITLESS, 1, ST_UINT, DN_NONE, 'I', 4, 0);

		//End of record book keeping
		OM->end_point_output();
		if (pointsoutput == 0) {
			OM->end_first_record();//only do this once		
		}
		pointsoutput++;
	};

	void write_result(const int& pointindex, const cLinearConstraint& C, const Vector& m, const Vector& m0) {
		if (C.alpha == 0.0) return;

		double phi = 0.0;
		if (C.alpha > 0.0) {
			phi = C.phi(m, m0);
		}

		OM->writefield(pointindex, C.alpha, C.alpha_field_name(), C.alpha_field_description(), UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pointindex, phi, C.phi_field_name(), C.phi_field_description(), UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
	}

	void write_result(const int& pointindex, const cNonLinearConstraint& C, const Vector& predicted) {
		if (C.alpha == 0.0)return;
		double phi = 0.0;
		if (C.alpha > 0.0) {
			phi = C.phi(predicted);
		}
		OM->writefield(pointindex, C.alpha, C.alpha_field_name(), C.description, UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
		OM->writefield(pointindex, phi, C.phi_field_name(), C.phi_field_description(), UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
	}

	void writeresult_emdata(const int& pointindex, const size_t& sysnum, const std::string& comp, const std::string& nameprefix, const std::string& descprefix, const char& form, const int& width, const int& decimals, const double& p, std::vector<double>& s, const bool& includeprimary, const std::string& units)
	{
		std::string DN_WINDOW = "em_window";
		std::string sysname = nameprefix + strprint("_EMSystem_%d_", (int)sysnum + 1);
		std::string sysdesc = descprefix + strprint(" EMSystem %d ", (int)sysnum + 1);		
		if (includeprimary) {
			std::string name = sysname + comp + "P";
			std::string desc = sysdesc + comp + "-component primary field";
			OM->writefield(pointindex,
				p, name, desc,units,
				1, ST_FLOAT, DN_NONE, form, width, decimals);
		}

		{
			std::string name = sysname + comp + "S";
			std::string desc = sysdesc + comp + "-component secondary field";
			OM->writefield(pointindex,
				s, name, desc, units,
				s.size(), ST_FLOAT, DN_WINDOW, form, width, decimals);
		}
	}
};

#endif
