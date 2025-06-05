/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include <stdio.h>
#include <sstream>
#include <vector>
#include <cstring>
#include <algorithm>
#include <iomanip>
#include <functional>
#include <variant>

#include "general_types.hpp"
#include "string_utils.hpp"
#include "vector_utils.hpp"

#include "airborne_types.hpp"
#include "inverter.hpp"
#include "tdemsystem.hpp"
#include "aemsysteminversioninfo.hpp"
#include "samplebunch.hpp"
#include <Eigen/Cholesky>
#include <Eigen/LU>

inline void function_trace(const char* filename, const char* functionname, int linenumber) {
	std::cout << filename << " " << functionname << " " << linenumber << std::endl;
};
//#define function_trace function_trace(__FILE__, __FUNCTION__, __LINE__);
#define function_trace

namespace AEM::INVERTER::SBSINVERTER {

	static constexpr const char SCALEFACTOR[] = "ScaleFactor";

	class cOutputOptions {

	private:
		fs::path DumpBasePath;

	public:
		fs::path LogFile;
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

		fs::path DumpPath(const size_t datafilerecord, const size_t iteration) const {
			fs::path p = DumpBasePath;
			p += strprint("rec_%07d", (int)datafilerecord + 1);
			p += fs::path::preferred_separator;
			p += strprint("it_%03d", (int)iteration);
			p += fs::path::preferred_separator;
			return p;
		};

		cOutputOptions() {};

		cOutputOptions(const cBlock& b) {
			if (b.getvalue("LogFile", LogFile) == false) {
				std::string s = "A LogFile was not specified in the Output block ... so creating ./output.log.\n";
				glog.warningmsg(s);
				LogFile = "./output.log";
			}

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

			if (b.getvalue("Dump", Dump)) {
				if (Dump) {
					if (b.getvalue("DumpPath", DumpBasePath)) {
						DumpBasePath += fs::path::preferred_separator;
						DumpBasePath.make_preferred();
						makedirectory(DumpBasePath);
					}
					else {
						std::string s = "If Dump = yes a DumpPath must be specified.\n";
						glog.errormsg(_SRC_, s);
					}
				}
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
			: cConstraint(_key, _allowed_methods, _initials, _description) {
		};

		bool parse(const std::string& controlstring) {

			std::istringstream is(controlstring);
			std::string inputkey;
			is >> inputkey;
			if (key != inputkey) return false;

			if (alreadyparsed == true) {
				std::ostringstream msg;
				msg << "Constraint has already been set: " << std::endl << controlstring << std::endl;
				glog.errormsg(_SRC_, msg.str());
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
			: cConstraint(_key, _allowed_methods, _subscript, _description) {
		};

		bool parse(const std::string& controlstring) {

			std::istringstream is(controlstring);
			std::string inputkey;
			is >> inputkey;
			if (key != inputkey) return false;

			if (alreadyparsed == true) {
				std::ostringstream msg;
				msg << "Constraint has already been set: " << std::endl << controlstring << std::endl;
				glog.errormsg(_SRC_, msg.str());
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

	// A non-class memeber to get a stmfile path from a SBSInverter control file
	fs::path get_stmpath(const fs::path& controlfile) {
		cBlock b(controlfile);
		cBlock sb = b.findblock("EMSystem");
		std::string stmfile;
		if (sb.getvalue("SystemFile", stmfile)) {
			return fs::path(stmfile);
		}
		else {
			glog.errormsg(_SRC_, "No 'SystemFile' has been specified.\n");
		}
		return fs::path(stmfile);
	};

	template<typename AEMSystemClass, typename RT>
	class SBSInverter : public Inverter {

	public:

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

		SBSInverter(const fs::path& controlfile, const int& size, const int& rank, const bool& usingopenmp, const std::string commandline)
			: Inverter(controlfile, size, rank, usingopenmp, commandline)
		{
			try {
				initialise(controlfile);
			}
			catch (const std::string& msg) {
				std::cerr << msg;
				glog.logmsg(msg);
			}
			catch (const std::exception& e) {
				std::cerr << e.what();
				glog.logmsg(std::string(e.what()));
			}
		};

		std::string rec_it_str() const {
			std::ostringstream os;
			os << "Record " << Bunch.master_record() << " It " << CIS.iteration;
			return os.str();
		};

	private:
		inline static const size_t NCOMP = 3;
		inline static const size_t XCOMP = 0;
		inline static const size_t YCOMP = 1;
		inline static const size_t ZCOMP = 2;

		inline static const size_t xzcomp  = 3;
		inline static const size_t xyzcomp = 4;

		size_t nIndexStorageComponents() const {
			return 5;//For X Y Z XZ XYZ
		};

		

		std::vector<std::vector<std::vector<std::vector<int>>>> _vindex_;
		std::vector<size_t> _si_;
		std::vector<size_t> _sysi_;
		std::vector<size_t> _ci_;
		std::vector<size_t> _wi_;
		std::vector<size_t> _ei_;

		Vector Observed_unculled; //Non-culled observed data vector
		Vector Error_unculled;//Non-culled noise or error vector

		double ErrorAddition = 0.0;
		int    BeginGeometrySolveIteration = 0;
		bool   FreeGeometry = false;
		Matrix Wr;//Composite reference model matrix

		size_t StartRecord = 0; // 1-based first record to be inverted
		size_t EndRecord = std::numeric_limits<size_t>::max();// 1-based last record to be inverted
		size_t Subsample = 0;
		size_t nSoundings = 0;
		size_t nBunchSubsample = 0;
		size_t nDataPerSounding = 0;
		size_t nAllData = 0;
		size_t nLayers = 0;
		size_t nParamPerSounding = 0;
		size_t nGeomParamPerSounding = 0;
		size_t nGGAOffsetParamPerSounding = 0;
		size_t nScalingParam = 0;

		const size_t nSystems() const { return SI.size(); };

		size_t nPointsOutput = 0;
		std::vector<EarthStore> EStore;
		std::vector<GeometryStore> GStore;
		std::vector<std::vector<GGAOffsetStore>> GGAStore;

		cOutputOptions OutputOpt;
		using SysInfo = SystemInversionInfo<AEMSystemClass, RT>;
		using CompInfo = ComponentInversionInfo<RT>;
		std::vector<SysInfo> SI;

		//Column definitions
		Parameter pConductivity;
		Parameter pThickness;

		//using IFDKeyMap = std::map<std::string, cInvertibleFieldDefinition, caseinsensetiveless<std::string>>;
		using IFDKeyMap = std::map<std::string, Parameter, caseinsensetiveless<std::string>>;
		IFDKeyMap ifdGMap;
		std::vector<cFieldDefinition> fdPFRGvec;

		//Sample instances
		cSampleBunch Bunch;
		std::vector<SampleId> Id;
		std::vector<cKeyVec<std::string, cFdVrnt, caseinsensetiveequal<std::string>>> AncFld;

		// Cull data vector to only the active (non-null) data
		Vector cull(const Vector& vall) const {
			assert(ActiveData.size() == nData);
			assert(vall.size() == nAllData);
			Vector vcull(nData);
			for (size_t i = 0; i < nData; i++) {
				vcull[i] = vall[ActiveData[i]];
			}
			return vcull;
		};

		// Expand the culled data vector to full size filling with nulls
		Vector un_cull(const Vector& vculled) const {
			assert(ActiveData.size() == nData);
			assert(vculled.size() == nData);
			if (nData == nAllData) return vculled;

			const double ud = undefinedvalue<double>();
			Vector vall(nAllData);
			vall.setConstant(ud);
			for (size_t i = 0; i < nData; i++) {
				vall[ActiveData[i]] = vculled[i];
			}
			return vall;
		};

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

		AEM::SystemType aem_system_type() const {
			return AEM::aem_system_type<AEMSystemClass>();
		};

		const size_t& nsoundings() {
			return nSoundings;
		}

		cLinearConstraint LCrefc;//Conductivity reference model constraint
		cLinearConstraint LCreft;//Thickness reference model constraint
		cLinearConstraint LCrefg;//Geometry reference model constraint
		cLinearConstraint LCrefo;//GGAoffset reference model constraint
		cLinearConstraint LCrefs;//Scaling factor reference model constraint

		cLinearConstraint LCvcsmth;//Vertical conductivity smoothness constraint
		cLinearConstraint LCvcsim;//Vertical Conductivity Similarity

		cLinearConstraint LClatc;//Lateral conductivity constraint	
		cLinearConstraint LClatg;//Lateral geometry constraint

		cNonLinearConstraint NLCbounds;//Log-barrier bounds constraint
		cNonLinearConstraint NLCcablen;//Cable length constraint

		void loadcontrolfile(const fs::path& filename)
		{
			glog.logmsg(0, "Loading control file %s\n", filename.string().c_str());
			Control = cBlock(filename);
			cBlock ob = Control.findblock("Output");
			cBlock ib = Control.findblock("Input");

			OutputOpt = cOutputOptions(ob);
			Verbose = ob.getboolvalue("verbose");

			if (Rank == 0) {
				makedirectory_for(OutputOpt.LogFile);
			}

#ifdef ENABLE_MPI
			cMpiEnv::world_barrier();
#endif

			std::string suffix = stringvalue(Rank, ".%04d");
			OutputOpt.LogFile = insert_after_filename(OutputOpt.LogFile, suffix);
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
				makedirectory_for(OM->datafilename());
			}
#ifdef ENABLE_MPI
			cMpiEnv::world_barrier();
#endif

			OM->opendatafile(IM->datafilename(), IM->subsamplerate());
		}

		const CompInfo& compinfo(const size_t& systemindex, const size_t& componentindex) const {
			return SI[systemindex].CI[componentindex];
		};

		bool solve_thickness() const {
			return pThickness.solve();
		};

		bool solve_conductivity() const {
			return pConductivity.solve();
		};

		bool solve_geometry_elementname(const std::string& gname) const {
			return ifdGMap.at(gname).solve();
		};

		bool solve_geometry_index(const size_t index) const {
			return ifdGMap.at(TDEmGeometry::element_name(index)).solve();
		};

		bool solve_geometry() const {
			if (nGeomParamPerSounding > 0) return true;
			return false;
		};

		bool solve_ggaoffsets() const {
			if (nGGAOffsetParamPerSounding > 0) return true;
			return false;
		};
		
		bool solve_scalingfactor(const size_t sysi, const size_t ci) const {
			//if (nScalingParam == 0) return false;
			return SI[sysi].CI[ci].scalingfactor.solve();
		};

		bool solve_scalingfactors() const {
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
			s << " nF= " << nForwards / CIS.iteration;
			s << " nJ= " << nJacobians;
			return s.str();
		}

		fs::path dumppath() const {
			const size_t& record = Bunch.master_record();
			fs::path p = OutputOpt.DumpPath(record, CIS.iteration);
			return p;
		};

		void dump_record_number() {
			const size_t& record = Bunch.master_record();
			std::ofstream of(dumppath() + "record.dat");
			of << "Record\t" << record << std::endl;
		}

		int cindex(const size_t& si, const size_t& li) const {
			if (solve_conductivity() == false) {
				glog.errormsg(_SRC_, "Out of boundes in cindex().");
			}
			return (int)(si * nParamPerSounding + pConductivity.get_poffset() + li);
			//return (int)(si * nParamPerSounding + cOffset + li);
		}

		int tindex(const size_t& si, const size_t& li) const {
			if (solve_thickness() == false) {
				glog.errormsg(_SRC_, "Out of boundes in tindex().");
			}
			return (int)(si * nParamPerSounding + pThickness.get_poffset() + li);
			//return (int)(si * nParamPerSounding + tOffset + li);
		}

		int gindex(const size_t& si, const std::string& gname) const {
			const Parameter& p = ifdGMap.at(gname);
			int poffset = p.get_poffset();
			if (poffset >= 0) return (int)(si * nParamPerSounding + poffset);
			return -1;
		};

		int gindex(const size_t& si, const size_t& gi) const {
			const std::string gname = TDEmGeometry::element_name(gi);
			return gindex(si,gname);
		};

		int ggaoffsetindex(const size_t& si, const size_t& sysi, const size_t& ci) const {
			const Parameter& p = SI[sysi].CI[ci].ggaoffset;
			int po = p.get_poffset();
			if (po >= 0) {
				const int pi = (int)(si * nParamPerSounding + po);
				return pi;
			}
			return -1;
		};

		int sfindex(const size_t& sysi, const size_t& ci) const {
			const Parameter& p = SI[sysi].CI[ci].scalingfactor;
			int po = p.get_poffset();
			return po;
		};

		void openlogfile() {
			glog.logmsg(0, "Opening log file %s\n", OutputOpt.LogFile.string().c_str());
			glog.open(OutputOpt.LogFile);
			glog.logmsg(0, "%s\n", CommandLine.c_str());
			glog.logmsg(0, "%s\n", versionstring(GAAEM_VERSION, __TIME__, __DATE__).c_str());
			glog.logmsg(0, "Working directory %s\n", getcurrentdirectory().c_str());
			if (UsingOpenMP && Size > 1) {
				glog.logmsg(0, "Using OpenMP threading Processes=%d\tRank=%d\n", Size, Rank);
			}
			else if (Size > 1) {
				glog.logmsg(0, "Using MPI Processes=%d\tRank=%d\n", Size, Rank);
			}
			else {
				glog.logmsg(0, "Standalone Processes=%d\tRank=%d\n", Size, Rank);
			}
			glog.logmsg(0, "Control file %s\n", Control.Filename.c_str());
			glog.log_to_file(Control.get_as_string());
			glog.flush();
		};

		void parse_options() {
			cBlock b = Control.findblock("Options");
			if (b.getvalue("StartRecord", StartRecord) == false) {
				StartRecord = 1;
			}

			if (b.getvalue("EndRecord", EndRecord) == false) {
				EndRecord = std::numeric_limits<size_t>::max();
			}

			if (b.getvalue("Subsample", Subsample) == false) {
				Subsample = 1;
			}

			if (b.getvalue("SoundingsPerBunch", nSoundings) == false) {
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

			NormType = NormType::L2;//default
			std::string nt = b.getstringvalue("NormType");
			if (!isdefined(nt)) {
				NormType = NormType::L2;
			}
			else if (strcasecmp(nt, "L1") == 0) {
				NormType = NormType::L1;
			}
			else if (strcasecmp(nt, "L2") == 0) {
				NormType = NormType::L2;
			}
			else {
				glog.errormsg(_SRC_, "Unknown NormType %s.", nt.c_str());
			}

			MaxIterations = b.getsizetvalue("MaximumIterations");
			MinimumPhiD = b.getdoublevalue("MinimumPhiD");
			MinimumImprovement = b.getdoublevalue("MinimumPercentageImprovement");
		};

		void parse_constraints(const cBlock& b) {

			LCrefc = cLinearConstraint("ConductivityReferenceModel", { }, "RefCon", "Conductivity reference model");
			LCreft = cLinearConstraint("ThicknessReferenceModel", { }, "RefThk", "Thickness reference model");
			LCrefg = cLinearConstraint("GeometryReferenceModel", { }, "RefGeom", "Geometry reference model");
			LCrefo = cLinearConstraint("GGAOffsetReferenceModel", { }, "RefGGAOffset", "GGA offset reference model");
			LCrefs = cLinearConstraint("ScalingFactorsReferenceModel", { }, "RefScalingFactors", "ScalingFactors reference model");

			LCvcsmth = cLinearConstraint("VerticalConductivity", { "Minimise1stDerivatives","Minimise2ndDerivatives" }, "VCsmth", "Vertical conductivity smoothness");
			LCvcsim = cLinearConstraint("VerticalConductivitySimilarity", { }, "VCsim", "Vertical conductivity similarity");
			LClatc = cLinearConstraint("LateralConductivity", { "Minimise1stDerivatives", "Minimise2ndDerivatives","Similarity" }, "LatCon", "Lateral conductivity smoothness");
			LClatg = cLinearConstraint("LateralGeometry", { "Similarity","MinimiseAccelerations","MinimiseAccelerationDerivatives","Minimise2ndDerivativesOfDifferenceFromReferenceModel" }, "LatGeom", "Lateral geometry smoothness");
			NLCcablen = cNonLinearConstraint("CableLength", { "Input","InputBunchMean","BunchSimilarity" }, "CabLength", "Cable length");
			NLCbounds = cNonLinearConstraint("Bounds", { }, "Bounds", "Log-barrier Bounds");

			for (size_t i = 0; i < b.Entries.size(); i++) {
				const std::string cstr = b.Entries[i];
				std::string key = cConstraint::get_key(cstr);
				if (key[0] == '/')continue;
				if (LCrefc.parse(cstr)) {
					LCrefc.set_operates_on_difference_from_reference_model();
				}
				else if (LCreft.parse(cstr)) {
					LCreft.set_operates_on_difference_from_reference_model();
				}
				else if (LCrefg.parse(cstr)) {
					LCrefg.set_operates_on_difference_from_reference_model();
				}
				else if (LCrefo.parse(cstr)) {
					LCrefo.set_operates_on_difference_from_reference_model();
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
				else if (NLCcablen.parse(cstr)) {
					//nothing to do
				}
				else if (NLCbounds.parse(cstr)) {
					//nothing to do
				}
				else {
					std::stringstream msg;
					msg << "Unknown constraint " << key << std::endl;
					glog.errormsg(_SRC_, msg.str());
				}
			}
			if (NLCbounds.alreadyparsed == false) {
				NLCbounds.alpha = 1.0;
			}
		};

		void set_field_definitions() {
			cBlock b = Control.findblock("Input.AncillaryFields");
			set_field_definitions_ancillary(b);

			if (nSoundings > 1) {
				if (AncFld[0].keyindex("line") < 0) {
					glog.errormsg(_SRC_, "Must specify a linenumber field.");
				}
			}

			b = Control.findblock("Input.Geometry");
			ifdGMap = get_field_definitions_geometry(b);
			fdPFRGvec = get_pfr_geometry_field_definitions(b, ifdGMap);

			b = Control.findblock("Input.Earth");
			bool status = b.getvalue("NumberOfLayers", nLayers);
			if (status == false) {
				std::stringstream msg;
				msg << "The NumberOfLayers must be specified in Input.Earth\n";
				glog.errormsg(_SRC_, msg.str());
			}

			pConductivity = Parameter(b, "Conductivity");
			if (nLayers > 1) {
				pThickness = Parameter(b, "Thickness");
			}
		};

		void set_field_definitions_ancillary(const cBlock& parent) {
			AncFld.resize(nSoundings);
			glog.logmsg(0, "Adding ancillary fields\n");
			const cBlock& b = parent;
			for (size_t i = 0; i < b.Entries.size(); i++) {
				std::string key = b.key(i);
				std::string value = b.value(i);
				if (key[0] == '/') continue;//Skip commented out fields
				if (key.size() == 0) continue;//Skip blank lines			
				glog.logmsg(0, "\tAdding ancillary field %s\n", value.c_str());
				cFieldDefinition fd(parent, key);
				cFdVrnt fdvrnt(fd, cVrnt());
				IM->set_variant_type(fd, fdvrnt.vnt);
				for (size_t si = 0; si < nSoundings; si++) {
					AncFld[si].add(key, fdvrnt);
				}
			}
		};

		IFDKeyMap get_field_definitions_geometry(const cBlock& parentblock) const {
			IFDKeyMap gmap;
			for (size_t i = 0; i < TDEmGeometry::NELEM; i++) {
				std::string gname = TDEmGeometry::element_name(i);
				if (gmap.find(gname) == gmap.end()) {
					gmap[gname] = Parameter(parentblock, gname);
				}
				else {
					std::string msg = strprint("Parameter %s has already been already added.", gname.c_str());
					glog.errormsg(_SRC_, msg);
				}
			}
			return gmap;
		};

		bool reconstruct_primary() const {
			for (size_t sysi = 0; sysi < SI.size(); sysi++) {
				if (SI[sysi].ReconstructPrimary == true) return true;
			}
			return false;
		};

		std::vector<cFieldDefinition> get_pfr_geometry_field_definitions(const cBlock& parent, const IFDKeyMap& gmap) const {
			std::vector<cFieldDefinition> fdvec;
			if (reconstruct_primary() == false) return fdvec;
			const size_t ng = TDEmGeometry::NELEM;
			fdvec.resize(ng);
			for (size_t gi = 0; gi < ng; gi++) {
				std::string gname = TDEmGeometry::element_name(gi);
				cBlock b = parent.findblock(gname);
				if (b.empty() == false) {
					//Check if there is a PFR field defined and if so define it
					cFieldDefinition fdtfr(b, "TFR");
					if (fdtfr.isinitialised()) {
						std::ostringstream msg;
						msg << "'TFR' is deprecated for specifing the primary-field-reconstruction geometry fields.	Please use 'PFR' instead.";
						glog.warningmsg(msg.str());
						fdvec[gi] = fdtfr;
					}

					cFieldDefinition fdpfr(b, "PFR");
					if (fdpfr.isinitialised()) {
						fdvec[gi] = fdpfr;
					}
					// Otherwiae leave it undefined and it will be defined to the input in read_geometry()
				}
			}
			return fdvec;
		};

		void setup_parameters() {
			Id.resize(nSoundings);
			EStore.resize(nSoundings);
			GStore.resize(nSoundings);

			nParamPerSounding = 0;
			nGeomParamPerSounding = 0;
			nGGAOffsetParamPerSounding = 0;

			if (solve_conductivity()) {
				pConductivity.set_poffset((int)nParamPerSounding);
				nParamPerSounding += nLayers;
			}

			if (solve_thickness()) {
				pThickness.set_poffset((int)nParamPerSounding);
				nParamPerSounding += nLayers - 1;
			}

			//Geometry params
			for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
				std::string gname = TDEmGeometry::element_name(gi);
				Parameter& p = ifdGMap.at(gname);
				if (p.solve()) {
					p.set_poffset((int)nParamPerSounding);
					nParamPerSounding++;
					nGeomParamPerSounding++;
				}
			};

			//GGA Offset
			for (size_t sysi = 0; sysi < SI.size(); sysi++) {
				for (size_t ci = 0; ci < NCOMP; ci++) {
					CompInfo& C = SI[sysi].CI[ci];
					Parameter& p = C.ggaoffset;
					if (p.solve()) {
						p.set_poffset((int)nParamPerSounding);
						nParamPerSounding++;
						nGGAOffsetParamPerSounding++;
					}
				}
			}

			//Scaling params
			for (size_t sysi = 0; sysi < SI.size(); sysi++) {
				for (size_t ci = 0; ci < NCOMP; ci++) {
					CompInfo& C = SI[sysi].CI[ci];
					Parameter& p = C.scalingfactor;
					if (p.solve()) {
						p.set_poffset((int)(nParamPerSounding * nSoundings + nScalingParam));
						nScalingParam++;
					}
				}
			}

			nParam = nParamPerSounding * nSoundings + nScalingParam;
			RefParam.resize(nParam);
			RefParamStd.resize(nParam);

			if (nGeomParamPerSounding == 0) {
				LCrefg.alpha = 0.0;
				NLCcablen.alpha = 0.0;
			}

			
			if (nGGAOffsetParamPerSounding == 0) {
				LCrefo.alpha = 0.0;
			}

			if (nScalingParam == 0) {
				LCrefs.alpha = 0.0;
			}

			if (nSoundings == 1) {
				LClatc.alpha = 0.0;
				LClatg.alpha = 0.0;
			}
		};

		void setup_parameter_bounds() {
			function_trace;
			const static double ud = undefinedvalue<double>();
			Param_Min.resize(nParam);
			Param_Max.resize(nParam);
			Param_Min.setConstant(ud);
			Param_Max.setConstant(ud);

			if (pConductivity.bound()) {
				for (size_t si = 0; si < nSoundings; si++) {
					const EarthStore& e = EStore[si];
					for (size_t li = 0; li < nLayers; li++) {
						const int pi = cindex(si, li);
						Param_Min[pi] = std::log10(e.minval.conductivity[li]);
						Param_Max[pi] = std::log10(e.maxval.conductivity[li]);
					}
				}
			}

			if (pThickness.bound()) {
				for (size_t si = 0; si < nSoundings; si++) {
					const EarthStore& e = EStore[si];
					for (size_t li = 0; li < nLayers - 1; li++) {
						const int pi = tindex(si, li);
						Param_Min[pi] = std::log10(e.minval.thickness[li]);
						Param_Max[pi] = std::log10(e.maxval.thickness[li]);
					}
				}
			}

			for (size_t si = 0; si < nSoundings; si++) {
				const GeometryStore& g = GStore[si];
				for (size_t i = 0; i < TDEmGeometry::NELEM; i++) {
					const std::string ename = TDEmGeometry::element_name(i);
					const Parameter& p = ifdGMap.at(ename);
					if (p.bound()) {
						const int pi = gindex(si, ename);
						Param_Min[pi] = g.minval[ename];
						Param_Max[pi] = g.maxval[ename];
					}
				}
			}

			for (size_t si = 0; si < nSoundings; si++) {
				for (size_t sysi = 0; sysi < SI.size(); sysi++) {
					for (size_t ci = 0; ci < NCOMP; ci++) {
						const CompInfo& C = SI[sysi].CI[ci];
						const Parameter& p = C.ggaoffset;
						if (p.bound()) {
							int pi = ggaoffsetindex(si, sysi, ci);
							Param_Min[pi] = C.ggaoffsetStore.minval[si];
							Param_Max[pi] = C.ggaoffsetStore.maxval[si];
						}
					}
				}
			}

			for (size_t sysi = 0; sysi < SI.size(); sysi++) {
				for (size_t ci = 0; ci < NCOMP; ci++) {
					const CompInfo& C  = SI[sysi].CI[ci];
					const Parameter& p = C.scalingfactor;
					if (p.bound()) {
						int pi = p.get_poffset();
						const CompInfo& C = SI[sysi].CI[ci];
						Param_Min[pi] = C.sfStore.minval;
						Param_Max[pi] = C.sfStore.maxval;
					}
				}
			}
		};

		void initialise_Wc() {
			cLinearConstraint& C = LCrefc;
			C.W = Matrix::Zero(nParam, nParam);
			if (solve_conductivity() == false)return;

			for (size_t si = 0; si < nSoundings; si++) {
				const EarthStore& e = EStore[si];
				std::vector<double> t(nLayers);
				if (nLayers == 1) {
					t[0] = 1;
				}
				else if (nLayers == 2) {
					t[0] = e.refval.thickness[0];
					t[1] = e.refval.thickness[0];
				}
				else {
					for (size_t i = 0; i < (nLayers - 1); i++) {
						t[i] = e.refval.thickness[i];
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
		};

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
		};

		void initialise_Wg() {
			cLinearConstraint& C = LCrefg;
			C.W = Matrix::Zero(nParam, nParam);
			if (solve_geometry() == false)return;

			double s = C.alpha / (double)(nGeomParamPerSounding * nSoundings);
			for (size_t si = 0; si < nSoundings; si++) {
				for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
					const int pi = gindex(si, gi);
					if (pi >= 0) {
						C.W(pi, pi) = s / (RefParamStd[pi] * RefParamStd[pi]);
					}
				}
			}
		};

		void initialise_Wo() {
			cLinearConstraint& C = LCrefo;
			C.W = Matrix::Zero(nParam, nParam);
			if (solve_ggaoffsets() == false)return;

			double s = C.alpha / (double)(nGGAOffsetParamPerSounding * nSoundings);
			for (size_t si = 0; si < nSoundings; si++) {
				for (size_t sysi = 0; sysi < SI.size(); sysi++) {
					for (size_t ci = 0; ci < NCOMP; ci++) {
						if (compinfo(sysi,ci).Use){
							const int pi = ggaoffsetindex(si, sysi, ci);
							if (pi >= 0) {
								C.W(pi, pi) = s / (RefParamStd[pi] * RefParamStd[pi]);
							}
						}
					}
				}
			}
		};

		void initialise_Ws() {
			cLinearConstraint& C = LCrefs;
			C.W = Matrix::Zero(nParam, nParam);
			if (solve_scalingfactors() == false)return;

			double s = C.alpha / (double)(nScalingParam);
			for (size_t sysi = 0; sysi < SI.size(); sysi++) {
				for (size_t ci = 0; ci < NCOMP; ci++) {
					if (compinfo(sysi, ci).Use) {
						const int pi = sfindex(sysi, ci);
						if (pi >= 0) {
							C.W(pi, pi) = s / (RefParamStd[pi] * RefParamStd[pi]);
						}
					}
				}
			}
		};
		
		void initialise_VC() {
			cLinearConstraint& C = LCvcsmth;
			C.W = Matrix::Zero(nParam, nParam);
			if (solve_conductivity() == false) return;
			if (C.alpha == 0) return;
			if (C.method == "Minimise1stDerivatives") {
				initialise_VC_1st_derivative(C);
			}
			else if (C.method == "Minimise2ndDerivatives") {
				initialise_VC_2nd_derivative(C);
			}
		};

		void initialise_VC_1st_derivative(cLinearConstraint& C) {
			if (nLayers < 3) return;
			Matrix L = Matrix::Zero(nSoundings * (nLayers - 1), nParam);
			size_t nrows = 0;
			for (size_t si = 0; si < nSoundings; si++) {
				const EarthStore& e = EStore[si];
				std::vector<double> t = e.refval.dummy_thickness();
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
		};

		void initialise_VC_2nd_derivative(cLinearConstraint& C) {
			if (nLayers < 3) return;
			Matrix L = Matrix::Zero(nSoundings * nLayers, nParam);
			size_t nrows = 0;
			for (size_t si = 0; si < nSoundings; si++) {
				const EarthStore& e = EStore[si];
				std::vector<double> t = e.refval.dummy_thickness();
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
				L(nrows, pi1) = 1.0 * s / t[0];
				nrows++;

				//Minimise 1st deriv on last interface
				pi0 = cindex(si, nLayers - 2);
				pi1 = cindex(si, nLayers - 1);
				s = std::sqrt(t[nLayers - 2] / tavg);
				L(nrows, pi0) = -1.0 * s / t[nLayers - 1];
				L(nrows, pi1) = 1.0 * s / t[nLayers - 1];
				nrows++;
			}
			C.W = L.transpose() * L;
			C.W *= (C.alpha / (double)(nrows));
		};

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
					L(nrows, psi) = 1.0 / std;
					for (size_t ri = 0; ri < nSoundings; ri++) {
						if (ri != si) {
							const int pri = cindex(ri, li);
							L(nrows, pri) = -1.0 / (std * (double)(nSoundings - 1));
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
				const EarthStore& e = EStore[si];
				std::vector<double> t = e.refval.dummy_thickness();
				double tavg = mean(t);
				//Loop over constraints equations
				for (size_t li = 0; li < nLayers; li++) {
					//double s = std::sqrt(t[li] / tavg);//sqrt because it gets squared in L'L									
					const int lpindex = cindex(si, li);
					L(nrows, lpindex) = 1.0;
					//Loop over layers for this equation
					for (size_t ki = 0; ki < nLayers; ki++) {
						const int kpindex = cindex(si, ki);
						if (li != ki) {
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
			for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
				if (solve_geometry_index(gi) == false)continue;
				for (size_t si = 1; si < nSoundings - 1; si++) {
					double d01 = std::hypot(Id[si].x - Id[si - 1].x, Id[si].y - Id[si - 1].y);
					double d12 = std::hypot(Id[si].x - Id[si + 1].x, Id[si].y - Id[si + 1].y);

					const int pi0 = gindex(si - 1, gi);
					const int pi1 = gindex(si, gi);
					const int pi2 = gindex(si + 1, gi);
					const double std = RefParamStd[pi1];
					L(nrows, pi0) = 1.0 / (d01 * std);
					L(nrows, pi1) = -1.0 / (d01 * std) - 1.0 / (d12 * std);
					L(nrows, pi2) = 1.0 / (d12 * std);
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
			for (size_t j = 0; j < nSoundings - 1; j++) {
				d += std::hypot(Id[j].x - Id[j + 1].x, Id[j].y - Id[j + 1].y);
			}
			d = d / (double)(nSoundings - 1);//average sample distance

			for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
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
			for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
				if (solve_geometry_index(gi) == false)continue;
				for (size_t si = 0; si < nSoundings; si++) {
					const int spi = gindex(si, gi);
					const double std = RefParamStd[spi];
					L(nrows, spi) = 1.0 / std;
					for (size_t ri = 0; ri < nSoundings; ri++) {
						const int rpi = gindex(ri, gi);
						if (si != ri) {
							L(nrows, rpi) = -1.0 / (std * (double)(nSoundings - 1));
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
					C.data[si] = GStore[si].input.txrx_dr();
				}
				else if (NLCcablen.method == "InputBunchMean") {
					C.data[si] = GStore[si].input.txrx_dr();
				}
				else if (NLCcablen.method == "BunchSimilarity") {
					C.data[si] = 0.0;
				}
				else {
					glog.errormsg(_SRC_, "Unknown NLCcablen.method.");
				}
				C.W(si, si) = s / (C.err[si] * C.err[si]);
			}

			if (NLCcablen.method == "InputBunchMean") {
				double mn = C.data.mean();
				for (size_t si = 0; si < nSoundings; si++) {
					C.data[si] = mn;
				}
			}
		}

		Vector CableLengths(const Vector& m) const {
			Vector cablelength(nSoundings);
			const std::vector<TDEmGeometry> gv = get_geometry(m);
			for (size_t si = 0; si < nSoundings; si++) {
				const TDEmGeometry& g = gv[si];
				const double dr = g.txrx_dr();
				cablelength[si] = dr;
			}
			return cablelength;
		}

		Vector CableLengthConstraint_forward(const Vector& m) const {
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

			const std::vector<TDEmGeometry> gv = get_geometry(m);
			double s = C.alpha / (double)(nSoundings);
			for (size_t si = 0; si < nSoundings; si++) {
				const TDEmGeometry& g = gv[si];
				const double dr = g.txrx_dr();
				const int pix = gindex(si, "txrx_dx");
				const int piy = gindex(si, "txrx_dy");
				const int piz = gindex(si, "txrx_dz");

				double f = 1.0;
				if (NLCcablen.method == "BunchSimilarity") {
					f = (double)(nSoundings - 1) / (double)nSoundings;
				}

				if (pix >= 0) {
					C.J(si, pix) = f * g.txrx_dx() / dr;
				}
				if (piy >= 0) {
					C.J(si, piy) = f * g.txrx_dy() / dr;
				}
				if (piz >= 0) {
					C.J(si, piz) = f * g.txrx_dz() / dr;
				}
			}
		}

		void initialise_BoundsConstraint() {
			cNonLinearConstraint& C = NLCbounds;
			C.W = Matrix::Zero(nParam, nParam);
			C.J = Matrix::Zero(nParam, nParam);
			C.err = Vector::Zero(nParam);
			C.data = Vector::Zero(nParam);
			if (C.alpha == 0.0) return;

			double s = C.alpha / (double)(nParam);
			C.data.resize(nParam);
			for (size_t pi = 0; pi < nParam; pi++) {
				C.err[pi] = 1.0;
				C.data[pi] = 0.0;
				C.W(pi, pi) = s / (C.err[pi] * C.err[pi]);
			}
		}

		static double log_barrier(const double& L, const double& U, const double& n, double& x) {
			//Make sure calculation is not subject to log(0) or divide by 0
			if (x <= L) x = L + (U - L) / 1000.0;
			else if (x >= U) x = U - (U - L) / 1000.0;

			const double v = -std::log(std::pow(((U - x) / (U - L)), n))
				- std::log(std::pow(((x - L) / (U - L)), n))
				+ 2.0 * std::log(1.0 / std::pow(2.0, n));
			return v;
		}

		static double log_barrier_deriv(const double& L, const double& U, const double& n, double x) {
			//Make sure calculation is not subject to log(0) or divide by 0
			if (x <= L) x = L + (U - L) / 1000.0;
			else if (x >= U) x = U - (U - L) / 1000.0;

			const double dv = (n * (L + U - 2.0 * x)) / ((L - x) * (U - x));
			return dv;
		}

		bool is_bound(const size_t pindex) const {
			if (isdefined<double>(Param_Min[pindex]) == false) return false;
			if (isdefined<double>(Param_Max[pindex]) == false) return false;
			return true;
		};

		Vector BoundsConstraint_forward(const Vector& m) const {			
			const cNonLinearConstraint& C = NLCbounds;
			Vector predicted = Vector::Zero(nParam);
			if (C.alpha == 0.0) return predicted;
			for (size_t pi = 0; pi < nParam; pi++) {				
				const double& L = Param_Min[pi];
				const double& U = Param_Max[pi];
				const double& N = 0.5;
				double x = m[pi];
				if(is_bound(pi)) {
					predicted[pi] = log_barrier(L, U, N, x);
				}
				else predicted[pi] = 0.0;
			}
			return predicted;
		}

		void BoundsConstraint_jacobian(const Vector& m) {
			cNonLinearConstraint& C = NLCbounds;
			if (C.alpha == 0.0) return;
			double s = C.alpha / (double)(nParam);
			for (size_t pi = 0; pi < nParam; pi++) {
				const double& L = Param_Min[pi];
				const double& U = Param_Max[pi];
				const double N = 0.5;
				double x = m[pi];
				if (is_bound(pi)) {
					C.J(pi, pi) = log_barrier_deriv(L, U, N, x);
				}
				else C.J(pi, pi) = 0.0;
			}
		}

		void initialise_Wr() {
			initialise_Wc();
			initialise_Wt();
			initialise_Wg();
			initialise_Wo();
			initialise_Ws();

			Wr = Matrix::Zero(nParam, nParam);
			if (LCrefc.alpha > 0.0) Wr += LCrefc.W;
			if (LCreft.alpha > 0.0) Wr += LCreft.W;
			if (LCrefg.alpha > 0.0) Wr += LCrefg.W;
			if (LCrefo.alpha > 0.0) Wr += LCrefo.W;
			if (LCrefs.alpha > 0.0) Wr += LCrefs.W;
		}

		void initialise_Wm() {
			initialise_Wr();
			initialise_VC();
			initialise_LC();
			initialise_LG();
			initialise_QC(LCvcsim);
			initialise_CableLengthConstraint();
			initialise_BoundsConstraint();
			Wm = Wr + LCvcsmth.W + LCvcsim.W + LClatc.W + LClatg.W;
		}

		void dump_W_matrices() {
			if (OutputOpt.Dump) {
				const std::string dp = dumppath().string();
				makedirectory(dp);
				writetofile(Wd, dp + "Wd.dat");
				writetofile(Wr, dp + "Wr.dat");

				LCrefc.write_W_matrix(dp);
				LCreft.write_W_matrix(dp);
				LCrefg.write_W_matrix(dp);
				LCrefo.write_W_matrix(dp);
				LCrefs.write_W_matrix(dp);
				LCvcsmth.write_W_matrix(dp);
				LCvcsim.write_W_matrix(dp);
				LClatc.write_W_matrix(dp);
				LClatg.write_W_matrix(dp);
				NLCcablen.write_W_matrix(dp);
				NLCbounds.write_W_matrix(dp);
				writetofile(Wm, dp + "Wm.dat");
			}
		}

		// Value index not data index 
		const int& vindex(const size_t& sampleindex, const size_t& systemindex, const size_t& componentindex, const size_t& windowindex) {
			const size_t& si = sampleindex;
			const size_t& sysi = systemindex;
			const size_t& ci = componentindex;
			const size_t& wi = windowindex;
			return _vindex_[si][sysi][ci][wi];
		};

		void set_fftw_lock() {
			#if defined _OPENMP
			//If OpenMP is being used set the thread lock while FFTW initialises
			if (UsingOpenMP) {
				omp_set_lock(&fftw_thread_lock);
			}
			#endif
		};

		void unset_fftw_lock() {
			#if defined _OPENMP
			if (UsingOpenMP) {
				omp_unset_lock(&fftw_thread_lock);
			}
			#endif
		};

		void initialise_systems() {
			if (aem_system_type() == AEM::SystemType::TimeDomain) {
				set_fftw_lock();
			}
			std::vector<cBlock> B = Control.findblocks("EMSystem");
			const size_t nsys = B.size();
			for (size_t sysi = 0; sysi < nsys; sysi++) {
				cBlock& b = B[sysi];
				fs::path stmfile;
				if (b.getvalue("SystemFile", stmfile)) {
					glog.logmsg(0, "Reading AEM system file %s\n", stmfile.string().c_str());
					SI.emplace_back(SystemInversionInfo<AEMSystemClass, RT>(b, nSoundings));
					SI[sysi].set_units(IM.get());
				}
				else glog.errormsg(_SRC_, "No AWM 'SystemFile' is specified.\n");
			}
			if (aem_system_type() == AEM::SystemType::TimeDomain) {
				unset_fftw_lock();
			}
		}

		void resize_reverse_vindex_arrays(const size_t& nalldata) {
			const size_t vsize = value_size<RT>();
			size_t n = nalldata / vsize;
			_si_.resize(n);
			_sysi_.resize(n);
			_ci_.resize(n);
			_wi_.resize(n);
		};

		void set_reverse_vindex(const size_t& vi, const size_t& si, const size_t& sysi, const size_t& ci, const size_t& wi) {
			_si_[vi]   = si;
			_sysi_[vi] = sysi;
			_ci_[vi]   = ci;
			_wi_[vi]   = wi;
		};

		bool unused_or_composite_component(const SysInfo& S, const size_t& ci) const {
			if (S.CI[ci].Use == false) return true;

			//Belongs to XZ or XYZ
			if (S.InvertXYZAmplitude) return true;
			else if (S.InvertXZAmplitude && ci == XCOMP) return true;
			else if (S.InvertXZAmplitude && ci == ZCOMP) return true;
			else return false;
		};

		void setup_data() {
			const size_t nc = nIndexStorageComponents();//5 because of possible xz and xyz inversions
			nAllData = 0;
			const size_t nsys = nSystems();
			const size_t vsize = value_size<RT>();
			_vindex_.resize(nSoundings);
			for (size_t si = 0; si < nSoundings; si++) {
				_vindex_[si].resize(nsys);
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					_vindex_[si][sysi].resize(nc);
					for (size_t ci = 0; ci < nc; ci++) {
						const size_t nw = SI[sysi].nWindows();
						_vindex_[si][sysi][ci].resize(nw);
						for (size_t wi = 0; wi < nw; wi++) {
							_vindex_[si][sysi][ci][wi] = -1;
						}
					}
				}
			}

			// First count the data
			for (size_t si = 0; si < nSoundings; si++) {
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					SysInfo& S = SI[sysi];
					if (S.InvertXYZAmplitude) {
						nAllData += S.nWindows() * vsize;
					}
					else if (S.InvertXZAmplitude) {
						nAllData += S.nWindows() * vsize;
					}

					for (size_t ci = 0; ci < NCOMP; ci++) {
						if (unused_or_composite_component(S, ci)) continue;
						CompInfo& C = S.CI[ci];
						if (C.Use) nAllData += S.nWindows() * vsize;
					}
				}
			}
			resize_reverse_vindex_arrays(nAllData);
			size_t firstcount = nAllData;

			int vi = 0;
			nAllData = 0;//Recount and assert later
			for (size_t si = 0; si < nSoundings; si++) {
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					SysInfo& S = SI[sysi];
					if (S.InvertXYZAmplitude) {
						nAllData += S.nWindows() * vsize;
						for (size_t wi = 0; wi < S.nWindows(); wi++) {
							_vindex_[si][sysi][xyzcomp][wi] = vi;
							set_reverse_vindex(vi,si,sysi,xyzcomp,wi);
							vi++;
						}
					}
					else if (S.InvertXZAmplitude) {
						nAllData += S.nWindows() * vsize;
						for (size_t wi = 0; wi < S.nWindows(); wi++) {
							_vindex_[si][sysi][xzcomp][wi] = vi;
							set_reverse_vindex(vi, si, sysi, xzcomp, wi);
							vi++;
						}
					}

					for (size_t ci = 0; ci < NCOMP; ci++) {
						if (unused_or_composite_component(S, ci)) continue;
						CompInfo& C = S.CI[ci];
						nAllData += S.nWindows() * vsize;
						for (size_t wi = 0; wi < S.nWindows(); wi++) {
							_vindex_[si][sysi][ci][wi] = vi;
							set_reverse_vindex(vi, si, sysi, ci, wi);
							vi++;
						}
					}
				}
			}
			assert(firstcount == nAllData);
		}

		double combine_noises(const double x, const double ex, const double z, const double ez) {
			double e;
			double r = std::hypot(x, z);
			if (r == 0) e = std::hypot(ex, ez);
			else e = std::hypot((x/r) * ex, (z/r) * ez);
			return e;
		};

		double combine_noises(const double x, const double ex, const double y, const double ey, const double z, const double ez) {
			double e;
			double r = std::hypot(x, y, z);
			if (r == 0) e = std::hypot(ex, ey, ez);
			else e = std::hypot((x / r) * ex, (y / r) * ey, (z / r) * ez);
			return e;
		};

		cdouble combine_noises(const cdouble x, const cdouble xe, const cdouble z, const cdouble ze) {
			const double& er = combine_noises(x.real(), xe.real(), z.real(), ze.real());
			const double& ei = combine_noises(x.imag(), xe.imag(), z.imag(), ze.imag());
			return cdouble(er,ei);
		};

		cdouble combine_noises(const cdouble x, const cdouble xe, const cdouble y, const cdouble ye, const cdouble z, const cdouble ze) {
			const double& er = combine_noises(x.real(), xe.real(), y.real(), ye.real(), z.real(), ze.real());
			const double& ei = combine_noises(x.imag(), xe.imag(), y.imag(), ye.imag(), z.imag(), ze.imag());
			return cdouble(er, ei);
		};

		bool initialise_bunch_data() {
			Observed_unculled.resize(nAllData);
			Error_unculled.resize(nAllData);
			
			const size_t nsys = nSystems();
			for (size_t si = 0; si < nSoundings; si++) {
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					SysInfo& S = SI[sysi];
					if (S.ReconstructPrimary) {
						TDEmVectorResponse<RT> P = S.System->forward_model_primary_field(GStore[si].pfr);
						if (S.CI[XCOMP].Use) S.CI[XCOMP].data[si].P = P[XCOMP];
						if (S.CI[YCOMP].Use) S.CI[YCOMP].data[si].P = P[YCOMP];
						if (S.CI[ZCOMP].Use) S.CI[ZCOMP].data[si].P = P[ZCOMP];
					}

					if (S.InvertXYZAmplitude) {
						for (size_t wi = 0; wi < S.nWindows(); wi++) {
							int vi = vindex(si, sysi, xyzcomp, wi);
							RT X, Z, Y;

							if (S.InvertPSI) {
								//PSI is stored in Total
								X = S.CI[XCOMP].data[si].T[wi];
								Y = S.CI[YCOMP].data[si].T[wi];
								Z = S.CI[ZCOMP].data[si].T[wi];
							}
							else if (S.InvertTotalField) {
								X = S.CI[XCOMP].data[si].T[wi];
								Y = S.CI[YCOMP].data[si].T[wi];
								Z = S.CI[ZCOMP].data[si].T[wi];
							}
							else {
								X = S.CI[XCOMP].data[si].S[wi];
								Y = S.CI[YCOMP].data[si].S[wi];
								Z = S.CI[ZCOMP].data[si].S[wi];
							}

							RT ov = AEM::ewise_hypot(X, Y, Z);
							set_data_vector(Observed_unculled, vi, ov);

							const RT& Xerr = S.CI[XCOMP].data[si].E[wi];
							const RT& Yerr = S.CI[YCOMP].data[si].E[wi];
							const RT& Zerr = S.CI[ZCOMP].data[si].E[wi];
							RT ev = combine_noises(X, Xerr, Y, Yerr, Z, Zerr);
							set_data_vector(Error_unculled, vi, ev);
						}
					}
					else if (S.InvertXZAmplitude) {
						for (size_t wi = 0; wi < S.nWindows(); wi++) {
							int vi = vindex(si, sysi, xzcomp, wi);
							RT X,Z;

							if (S.InvertPSI) {
								//PSI is stored in Total
								X = S.CI[XCOMP].data[si].T[wi];
								Z = S.CI[ZCOMP].data[si].T[wi];
							}
							else if (S.InvertTotalField) {
								X = S.CI[XCOMP].data[si].T[wi];
								Z = S.CI[ZCOMP].data[si].T[wi];
							}
							else {
								X = S.CI[XCOMP].data[si].S[wi];
								Z = S.CI[ZCOMP].data[si].S[wi];
							}

							RT ov = AEM::ewise_hypot(X, Z);
							set_data_vector(Observed_unculled, vi, ov);

							const RT& Xerr = S.CI[XCOMP].data[si].E[wi];
							const RT& Zerr = S.CI[ZCOMP].data[si].E[wi];
							RT ev = combine_noises(X,Xerr,Z,Zerr);
							set_data_vector(Error_unculled, vi, ev);
						}
					}

					for (size_t ci = 0; ci < NCOMP; ci++) {
						if (unused_or_composite_component(S, ci)) continue;
						if (S.CI[ci].Use) {
							const SoundingData<RT>& d = S.CI[ci].data[si];
							for (size_t wi = 0; wi < S.nWindows(); wi++) {
								int vi = vindex(si, sysi, ci, wi);
								RT ov;
								if (S.InvertPSI) ov = d.T[wi];
								else if (S.InvertTotalField) ov = d.T[wi];
								else ov = d.S[wi];
								set_data_vector(Observed_unculled, vi, ov);

								const RT& ev = d.E[wi];
								set_data_vector(Error_unculled, vi, ev);
							}
						}
					}
				}
			}

			if (ErrorAddition > 0.0) {
				for (size_t k = 0; k < Error_unculled.size(); k++) {
					Error_unculled[k] = Error_unculled[k] + ErrorAddition;
				}
			}

			//Work out indices to be culled
			ActiveData.clear();
			for (size_t i = 0; i < nAllData; i++) {
				if (!isnull(Observed_unculled[i]) && !isnull(Error_unculled[i])) ActiveData.push_back(i);
			}
			nData = ActiveData.size();

			if (nData != nAllData) {
				size_t ncull = nAllData - nData;
				OutputMessage += strprint(", %d null data/noise were culled", (int)ncull);
			}
			Err = cull(Error_unculled);
			Obs = cull(Observed_unculled);

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
		};

		void initialise_bunch_parameters() {

			const size_t nsys = nSystems();
			for (size_t si = 0; si < nSoundings; si++) {
				const EarthStore& e = EStore[si];
				const GeometryStore& g = GStore[si];
				if (solve_conductivity()) {
					for (size_t li = 0; li < nLayers; li++) {
						RefParam[cindex(si, li)] = log10(e.refval.conductivity[li]);
						RefParamStd[cindex(si, li)] = e.refvalstd.conductivity[li];
					}
				}

				if (solve_thickness()) {
					for (size_t li = 0; li < nLayers - 1; li++) {
						RefParam[tindex(si, li)] = log10(e.refval.thickness[li]);
						RefParamStd[tindex(si, li)] = e.refvalstd.thickness[li];
					}
				}

				for (int gi = 0; gi < TDEmGeometry::NELEM; gi++) {
					const std::string gname = TDEmGeometry::element_name(gi);
					const int pi = gindex(si, gi);
					if (pi >= 0) {
						RefParam[pi] = g.refval[gname];
						RefParamStd[pi] = g.refvalstd[gname];
					}
				}

				//ggaoffsets
				if (true) {
					for (size_t sysi = 0; sysi < nsys; sysi++) {
						for (size_t ci = 0; ci < NCOMP; ci++) {
							const CompInfo& C = compinfo(sysi, ci);
							if (C.Use) {
								const int pi = ggaoffsetindex(si, sysi, ci);
								if (pi >= 0) {
									RefParam[pi] = C.ggaoffsetStore.refval[si];
									RefParamStd[pi] = C.ggaoffsetStore.refvalstd[si];
								}
							}
						}
					}
				}
			}

			//Scaling params
			if (solve_scalingfactors()) {
				for (size_t sysi = 0; sysi < SI.size(); sysi++) {
					for (size_t ci = 0; ci < NCOMP; ci++) {
						const CompInfo& C = compinfo(sysi, ci);
						if (C.Use) {
							const int pi = sfindex(sysi, ci);
							if (pi >= 0) {
								RefParam[pi] = C.sfStore.refval;
								RefParamStd[pi] = C.sfStore.refvalstd;
							}
						}
					}
				}
			}
		};

		double max_step_fraction_to_bound(const Vector& m_old, const Vector& m_new, const Vector& dm) {
			static const double ud = undefinedvalue<double>();
			double maxf = 1.0;
			for (size_t pi = 0; pi < nParam; pi++) {
				if (Param_Min[pi] == ud) continue;

				if (m_new[pi] <= Param_Min[pi]) {
					double f = (Param_Min[pi] - m_old[pi]) / dm[pi];
					maxf = std::min(f, maxf);
				}
				else if (m_new[pi] >= Param_Max[pi]) {
					double f = (Param_Max[pi] - m_old[pi]) / dm[pi];
					maxf = std::min(f, maxf);
				}
			}
			return maxf;
		}

		Vector parameter_change(const double& lambda, const Vector& m_old, const Vector& pred) {
			const Vector m_new = solve_linear_system(lambda, m_old, pred);
			//std::cout << m_old;
			Vector dm = m_new - m_old;
			//std::cout << dm;
			double maxf = max_step_fraction_to_bound(m_old, m_new, dm);
			if (maxf < 1.0) return dm *= maxf;
			return dm;
		};

		std::vector<Earth1D> get_earth(const Vector& parameters)
		{
			std::vector<Earth1D> ev(nSoundings);;
			for (size_t si = 0; si < nSoundings; si++) {
				ev[si] = EStore[si].refval;
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

		std::vector<TDEmGeometry> get_geometry(const Vector& parameters) const {
			std::vector<TDEmGeometry> gv(nSoundings);
			for (size_t si = 0; si < nSoundings; si++) {
				gv[si] = GStore[si].input;
				for (int gi = 0; gi < TDEmGeometry::NELEM; gi++) {
					const std::string gname = TDEmGeometry::element_name(gi);
					const int pi = gindex(si, gname);
					if (pi >= 0) {
						gv[si][gname] = parameters[pi];
					}
				}
			}
			return gv;
		}

		Vec3d get_gga_offsets(const size_t& sysi, const size_t& si, const Vector& parameters) const {
			Vec3d v(0.0, 0.0, 0.0);
			const SysInfo& S = SI[sysi];
			for (size_t ci = 0; ci < NCOMP; ci++) {
				const CompInfo& C = S.CI[ci];
				if (C.Use) {
					v[ci] = C.ggaoffsetStore.input[si];
					const int pi = ggaoffsetindex(si, sysi, ci);
					if (pi >= 0) v[ci] = parameters[pi];
				}
			}
			//std::cout << v;
			//std::cout << parameters;
			return v;
		}

		Vec3d get_scalefactors(const size_t sysi, const Vector& parameters) const {
			Vec3d sf(1.0, 1.0, 1.0);
			const SysInfo& S = SI[sysi];
			for (int ci = 0; ci < NCOMP; ci++) {				
				const CompInfo& C = S.CI[ci];
				if (C.Use) {
					//sf[ci] = C.sfStore.input;
					const int pi = sfindex(sysi, ci);
					if (pi >= 0) {
						sf[ci] = parameters[pi];
					}
				}
			}
			return sf;
		};

		// Forward models and derivatives
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

		void set_data_vector(Vector& vec, const int& value_index, const double& value) const {
			vec[value_index] = value;
		};

		void set_data_vector(Vector& vec, const int& value_index, const cdouble& value) const {
			vec[2 * value_index] = value.real();
			vec[1 + 2 * value_index] = value.imag();
		};

		void set_mat(Matrix& mat, const int& value_index, const size_t& parameter_index, const double& value) {
			mat(value_index, parameter_index) = value;
		};

		void set_mat(Matrix& mat, const int& value_index, const size_t& parameter_index, const cdouble& value) {
			mat(2 * value_index, parameter_index) = value.real();
			mat(1 + 2 * value_index, parameter_index) = value.imag();
		};

		void forwardmodel_impl(const Vector& parameters, Vector& predicted, Matrix& jacobian, bool computederivatives) {
			
			//D = dgga + (P + S)/ga 
			
			Vector pred_all(nAllData);
			Matrix J_all;
			if (computederivatives) {
				J_all.resize(nAllData, nParam);
				J_all.setZero();
			}

			//bookmark
			const size_t nsys = nSystems();
			std::vector<Earth1D> ev = get_earth(parameters);
			std::vector<TDEmGeometry> gv = get_geometry(parameters);
			for (size_t sysi = 0; sysi < nsys; sysi++) {
				SysInfo& S = SI[sysi];
				AEMSystem<RT>& AEMSystem = S.sys();
				const size_t nw = AEMSystem.nWindows();

				//Vec3d scalefactors(1.0, 1.0, 1.0);
				bool solvesf = solve_scalingfactors();
				//if (solvesf)
				const Vec3d scalefactors = get_scalefactors(sysi, parameters);
				//std::cout << scalefactors << std::endl;
								
				for (size_t si = 0; si < nSoundings; si++) {
					const Earth1D& E = ev[si];
					const TDEmGeometry& G = gv[si];

					TDEmResponse<RT> R = AEMSystem.forward_model(E, G);

					TDEmVectorResponse<RT> FM;
					if (S.InvertTotalField) FM = R.totalfield();
					else FM = R.S;
					
					//if (solvesf)
					FM.scale_components(scalefactors);

					Vec3d ga_scaled;
					if (S.InvertPSI) {
						ga_scaled = S.get_scaled_ga(si);
						Vec3d ggaoffsets = get_gga_offsets(sysi, si, parameters);
						FM.divide_components(ga_scaled);
						FM.plus_components(ggaoffsets);
					};

					// Predicted
					if (S.InvertXYZAmplitude) {
						TDEmScalarResponse<RT> AMPFM = FM.xyzamp();
						for (size_t wi = 0; wi < nw; wi++) {
							int di = vindex(si, sysi, xyzcomp, wi);
							set_data_vector(pred_all, di, AMPFM[wi]);
						}
					}
					else if (S.InvertXZAmplitude) {
						TDEmScalarResponse<RT> AMPFM = FM.xzamp();
						for (size_t wi = 0; wi < nw; wi++) {
							int vi = vindex(si, sysi, xzcomp, wi);
							set_data_vector(pred_all, vi, AMPFM[wi]);
						}
					}

					for (size_t ci = 0; ci < NCOMP; ci++) {
						if (unused_or_composite_component(S, ci)) continue;
						for (size_t wi = 0; wi < nw; wi++) {
							int vi = vindex(si, sysi, ci, wi);
							set_data_vector(pred_all, vi, FM[ci][wi]);
						}
					}

					// Jacobian
					if (computederivatives) {
						TDEmVectorResponse<RT> DRV(nw);
						
						// Scale factor derivatives
						if (solvesf) {
							for (size_t ci = 0; ci < NCOMP; ci++) {
								if (S.CI[ci].Use) {
									const int pindex = sfindex(sysi, ci);
									if (pindex >= 0) {
										DRV = FM;
										// Here filling DRV with the forward model itself as the derivative w.r.t to scale factor param is the forward model itself
										// But zero for other components
										Vec3d f(0, 0, 0);
										f[ci] = 1.0;
										DRV.scale_components(f);
										if (S.InvertPSI) DRV.divide_components(ga_scaled);
										fillMatrixColumn(J_all, pindex, si, sysi, FM, DRV);
									}
								}
							}
						}

						// GGA Offset derivatives
						if (solve_ggaoffsets()) {
							for (size_t ci = 0; ci < NCOMP; ci++) {
								if (S.CI[ci].Use) {
									const int pindex = ggaoffsetindex(si,sysi,ci);
									if (pindex >= 0) {
										// Here filling with 1's for real of component ci
										//      and 0's for other components and for imahinary
										Vec3d v(0, 0, 0);
										v[ci] = 1.0;
										DRV.set_values(v);
										//std::cout << DRV << std::endl;
										fillMatrixColumn(J_all, pindex, si, sysi, FM, DRV);
									}
								}
							}
						}

						// Conductivity derivatives
						if (solve_conductivity()) {
							for (size_t li = 0; li < nLayers; li++) {
								const int pindex = cindex(si, li);
								R = AEMSystem.derivative(G,CalculationType(CMode::DC, li));
								if (S.InvertTotalField) DRV = R.totalfield();
								else DRV = R.S;
								//multiply by natural log(10) as parameters are in logbase10 units
								const double f = log(10.0) * E.conductivity[li];
								DRV *= f;
								//if(solvesf) 
								DRV.scale_components(scalefactors);
								if (S.InvertPSI) DRV.divide_components(ga_scaled);
								fillMatrixColumn(J_all, pindex, si, sysi, FM, DRV);
							}
						}

						// Thickness derivatives
						if (solve_thickness()) {
							for (size_t li = 0; li < nLayers - 1; li++) {
								const int pindex = tindex(si, li);
								R = AEMSystem.derivative(G,CalculationType(CMode::DT, li));
								if (S.InvertTotalField) DRV = R.totalfield();
								else DRV = R.S;
								//multiply by natural log(10) as parameters are in logbase10 units
								double f = log(10.0) * E.thickness[li];
								DRV *= f;
								//if (solvesf) 
								DRV.scale_components(scalefactors);
								if (S.InvertPSI) DRV.divide_components(ga_scaled);
								fillMatrixColumn(J_all, pindex, si, sysi, FM, DRV);
							}
						}

						// Geometry derivatives
						if (FreeGeometry) {
							const size_t ng = G.nelem();
							for (size_t gi = 0; gi < ng; gi++) {
								const std::string gname = G.element_name(gi);
								if (solve_geometry_elementname(gname)) {
									const CMode cmode = G.derivative_mode(gi);
									const CalculationType ctype(cmode);
									const size_t pindex = gindex(si, gname);
									R = AEMSystem.derivative(G, ctype);
									if (S.InvertTotalField) DRV = R.totalfield();
									else DRV = R.S;
									//if (solvesf)
									DRV.scale_components(scalefactors);
									if (S.InvertPSI) DRV.divide_components(ga_scaled);
									//std::cout << DRV << std::endl;
									fillMatrixColumn(J_all, pindex, si, sysi, FM, DRV);
								}
							}
						}
					}//Derivatives block
				}//sounding loop
			}//system loop
			predicted = cull(pred_all);
			if (computederivatives) jacobian = cull(J_all);

			if (Verbose && computederivatives) {
				//std::cerr << "\n-----------------\n";
				//std::cerr << "J_all: It " << CIS.iteration + 1 << std::endl;			
				//std::cerr << J_all;			
				//std::cerr << "\n-----------------\n";
			}

			if (OutputOpt.Dump && computederivatives) {
				std::string dp = dumppath().string() + "J.dat";
				writetofile(J_all, dp);
			}
		}

		void set_predicted(const Vector& parameters) {
			std::vector<Earth1D> ev = get_earth(parameters);
			std::vector<TDEmGeometry> gv = get_geometry(parameters);
			const size_t nsys = nSystems();
			for (size_t sysi = 0; sysi < nsys; sysi++) {
				SysInfo& S = SI[sysi];
				S.predicted.resize(nSoundings);

				AEMSystem<RT>& A = S.sys();
				const size_t& nw = A.nWindows();
				for (size_t si = 0; si < nSoundings; si++) {
					const Earth1D& e = ev[si];
					const TDEmGeometry& g = gv[si];
					S.predicted[si] = A.forward_model(e, g);
				}
			}
		};

		double amplitude_derivative(const double& x, const double& dxdp, const double& y, const double& dydp) const {
			// f(x,y) = (x*x + y*y)^0.5
			// df/dx = 1/2*2x*(x*x + y*y)^-0.5=x/f
			// df/dy = 1/2*2y*(x*x + y*y)^-0.5=y/f
			// df/dp = df/dx * dx/dp + df/dy * dy/dp
			//       = x/f * dx/dp + y/f * dy/dp
			const double f = std::hypot(x, y);
			const double d = (x * dxdp + y * dydp) / f;
			return d;
		};

		double amplitude_derivative(const double& x, const double& dxdp, const double& y, const double& dydp, const double& z, const double& dzdp) const {
			const double f = std::hypot(x, y, z);
			const double d = (x * dxdp + y * dydp + z * dzdp) / f;
			return d;
		};

		cdouble amplitude_derivative(const cdouble& x, const cdouble& dxdp, const cdouble& y, const cdouble& dydp) const {
			const double r = amplitude_derivative(x.real(), dxdp.real(), y.real(), dydp.real());
			const double i = amplitude_derivative(x.imag(), dxdp.imag(), y.imag(), dydp.imag());
			return cdouble(r, i);
		};

		cdouble amplitude_derivative(const cdouble& x, const cdouble& dxdp, const cdouble& y, const cdouble& dydp, const cdouble& z, const cdouble& dzdp) const {
			const double r = amplitude_derivative(x.real(), dxdp.real(), y.real(), dydp.real(), z.real(), dzdp.real());
			const double i = amplitude_derivative(x.imag(), dxdp.imag(), y.imag(), dydp.imag(), z.imag(), dzdp.imag());
			return cdouble(r, i);
		};

		TDEmScalarResponse<RT> xz_amplitude_derivative(const TDEmVectorResponse<RT>& FM, const TDEmVectorResponse<RT>& DRV) const {
			const size_t nw = FM.nWindows();
			TDEmScalarResponse<RT> SD(nw);
			for (size_t wi = 0; wi < nw; wi++) {
				const RT& x    = FM[XCOMP][wi];
				const RT& dxdp = DRV[XCOMP][wi];
				const RT& z    = FM[ZCOMP][wi];
				const RT& dzdp = DRV[ZCOMP][wi];
				SD[wi] = amplitude_derivative(x, dxdp, z, dzdp);
			}
			return SD;
		};

		TDEmScalarResponse<RT> xyz_amplitude_derivative(const TDEmVectorResponse<RT>& FM, const TDEmVectorResponse<RT>& DRV) const {
			const size_t nw = FM.nWindows();
			TDEmScalarResponse<RT> SD(nw);
			for (size_t wi = 0; wi < nw; wi++) {
				const RT& x    = FM[XCOMP][wi];
				const RT& dxdp = DRV[XCOMP][wi];
				const RT& y    = FM[YCOMP][wi];
				const RT& dydp = DRV[YCOMP][wi];
				const RT& z    = FM[ZCOMP][wi];
				const RT& dzdp = DRV[ZCOMP][wi];
				SD[wi] = amplitude_derivative(x, dxdp, y, dydp, z, dzdp);
			}
			return SD;
		};

		void fillMatrixColumnComponent(Matrix& M, const size_t& pindex, const size_t& si, const size_t& sysi, const size_t& ci, const TDEmScalarResponse<RT>& SD) {
			const size_t nw = SD.size();
			for (size_t wi = 0; wi < nw; wi++) {
				const int vi = vindex(si, sysi, ci, wi);
				set_mat(M, vi, pindex, SD[wi]);
			}
		};

		void fillMatrixColumn(Matrix& M, const size_t& pindex, const size_t& si, const size_t& sysi, const TDEmVectorResponse<RT>& FM, const TDEmVectorResponse<RT>& DRV) {
			const SysInfo& S = SI[sysi];
			TDEmScalarResponse<RT> SD;
			if (S.InvertXYZAmplitude) {
				SD = xyz_amplitude_derivative(FM, DRV);
				//std::cout << SD << std::endl;
				fillMatrixColumnComponent(M, pindex, si, sysi, xyzcomp, SD);
			}
			else if (S.InvertXZAmplitude) {
				SD = xz_amplitude_derivative(FM, DRV);
				fillMatrixColumnComponent(M, pindex, si, sysi, xzcomp, SD);
			}

			for (size_t ci = 0; ci < NCOMP; ci++) {
				if (unused_or_composite_component(S, ci)) continue;
				SD = DRV.component(ci);
				//std::cout << SD << std::endl;
				fillMatrixColumnComponent(M, pindex, si, sysi, ci, SD);
			}
		};

		// Etc
		void save_iteration_file(const cIterationState& S) const {
			std::ofstream ofs(dumppath().string() + "iteration.dat");
			ofs << S.info_string();
		};

		bool read_bunch(const size_t& record) {
			bool bunchstatus = false;
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

			//if (solve_scalingfactors()) {
				for (size_t sysi = 0; sysi < SI.size(); sysi++) {
					for (size_t ci = 0; ci < NCOMP; ci++) {
						CompInfo& C = SI[sysi].CI[ci];
						if (C.Use) {
							const Parameter& p = C.scalingfactor;
							if (p.initialised()) {
								ScaleFactorsStore& s = C.sfStore;
								IM->read(p.get_fd(Parameter::INPUT), s.input, 1);
								//if (p.solve()) {
									IM->read(p.get_fd(Parameter::REF), s.refval, 1);
									IM->read(p.get_fd(Parameter::STD), s.refvalstd, 1);
									if (p.bound()) {
										IM->read(p.get_fd(Parameter::MIN), s.minval, 1);
										IM->read(p.get_fd(Parameter::MAX), s.maxval, 1);
									}
								//}
							}
						}
					}
				}
			//}

			return true;
		}

		bool read_record(const size_t& bunchsoundingindex) {
			const size_t& si = bunchsoundingindex;
			bool readstatus = true;
			EarthStore& e = EStore[si];
			GeometryStore& g = GStore[si];
			
			if (IM->parse_record() == false) return false;

			bool status;
			Id[si].uniqueid = (int)IM->record();

			status = read_ancillary_fields(si);
			status = read_geometry(si, ifdGMap);

			status = IM->read(pConductivity.get_fd(Parameter::INPUT), e.refval.conductivity, nLayers); if (status == false) readstatus = false;
			if (pConductivity.solve()) {
				status = IM->read(pConductivity.get_fd(Parameter::REF), e.refval.conductivity, nLayers); if (status == false) readstatus = false;
				status = IM->read(pConductivity.get_fd(Parameter::STD), e.refvalstd.conductivity, nLayers); if (status == false) readstatus = false;
				if (pConductivity.bound()) {
					status = IM->read(pConductivity.get_fd(Parameter::MIN), e.minval.conductivity, nLayers); if (status == false) readstatus = false;
					status = IM->read(pConductivity.get_fd(Parameter::MAX), e.maxval.conductivity, nLayers); if (status == false) readstatus = false;
				}
			}

			status = IM->read(pThickness.get_fd(Parameter::INPUT), e.refval.thickness, nLayers-1); if (status == false) readstatus = false;
			if (pThickness.solve()) {
				status = IM->read(pThickness.get_fd(Parameter::REF), e.refval.thickness, nLayers-1); if (status == false) readstatus = false;
				status = IM->read(pThickness.get_fd(Parameter::STD), e.refvalstd.thickness, nLayers-1); if (status == false) readstatus = false;
				if (pThickness.bound()) {
					status = IM->read(pThickness.get_fd(Parameter::MIN), e.minval.thickness, nLayers-1); if (status == false) readstatus = false;
					status = IM->read(pThickness.get_fd(Parameter::MAX), e.maxval.thickness, nLayers-1); if (status == false) readstatus = false;
				}
			}
			else {
				//This is purely only so the the refvalstd 
				e.refvalstd.thickness = std::vector<double>(nLayers - 1,0.0);
			}
			e.sanity_check();

			const size_t nsys = nSystems();
			for (size_t sysi = 0; sysi < nsys; sysi++) {
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
					glog.errormsg(_SRC_, "Bad variant type\n");
				}

				return true;
			}
			return false;
		}

		bool read_geometry(const size_t& bunchindex, IFDKeyMap& map) {
			bool status = true;
			const size_t si = bunchindex;
			GeometryStore& gstore = GStore[si];
			for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
				std::string ename = TDEmGeometry::element_name(gi);
				const Parameter& p = map.at(ename);
				bool inpstatus = IM->read(p.get_fd(Parameter::INPUT), gstore.input[gi]);
				bool refstatus = IM->read(p.get_fd(Parameter::REF), gstore.refval[gi]);

				if (refstatus == false && inpstatus == true) {
					gstore.refval[gi] = gstore.input[gi];
					refstatus = true;
				}
				else if (inpstatus == false && refstatus == true) {
					gstore.input[gi] = gstore.refval[gi];
					inpstatus = true;
				}

				if (inpstatus == false) {
					std::ostringstream msg;
					msg << "No 'Input or Ref' defined for " << ename << std::endl;
					glog.errormsg(_SRC_, msg.str());
				}

				if (refstatus == false) {
					std::ostringstream msg;
					msg << "No 'Ref or Input' defined for " << ename << std::endl;
					glog.errormsg(_SRC_, msg.str());
				}

				if (p.solve()) {
					bool stdstatus = IM->read(p.get_fd(Parameter::STD), gstore.refvalstd[gi]);
					if (stdstatus == false) {
						std::ostringstream msg;
						msg << "No 'Std' defined for " << ename << std::endl;
						glog.errormsg(_SRC_, msg.str());
					}

					if (p.bound()) {
						bool minstatus = IM->read(p.get_fd(Parameter::MIN), gstore.minval[gi]);
						bool maxstatus = IM->read(p.get_fd(Parameter::MAX), gstore.maxval[gi]);
					}
				}

				if (reconstruct_primary()) {
					const cFieldDefinition& fd = fdPFRGvec[gi];
					if (fd.isinitialised()) {
						bool pfrstatus = IM->read(fd, gstore.pfr[gi]);
						if (pfrstatus == false) {
							std::ostringstream msg;
							msg << "Could not read the specified 'PFR' field for geometry element " << ename << std::endl;
							glog.errormsg(_SRC_, msg.str());
						}
					}
					else {
						//PFR geometry is not defined so make it equal to input geometry
						gstore.pfr[gi] = gstore.input[gi];
					}
				}
			}
			return status;
		};

		void read_system_data(size_t& sysindex, const size_t& soundingindex) {
			SysInfo& S = SI[sysindex];
			for (size_t ci = 0; ci < NCOMP; ci++) {
				CompInfo& C = S.CI[ci];
				if (C.Use) {
					C.readdata(IM, soundingindex);
					C.read_ggaoffset_data(IM, soundingindex);
					if (C.EstimateNoiseFromModel) {
						C.estimate_noise_from_model(soundingindex, S.InvertPSI, S.InvertTotalField);
					}
				}
			}
		};

		void dump_first_iteration() {

			const std::string dp = dumppath().string();
			makedirectory(dumppath());

			const size_t si = Bunch.master_index();
			GeometryStore& g = GStore[si];
			EarthStore& e = EStore[si];
			SampleId& id = Id[si];

			dump_id_info(id, dp + "Id_Info.dat");
			dump_system_info(dp + "System_Info.dat");
			dump_component_info(dp + "Component_Info.dat");
			if (SI[0].InvertPSI) {
				dump_gga_info(dp + "gga_Info.dat");
			}

			write(Obs, dp + "observed.dat");
			write(Err, dp + "observed_std.dat");

			g.refval.write(dp + "geometry_start.dat");
			e.refval.write(dp + "earth_start.dat");

			g.refval.write(dp + "geometry_ref.dat");
			e.refval.write(dp + "earth_ref.dat");

			g.refvalstd.write(dp + "geometry_std.dat");
			e.refvalstd.write(dp + "earth_std.dat");

		};

		void dump_id_info(const SampleId& id, const fs::path& path) {
			std::ofstream of(path);
			of << "UniquwId " << id.uniqueid << std::endl;
			of << "Survey " << id.survey << std::endl;
			of << "Date " << id.date << std::endl;
			of << "Flight " << id.flight << std::endl;
			of << "Line " << id.line << std::endl;
			of << "Fiducial " << id.fiducial << std::endl;
			of << "X " << id.x << std::endl;
			of << "Y " << id.y << std::endl;
			of << "Elevation " << id.elevation << std::endl; \
		};

		void dump_system_info(const fs::path& path) {
			std::ofstream of(path);
			std::string comma = ",";
			of << "Index" << comma
				<< "Name" << comma
				<< "SystemTyps" << comma
				<< "Units" << comma
				<< "InvertTotalField" << comma
				<< "InvertPSI" << comma
				<< "InvertXZAmplitude" << comma
				<< "InvertXYZAmplitude" << comma
				<< "ReconstructPrimary" << std::endl;

			for (size_t sysi = 0; sysi < nSystems(); sysi++) {
				const SysInfo& S = SI[sysi];
				of	<< sysi << comma
					<< S.sys().name() << comma
					<< S.sys().type_string() << comma
					<< S.Units << comma
					<< S.InvertTotalField << comma
					<< S.InvertPSI << comma
					<< S.InvertXZAmplitude << comma
					<< S.InvertXYZAmplitude << comma
					<< S.ReconstructPrimary << std::endl;
			}
		};

		void dump_component_info(const fs::path& path) {
			std::ofstream of(path);
			of	<< "SystemIndex" << ","
				<< "ComponentIndex" << ","
				<< "Name" << ","
				<< "Use" << ","
				<< "nWindows" << ","
				<< "nElements" << ","
				<< "nChannels" << std::endl;

			for (size_t sysi = 0; sysi < nSystems(); sysi++) {
				const SysInfo& S = SI[sysi];
				if (S.InvertXYZAmplitude) {
					const CompInfo& C = S.CI[0];
					of << sysi << ","
						<< 0 << ","
						<< "XYZAMP" << ","
						<< C.Use << ","
						<< C.nWindows() << ","
						<< C.nElements() << ","
						<< C.nChannels() << std::endl;
				}
				else if (S.InvertXZAmplitude) {
					const CompInfo& C = S.CI[0];
					of << sysi << ","
						<< 0 << ","
						<< "XZAMP" << ","
						<< C.Use << ","
						<< C.nWindows() << ","
						<< C.nElements() << ","
						<< C.nChannels() << std::endl;
				}
				else {
					for (size_t ci = 0; ci < NCOMP; ci++) {
						const CompInfo& C = S.CI[ci];
						of << sysi << ","
							<< ci << ","
							<< C.Name << ","
							<< C.Use << ","
							<< C.nWindows() << ","
							<< C.nElements() << ","
							<< C.nChannels() << std::endl;
					}
				}
			}
		};

		void dump_gga_info(const fs::path& path) {
			std::ofstream of(path);
			of  << "SystemIndex" << ","
				<< "ComponentIndex" << ","
				<< "SampleIndex" << ","
				<< "gga" << std::endl;

			for (size_t sysi = 0; sysi < nSystems(); sysi++) {
				const SysInfo& S = SI[sysi];
				for (size_t ci = 0; ci < NCOMP; ci++) {
					const CompInfo& C = S.CI[ci];
					for (size_t si = 0; si < nSoundings; si++) {
						double gga = 0;
						if (C.Use) gga = C.get_gga(si);
						of << sysi << ","
							<< ci << ","
							<< si << ","
							<< gga << std::endl;
					}
				}
			}
		};

		void dump_iteration(const cIterationState& state) {
			const std::string dp = dumppath().string();
			makedirectory(dp);
			writetofile(Obs, dp + "d.dat");
			writetofile(Err, dp + "e.dat");
			writetofile(state.param, dp + "m.dat");
			writetofile(state.pred, dp + "g.dat");
			std::vector<Earth1D> e = get_earth(state.param);
			std::vector <TDEmGeometry> g = get_geometry(state.param);
			e[Bunch.master_index()].write(dp + "earth_inv.dat");
			g[Bunch.master_index()].write(dp + "geometry_inv.dat");
			dump_earth_all(e, dp + "earth_all.dat");
			dump_geometry_all(g, dp + "geometry_all.dat");
			save_iteration_file(state);
		}

		void dump_earth_all(const std::vector <Earth1D> e, const std::string& path) {
			std::ofstream of(path);
			for (size_t si = 0; si < e.size(); si++) {
				const std::vector<double>& c = e[si].conductivity;
				const std::vector<double>t = e[si].dummy_thickness();
				for (size_t li = 0; li < e[si].nlayers(); li++) {
					of << t[li] << " " << c[li] << std::endl;
				}
			}
		}

		void dump_geometry_all(const std::vector <TDEmGeometry> g, const std::string& path) {
			std::ofstream of(path);
			for (size_t si = 0; si < g.size(); si++) {
				for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
					of << g[si][gi] << std::endl;
				}
			}
		};

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
			function_trace;
			//std::cout << "iterate"  <<std::endl;
			setup_parameter_bounds();
			CIS.iteration = 0;
			//CIS.lambda = 1e8;
			CIS.param = RefParam;
			forwardmodel(CIS.param, CIS.pred);
			CIS.phid = phiData(CIS.pred);
			CIS.targetphid = CIS.phid;
			CIS.phim = phiModel(CIS.param);
			TerminationReason = "Has not terminated";

			double percentimprovement = 100.0;
			bool   keepiterating = true;
			while (keepiterating == true) {
				if (Verbose && nScalingParam > 0) {
					Vec3d scalefactors = get_scalefactors(0, CIS.param);
					std::cout << "Scaling Factors ";
					for (size_t ci = 0; ci < NCOMP; ci++) {
						std::cout << scalefactors[ci] << " ";
					}
					std::cout << std::endl;
				}

				if (CIS.iteration >= MaxIterations) {
					keepiterating = false;
					TerminationReason = "Too many iterations";
				}
				else if (CIS.iteration > 0 && CIS.phid <= MinimumPhiD) {
					keepiterating = false;
					TerminationReason = "Reached minimum";
				}
				else if (percentimprovement < 0) {
					keepiterating = false;
					TerminationReason = "No improvement";
				}
				else if (CIS.iteration > 10 && percentimprovement < MinimumImprovement) {
					keepiterating = false;
					TerminationReason = "Small % improvement";
				}
				else {
					if (Verbose) std::cerr << CIS.info_string();
					if (CIS.iteration + 1 >= BeginGeometrySolveIteration) FreeGeometry = true;
					else FreeGeometry = false;

					const double targetphid = std::max(CIS.phid * 0.7, MinimumPhiD);
					
					Vector g;
					forwardmodel_and_jacobian(CIS.param, g, J);
					if (CIS.iteration == 0) {
						CIS.lambda = 1e8;
						if (OutputOpt.Dump) {
							dump_first_iteration();
							dump_iteration(CIS);
						}
					}

					const cTrial t = lambda_search_targetphid(CIS.lambda, targetphid);
					const Vector dm = parameter_change(t.lambda, CIS.param, CIS.pred);
					const Vector m = CIS.param + (t.stepfactor * dm);

					forwardmodel(m, g);
					const double phid = phiData(g);
					percentimprovement = 100.0 * (CIS.phid - phid) / (CIS.phid);

					if (phid <= CIS.phid) {
						CIS.iteration++;
						CIS.param = m;
						CIS.pred = g;
						CIS.targetphid = targetphid;
						CIS.phid = phid;
						CIS.lambda = t.lambda;
						CIS.phim = phiModel(CIS.param);
						if (OutputOpt.Dump) dump_iteration(CIS);
					}
				}
			}

			std::vector<Earth1D> ev = get_earth(CIS.param);
			std::vector<TDEmGeometry> gv = get_geometry(CIS.param);
			for (size_t si = 0; si < nSoundings; si++) {
				EStore[si].invmodel = ev[si];
				GStore[si].invmodel = gv[si];
			}

			set_predicted(CIS.param);
			forwardmodel_and_jacobian(CIS.param, CIS.pred, J);
			ParameterSensitivity = compute_parameter_sensitivity();
			ParameterUncertainty = compute_parameter_uncertainty();
		}

		int execute() {
			function_trace;
			bool readstatus = true;
			int paralleljob = 0;
			do {
				int record = ((int)StartRecord - 1) + paralleljob * (int)IM->subsamplerate();
				if (record > (EndRecord - 1))break;
				if ((paralleljob % Size) == Rank) {
					std::ostringstream s;
					//std::cout << "read_bunch" << std::endl;
					if ((readstatus = read_bunch(record))) {
						s << bunch_id();
						//std::cout << "initialise_bunch" << std::endl;
						if (initialise_bunch()) {
							double t1 = gettime();
							//std::cout << "iterate" << std::endl;
							iterate();
							double t2 = gettime();
							double etime = t2 - t1;
							//std::cout << "write_result" << std::endl;
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
						glog.logmsg(0, s.str());
					}
				}
				paralleljob++;
			} while (readstatus == true);
			glog.close();
			return 0;
		}

		double phiModel(const Vector& m) const {
			double v = 0.0;
			v += LCrefc.phi(m, RefParam);
			v += LCreft.phi(m, RefParam);
			v += LCrefg.phi(m, RefParam);
			v += LCrefo.phi(m, RefParam);
			v += LCrefs.phi(m, RefParam);
			v += LCvcsmth.phi(m, RefParam);
			v += LCvcsim.phi(m, RefParam);
			v += LClatc.phi(m, RefParam);
			v += LClatg.phi(m, RefParam);

			Vector clfwd = CableLengthConstraint_forward(m);
			v += NLCcablen.phi(clfwd);

			Vector bndfwd = BoundsConstraint_forward(m);
			v += NLCbounds.phi(bndfwd);

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
			double elambda = 1000.0 * (s0[0] / s1[0]);
			return elambda;
		}

		Vector solve_linear_system(const double& lambda, const Vector& param, const Vector& pred) {
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

			
			//std::cout << m0;
			//std::cout << RefParamStd;

			//std::cout << Wr;


			Matrix V = Wd;
			if (NormType == NormType::L1) {
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
			//std::cout << b;

			if (LClatg.operates_on_difference_from_reference_model()) {
				b += lambda * (LClatg.W * m0);
			}

			if (NLCcablen.alpha > 0) {
				cNonLinearConstraint& C = NLCcablen;
				CableLengthConstraint_jacobian(m);
				Vector predicted = CableLengthConstraint_forward(m);
				A += C.J.transpose() * C.W.transpose() * C.J;
				b += C.J.transpose() * C.W.transpose() * (C.data - predicted + C.J * m);
			}

			if (NLCbounds.alpha > 0) {
				cNonLinearConstraint& C = NLCbounds;
				BoundsConstraint_jacobian(m);
				Vector predicted = BoundsConstraint_forward(m);
				A += C.J.transpose() * C.W.transpose() * C.J;
				b += C.J.transpose() * C.W.transpose() * (C.data - predicted + C.J * m);
			}
			
			const Eigen::LLT<Matrix> lltOfA(A);
			if (lltOfA.info() == Eigen::NumericalIssue) {
				std::cerr << "\nAt " << bunch_id() << ": The matrix A is possibly non semi - positive definite" << std::endl << A << std::endl;
			}
			Vector x = lltOfA.solve(b);
			return x;
		}

		size_t fill_dataspace_vector(const Vector& vec, const size_t& si, const size_t& sysi, const size_t& ci, std::vector<double>& vout) {
			const size_t nw = SI[sysi].nWindows();
			const size_t nv = _si_.size();
			vout.resize(nw);
			size_t nset = 0;
			for (size_t vi = 0; vi < nv; vi++) {
				if (_si_[vi] == si && _sysi_[vi] == sysi && _ci_[vi] == ci) {
					const size_t& wi = _wi_[vi];
					vout[wi] = vec[vi];
					nset++;
				}
			}
			assert(nset == nw);
			return nset;
		};

		size_t fill_dataspace_vector(const Vector& vec, const size_t& si, const size_t& sysi, const size_t& ci, std::vector<cdouble>& vout) {
			const size_t nw = SI[sysi].nWindows();
			const size_t nv = _si_.size();
			vout.resize(nw);
			size_t nset = 0;
			for (size_t vi = 0; vi < nv; vi++) {
				if (_si_[vi] == si && _sysi_[vi] == sysi && _ci_[vi] == ci) {
					const size_t& wi = _wi_[vi];
					vout[wi] = std::complex<double>(vec[vi*2],vec[1 + vi*2]);
					nset++;
				}
			}
			assert(nset == nw);
			return nset;
		}

		std::vector<RT> get_dataspace_vector(const Vector& vec, const size_t& si, const size_t& sysi, const size_t& ci){
			std::vector<RT> vout;
			fill_dataspace_vector(vec, si, sysi, ci, vout);
			return vout;
		}

		void write_result(const int& pointindex) {
			//bookmark
			const Vector& m = CIS.param;
			const Vector& m0 = RefParam;
			const Vector pred_unculled = un_cull(CIS.pred);


			const int& pi = (int)Bunch.master_record();
			const int& si = (int)Bunch.master_index();
			OM->begin_point_output();

			//Ancillary	
			OM->writefield(pi, Id[si].uniqueid, "uniqueid", "Inversion sequence number", UNITLESS, 1, ST_UINT, DN_NONE, 'I', 12, 0);
			for (size_t fi = 0; fi < AncFld[si].size(); fi++) {
				cFdVrnt& fdv = AncFld[si][fi].second;
				cAsciiColumnField c;
				IM->get_acsiicolumnfield(fdv.fd, c);
				if (fi == 0) c.width++;//First column after unique_id needs to have extra space in case it fills the width to enfore space between columns
				OM->writevrnt(pi, fdv.vnt, c);
			}

			//Geometry Input
			bool invertedfieldsonly = false;
			for (size_t i = 0; i < TDEmGeometry::NELEM; i++) {
				if (invertedfieldsonly && solve_geometry_index(i) == false)continue;
				OM->writefield(pi, GStore[si].input[i], "input_" + GStore[si].input.element_name(i), "Input " + GStore[si].input.description(i), GStore[si].input.units(i), 1, ST_FLOAT, DN_NONE, 'F', 9, 2);
			}

			// Modelled geometry parameters
			const TDEmGeometry& g = GStore[si].invmodel;
			invertedfieldsonly = true;
			for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
				if (invertedfieldsonly && solve_geometry_index(gi) == false)continue;
				OM->writefield(pi, g[gi], "inverted_" + g.element_name(gi), "Inverted " + g.description(gi), g.units(gi), 1, ST_FLOAT, DN_NONE, 'F', 9, 2);
			};

			//Scaling factor parameters
			if (solve_scalingfactors()) {
				const size_t nsys = nSystems();
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					Vec3d sf = get_scalefactors(sysi, m);
					for (size_t ci = 0; ci < NCOMP; ci++) {
						const CompInfo& C = SI[sysi].CI[ci];
						if (C.Use && solve_scalingfactor(sysi, ci)) {
							std::string comp = C.Name;
							std::string fname = "scalingfactor" + strprint("_EMSystem_%d_", (int)sysi + 1) + comp;
							std::string fdesc = "Scaling factor" + strprint(" EMSystem %d ", (int)sysi + 1) + comp + "-component";
							OM->writefield(pi, sf[ci], fname, fdesc, UNITLESS, 1, ST_FLOAT, DN_NONE, 'F', 8, 4);
						}
					}
				}
			}

			//Delta gga	parameters
			if (solve_ggaoffsets()) {
				const size_t nsys = nSystems();
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					Vec3d dgga = get_gga_offsets(sysi, si, m);
					for (size_t ci = 0; ci < NCOMP; ci++) {
						const CompInfo& C = SI[sysi].CI[ci];
						if (C.Use && C.ggaoffset.solve()) {
							std::string comp = C.Name;
							std::string fname = "deltagga" + strprint("_EMSystem_%d_", (int)sysi + 1) + comp;
							std::string fdesc = "Coupling ratio offset" + strprint(" EMSystem %d ", (int)sysi + 1) + comp + "-component";
							OM->writefield(pi, dgga[ci], fname, fdesc, UNITLESS, 1, ST_FLOAT, DN_NONE, 'F', 8, 4);
						}
					}
				}
			}

			//ndata
			OM->writefield(pi, nData, "ndata", "Number of data in inversion", UNITLESS, 1, ST_UINT, DN_NONE, 'I', 4, 0);

			//Earth	parameters
			const Earth1D& e = EStore[si].invmodel;
			OM->writefield(pi, nLayers, "nlayers", "Number of layers ", UNITLESS, 1, ST_UINT, DN_NONE, 'I', 4, 0);
			OM->writefield(pi, e.conductivity, "conductivity", "Layer conductivity", "S/m", e.conductivity.size(), ST_FLOAT, DN_LAYER, 'E', 15, 6);

			if (nLayers > 1) {
				double bottomlayerthickness = 100.0;
				if (solve_thickness() == false && nLayers > 1) {
					bottomlayerthickness = e.thickness[nLayers - 2];
				}
				std::vector<double> thickness = e.thickness;
				thickness.push_back(bottomlayerthickness);
				OM->writefield(pi, thickness, "thickness", "Layer thickness", "m", thickness.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);

				if (OutputOpt.PositiveLayerTopDepths) {
					std::vector<double> dtop = e.layer_top_depth();
					OM->writefield(pi, dtop, "depth_top", "Depth to top of layer", "m", dtop.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);
				}

				if (OutputOpt.NegativeLayerTopDepths) {
					std::vector<double> ndtop = -1.0 * e.layer_top_depth();
					OM->writefield(pi, ndtop, "depth_top_negative", "Negative of depth to top of layer", "m", ndtop.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);
				}

				if (OutputOpt.PositiveLayerBottomDepths) {
					std::vector<double> dbot = e.layer_bottom_depth();
					OM->writefield(pi, dbot, "depth_bottom", "Depth to bottom of layer", "m", dbot.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);
				}

				if (OutputOpt.NegativeLayerBottomDepths) {
					std::vector<double> ndbot = -1.0 * e.layer_bottom_depth();
					OM->writefield(pi, ndbot, "depth_bottom_negative", "Negative of depth to bottom of layer", "m", ndbot.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);
				}

				if (OutputOpt.InterfaceElevations) {
					std::vector<double> etop = e.layer_top_depth();
					etop += Id[si].elevation;
					OM->writefield(pi, etop, "elevation_interface", "Elevation of interface", "m", etop.size(), ST_FLOAT, DN_LAYER, 'F', 9, 2);
				}
			}

			if (OutputOpt.ParameterSensitivity) {
				std::vector<double> ps = copy(ParameterSensitivity);
				if (solve_conductivity()) {
					std::vector<double> v(ps.begin() + cindex(si, 0), ps.begin() + cindex(si, 0) + nLayers);
					OM->writefield(pi, v, "conductivity_sensitivity", "Conductivity parameter sensitivity", UNITLESS, v.size(), ST_FLOAT, DN_LAYER, 'E', 15, 6);
				}

				if (solve_thickness()) {
					std::vector<double> v(ps.begin() + tindex(si, 0), ps.begin() + tindex(si, 0) + nLayers - 1);
					v.push_back(0.0);//halfspace layer not a parameter
					OM->writefield(pi, v, "thickness_sensitivity", "Thickness parameter sensitivity", UNITLESS, v.size(), ST_FLOAT, DN_LAYER, 'E', 15, 6);
				}

				const TDEmGeometry& g = GStore[si].input;
				for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
					if (solve_geometry_index(gi) == true) {
						const std::string gname = g.element_name(gi);
						std::string name = "inverted_" + gname + "_sensitivity";
						std::string desc = g.description(gi) + " parameter sensitivity";
						OM->writefield(pi, ps[gindex(si, gname)], name, desc, UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
					}
				}
			}

			if (OutputOpt.ParameterUncertainty) {
				std::vector<double> pu = copy(ParameterUncertainty);
				if (solve_conductivity()) {
					std::vector<double> v(pu.begin() + cindex(si, 0), pu.begin() + cindex(si, 0) + nLayers);
					OM->writefield(pi, v, "conductivity_uncertainty", "Conductivity parameter uncertainty", "log10(S/m)", v.size(), ST_FLOAT, DN_LAYER, 'E', 15, 6);
				}

				if (solve_thickness()) {
					std::vector<double> v(pu.begin() + tindex(si, 0), pu.begin() + tindex(si, 0) + nLayers - 1);
					v.push_back(0.0);//halfspace layer not a parameter
					OM->writefield(pi, v, "thickness_uncertainty", "Thickness parameter uncertainty", "log10(m)", v.size(), ST_FLOAT, DN_LAYER, 'E', 15, 6);
				}

				const TDEmGeometry& g = GStore[si].input;
				for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
					if (solve_geometry_index(gi) == false) continue;
					const std::string gname = g.element_name(gi);
					std::string name = "inverted_" + gname + "_uncertainty";
					std::string desc = g.description(gi) + " parameter uncertainty";
					OM->writefield(pi, pu[gindex(si, gname)], name, desc, g.units(gi), 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
				}
			}

			//bookmark

			// Output EM data
			write_emdata_version(Observed_unculled, "Observed", "observed data", pi, si);
			write_emdata_version(Error_unculled, "Noise", "estimated data noise", pi, si);
			write_emdata_version(pred_unculled, "Predicted", "predicted data", pi, si);

			//Inversion parameters and norms
			write_result(pi, LCrefc, m, m0);
			write_result(pi, LCreft, m, m0);
			write_result(pi, LCrefg, m, m0);
			write_result(pi, LCrefo, m, m0);
			write_result(pi, LCrefs, m, m0);
			write_result(pi, LCvcsmth, m, m0);
			write_result(pi, LCvcsim, m, m0);
			write_result(pi, LClatc, m, m0);
			write_result(pi, LClatg, m, m0);

			Vector clfwd = CableLengthConstraint_forward(m);
			write_result(pi, NLCcablen, clfwd);

			Vector bndfwd = BoundsConstraint_forward(m);
			write_result(pi, NLCbounds, bndfwd);

			OM->writefield(pi, CIS.phid, "PhiD", "Normalised data misfit", UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
			OM->writefield(pi, CIS.phim, "PhiM", "Combined model norm", UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
			OM->writefield(pi, CIS.lambda, "Lambda", "Lambda regularization parameter", UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
			OM->writefield(pi, CIS.iteration, "Iterations", "Number of iterations", UNITLESS, 1, ST_UINT, DN_NONE, 'I', 4, 0);

			//End of record book keeping
			OM->end_point_output();
			if (nPointsOutput == 0) {
				OM->end_first_record();//only do this once
			}
			nPointsOutput++;
		};

		void write_result(const int& pointindex, const cLinearConstraint& C, const Vector& m, const Vector& m0) {
			if (C.alpha == 0.0) return;

			double phi = 0.0;
			if (C.alpha > 0.0) {
				phi = C.phi(m, m0);
			}

			OM->writefield(pointindex, C.alpha, C.alpha_field_name(), C.alpha_field_description(), UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
			OM->writefield(pointindex, phi, C.phi_field_name(), C.phi_field_description(), UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
		};

		void write_result(const int& pointindex, const cNonLinearConstraint& C, const Vector& predicted) {
			if (C.alpha == 0.0)return;
			double phi = 0.0;
			if (C.alpha > 0.0) {
				phi = C.phi(predicted);
			}
			OM->writefield(pointindex, C.alpha, C.alpha_field_name(), C.description, UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
			OM->writefield(pointindex, phi, C.phi_field_name(), C.phi_field_description(), UNITLESS, 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
		};

		template <typename T>
		void write_plain(const int& pointindex, const cOutputField& of, const std::vector<T>& values) {
			OM->writefield(pointindex, values, of);
		}

		template <typename T>
		void write_plain_or_complex(const int& pointindex, const cOutputField& of, const std::vector<T>& values) {
			write_plain(pointindex, of, values);
		};

		void write_plain_or_complex(const int& pointindex, const cOutputField& of, const std::vector<cdouble>& values) {
			cOutputField ofr = of;
			ofr.name += "_Real";
			std::string& s = ofr.atts.refval(cAsciiColumnField::DESC);
			s += " (real part)";
			std::string s1 = ofr.atts.refval(cAsciiColumnField::DESC);
			write_plain(pointindex, ofr, real(values));

			cOutputField ofi = of;
			ofi.name += "_Imag";
			ofi.atts.refval(cAsciiColumnField::DESC) += " (imaginary part)";
			write_plain(pointindex, ofi, imaginary(values));
		};

		void write_emdata_version(const Vector& vector_unculled, const std::string& version_name, const std::string& version_desc, const size_t& pointindex, const size_t& sampleindex) {
			const std::string& vname = version_name;
			const std::string& vdesc = version_desc;
			const size_t nsys = nSystems();
			cAsciiColumnFormat emfmt('E', 15, 6);

			std::string qname, qdesc;
			std::string cname, cdesc;
			for (size_t sysi = 0; sysi < nsys; sysi++) {
				const SysInfo& S = SI[sysi];
				if (S.InvertPSI) {
					qname = "PSI"; qdesc = "PSI";
				}
				else if (S.InvertTotalField) {
					qname = "TotalField"; qdesc = "total field";
				}
				else {
					qname = "SecondaryField"; qdesc = "secondary field";
				}

				for (size_t ci = 0; ci < NCOMP; ci++) {
					if (unused_or_composite_component(S, ci)) continue;
					const CompInfo& C = S.CI[ci];
					cname = C.name(); cdesc = C.longname();
					std::vector<RT> v = get_dataspace_vector(vector_unculled, sampleindex, sysi, ci);
					writeresult_emdata_array(pointindex, sysi, vname, qname, cname, vdesc, qdesc, cdesc, S.Units, emfmt, v);
				}

				if (S.InvertXYZAmplitude) {
					cname = "XYZAMP"; cdesc = "XYZ-amplitude";
					std::vector<RT> v = get_dataspace_vector(vector_unculled, sampleindex, sysi, xyzcomp);
					writeresult_emdata_array(pointindex, sysi, vname, qname, cname, vdesc, qdesc, cdesc, S.Units, emfmt, v);
				}
				if (S.InvertXZAmplitude) {
					cname = "XZAMP"; cdesc = "XZ-amplitude";
					std::vector<RT> v = get_dataspace_vector(vector_unculled, sampleindex, sysi, xzcomp);
					writeresult_emdata_array(pointindex, sysi, vname, qname, cname, vdesc, qdesc, cdesc, S.Units, emfmt, v);
				}
			}
		};

		//bookmark
		void writeresult_emdata_array(const int& pointindex,
			const size_t& sysindex,
			const std::string& vname, //Observed Predicted Noise
			const std::string& qname, //Primary Secondary Total PSI
			const std::string& cname, //X Y Z XZAMP XYZAMP
			const std::string& vdesc, //Observed Predicted Noise
			const std::string& qdesc, //Primary Secondary Total PSI
			const std::string& cdesc, //X Y Z XZAMP XYZAMP
			const std::string& units,
			const cAsciiColumnFormat& fmt,
			const std::vector<RT>& array)
		{
			const SystemInversionInfo<AEMSystemClass, RT>& S = SI[sysindex];
			const BinaryStorageType btype = ST_FLOAT;
			std::string dimensionname = "em_window";
			std::string sysname = strprint("EMSystem_%d_", (int)sysindex + 1);
			std::string sysdesc = strprint("EMSystem %d ", (int)sysindex + 1);
			const int nbands = array.size();
			
			std::string name = sysname + vname + "_" + qname + "_" + cname;
			std::string desc = sysdesc + vdesc + " " + qdesc + " " + cdesc;

			cOutputField of(name, desc, units, nbands, btype, dimensionname, fmt);
			write_plain_or_complex(pointindex, of, array);
		};
	};
};//End namespace


