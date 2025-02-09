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
		inline static const size_t XZAMP = 3;
		std::vector<std::vector<std::vector<std::vector<int>>>> _vindex_;

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
		size_t nScalingParam = 0;

		const size_t nSystems() const { return SysInfo.size(); };

		size_t nPointsOutput = 0;
		std::vector<GeometryStore> GStore;
		std::vector<EarthStore> EStore;
		cOutputOptions OutputOpt;
		std::vector<AEMSystemInversionInfo<AEMSystemClass, RT>> SysInfo;

		//Column definitions
		cInvertibleFieldDefinition fdC;
		cInvertibleFieldDefinition fdT;
		using IFDKeyMap = std::map<std::string, cInvertibleFieldDefinition, caseinsensetiveless<std::string>>;		
		IFDKeyMap ifdGMap;
		std::vector<cFieldDefinition> fdPFRGvec;

		//Sample instances
		cSampleBunch Bunch;
		std::vector<SampleId> Id;
		std::vector<cKeyVec<std::string, cFdVrnt, caseinsensetiveequal<std::string>>> AncFld;

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

		template <typename T>
		size_t value_size() {};
		template<> size_t value_size<double>() { return 1; };
		template<> size_t value_size<cdouble>() { return 2; };

		const AEM::SystemType& aem_system_type() const {
			return AEM::aem_system_type<AEMSystemClass>();
		};

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

		bool solve_thickness() const
		{
			return fdT.solve;
		};

		bool solve_conductivity() const {
			return fdC.solve;
		};

		bool solve_geometry_elementname(const std::string& gname) const {
			return ifdGMap.at(gname).solve;
		};

		bool solve_geometry_index(const size_t index) const {
			return ifdGMap.at(TDEmGeometry::element_name(index)).solve;
		};

		bool solve_geometry() const {
			if (nGeomParamPerSounding > 0) return true;
			return true;
		};

		bool solve_scalingfactor(const size_t sysi, const size_t ci) const {
			const cInvertibleFieldDefinition& ifd = get_ifd(sysi, ci, SCALEFACTOR);
			return ifd.solve;
		};

		bool solve_scalingfactors() const {
			if (nScalingParam > 0) return true;
			else return false;
		}

		cInvertibleFieldDefinition& get_ifd(const size_t sysi, const size_t ci, const std::string& key) {
			AEMComponentInversionInfo<RT>& C = SysInfo[sysi].CompInfo[ci];
			return C.get_ifd(key);
		};

		const cInvertibleFieldDefinition& get_ifd(const size_t sysi, const size_t ci, const std::string& key) const {
			const AEMComponentInversionInfo<RT>& C = SysInfo[sysi].CompInfo[ci];
			return C.get_ifd(key);
		};

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
			return (int)(si * nParamPerSounding + fdC.get_poffset() + li);
			//return (int)(si * nParamPerSounding + cOffset + li);
		}

		int tindex(const size_t& si, const size_t& li) const {
			if (solve_thickness() == false) {
				glog.errormsg(_SRC_, "Out of boundes in tindex().");
			}
			return (int)(si * nParamPerSounding + fdT.get_poffset() + li);
			//return (int)(si * nParamPerSounding + tOffset + li);
		}

		int gindex(const size_t& si, const std::string& gname) const {
			const cInvertibleFieldDefinition& ifd = ifdGMap.at(gname);
			int poffset = ifd.get_poffset();
			if (poffset >= 0) return (int)(si * nParamPerSounding + poffset);
			return -1;
		};

		int gindex(const size_t& si, const size_t& gi) const {
			const std::string gname = TDEmGeometry::element_name(gi);
			return gindex(si,gname);
		};

		int sfindex(const size_t& sysi, const size_t& ci) const {
			return get_ifd(sysi, ci, SCALEFACTOR).get_poffset();
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

			fdC = cInvertibleFieldDefinition(b, "Conductivity");
			if (nLayers > 1) {
				fdT = cInvertibleFieldDefinition(b, "Thickness");
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

		IFDKeyMap get_field_definitions_geometry(const cBlock& parent) const {
			IFDKeyMap gmap;
			for (size_t i = 0; i < TDEmGeometry::NELEM; i++) {
				std::string gname = TDEmGeometry::element_name(i);
				if (gmap.count(gname) == 0) {
					gmap[gname] = cInvertibleFieldDefinition(parent, gname);
				}
				else {
					std::string msg = strprint("Parameter %s has already been already added.", gname.c_str());
					glog.errormsg(_SRC_, msg);
				}
			}
			return gmap;
		};

		bool reconstruct_primary() const {
			for (size_t sysi = 0; sysi < SysInfo.size(); sysi++) {
				if (SysInfo[sysi].ReconstructPrimary == true) return true;
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
			//cOffset = 0;
			//tOffset = 0;

			//int poffset = 0;
			if (solve_conductivity()) {
				fdC.set_poffset((int)nParamPerSounding);
				nParamPerSounding += nLayers;
			}

			if (solve_thickness()) {
				fdT.set_poffset((int)nParamPerSounding);
				nParamPerSounding += nLayers - 1;
			}

			//Geometry params
			for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
				std::string gname = TDEmGeometry::element_name(gi);
				cInvertibleFieldDefinition& ifd = ifdGMap.at(gname);
				if (ifd.solve) {
					ifd.set_poffset((int)nParamPerSounding);
					nParamPerSounding++;
					nGeomParamPerSounding++;
				}
				else {
					ifd.set_poffset(-1);
				}
			};

			//Scaling params
			for (size_t sysi = 0; sysi < SysInfo.size(); sysi++) {
				for (size_t ci = 0; ci < NCOMP; ci++) {
					cInvertibleFieldDefinition& ifd = get_ifd(sysi, ci, SCALEFACTOR);
					if(ifd.solve) {
						ifd.set_poffset((int)(nParamPerSounding * nSoundings + nScalingParam));
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

			if (nScalingParam == 0) {
				LCrefs.alpha = 0.0;
			}

			if (nSoundings == 1) {
				LClatc.alpha = 0.0;
				LClatg.alpha = 0.0;
			}
		};

		void setup_parameter_bounds() {
			const static double ud = undefinedvalue<double>();
			Param_Min.resize(nParam);
			Param_Max.resize(nParam);
			Param_Min.setConstant(ud);
			Param_Max.setConstant(ud);

			if (fdC.bound()) {
				for (size_t si = 0; si < nSoundings; si++) {
					const EarthStore& e = EStore[si];
					for (size_t li = 0; li < nLayers; li++) {
						const int pi = cindex(si, li);
						Param_Min[pi] = std::log10(e.minval.conductivity[li]);
						Param_Max[pi] = std::log10(e.maxval.conductivity[li]);
					}
				}
			}

			if (fdT.bound()) {
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
					const cInvertibleFieldDefinition& e = ifdGMap.at(ename);
					if (e.bound()) {
						const int pi = gindex(si, ename);
						Param_Min[pi] = g.minval[ename];
						Param_Max[pi] = g.maxval[ename];
					}
				}
			}

			for (size_t sysi = 0; sysi < SysInfo.size(); sysi++) {
				for (size_t ci = 0; ci < NCOMP; ci++) {
					const cInvertibleFieldDefinition& ifd = get_ifd(sysi, ci, SCALEFACTOR);
					if (ifd.bound()) {
						int pi = ifd.get_poffset();
						if (ifd.solve) {
							const AEMComponentInversionInfo<RT>& C = SysInfo[sysi].CompInfo[ci];
							Param_Min[pi] = C.sfStore.minval;
							Param_Max[pi] = C.sfStore.maxval;
						}
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
			if (nGeomParamPerSounding <= 0)return;

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

		void initialise_Ws() {
			cLinearConstraint& C = LCrefs;
			C.W = Matrix::Zero(nParam, nParam);
			if (solve_scalingfactors() == false)return;

			double s = C.alpha / (double)(nScalingParam);
			for (size_t sysi = 0; sysi < SysInfo.size(); sysi++) {
				for (size_t ci = 0; ci < NCOMP; ci++) {
					const int pi = sfindex(sysi, ci);
					if (pi >= 0) {
						C.W(pi, pi) = s / (RefParamStd[pi] * RefParamStd[pi]);
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

		Vector CableLengths(const Vector& m) {
			Vector cablelength(nSoundings);
			const std::vector<TDEmGeometry> gv = get_geometry(m);
			for (size_t si = 0; si < nSoundings; si++) {
				const TDEmGeometry& g = gv[si];
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

		Vector BoundsConstraint_forward(const Vector& m) {
			cNonLinearConstraint& C = NLCbounds;
			Vector predicted = Vector::Zero(nParam);
			if (C.alpha == 0.0) return predicted;
			for (size_t pi = 0; pi < nParam; pi++) {
				const double& L = Param_Min[pi];
				const double& U = Param_Max[pi];
				const double& N = 0.5;
				double x = m[pi];
				predicted[pi] = log_barrier(L, U, N, x);
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
				C.J(pi, pi) = log_barrier_deriv(L, U, N, x);
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
			if (aem_system_type() == AEM::SystemType::TimeDomain) set_fftw_lock();
			std::vector<cBlock> B = Control.findblocks("EMSystem");
			const size_t nsys = B.size();
			for (size_t sysi = 0; sysi < nsys; sysi++) {
				cBlock& b = B[sysi];
				fs::path stmfile;
				if (b.getvalue("SystemFile", stmfile)) {
					glog.logmsg(0, "Reading AEM system file %s\n", stmfile.string().c_str());
					SysInfo.emplace_back(AEMSystemInversionInfo<AEMSystemClass, RT>(b, nSoundings));
					SysInfo[sysi].set_units(IM.get());
				}
				else glog.errormsg(_SRC_, "No AWM 'SystemFile' is specified.\n");
			}
			if (aem_system_type() == AEM::SystemType::TimeDomain) unset_fftw_lock();
		}

		void setup_data() {
			nAllData = 0;
			const size_t nsys = nSystems();
			const size_t vsize = value_size<RT>();
			_vindex_.resize(nSoundings);
			for (size_t si = 0; si < nSoundings; si++) {
				_vindex_[si].resize(nsys);
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					_vindex_[si][sysi].resize(XZAMP+1);//4 because of xzinversion				
					for (size_t ci = 0; ci < (XZAMP+1); ci++) {
						_vindex_[si][sysi][ci].resize(SysInfo[sysi].nwindows);
						for (size_t wi = 0; wi < SysInfo[sysi].nwindows; wi++) {
							_vindex_[si][sysi][ci][wi] = -1;
						}
					}
				}
			}

			int vi = 0;
			for (size_t si = 0; si < nSoundings; si++) {
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					auto& S = SysInfo[sysi];
					if (S.InvertXZAmplitude) {
						nAllData += S.nwindows * vsize;
						for (size_t wi = 0; wi < S.nwindows; wi++) {
							_vindex_[si][sysi][XZAMP][wi] = vi;
							vi++;
						}

						if (S.CompInfo[YCOMP].Use) {
							nAllData += S.nwindows * vsize;
							for (size_t wi = 0; wi < S.nwindows; wi++) {
								_vindex_[si][sysi][YCOMP][wi] = vi;
								vi++;
							}
						}
					}
					else {
						for (size_t ci = 0; ci < NCOMP; ci++) {
							auto& c = S.CompInfo[ci];
							if (c.Use) {
								nAllData += S.nwindows * vsize;
								for (size_t wi = 0; wi < S.nwindows; wi++) {
									_vindex_[si][sysi][ci][wi] = vi;
									vi++;
								}
							}
						}
					}
				}
			}
		}

		bool initialise_bunch_data() {
			Vector obs(nAllData);
			Vector err(nAllData);
			Vector pred(nAllData);

			RT ov, ev;
			const size_t nsys = nSystems();
			for (size_t si = 0; si < nSoundings; si++) {
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					auto& S = SysInfo[sysi];
					AEMSystem<RT>& A = *S.System;
					if (S.ReconstructPrimary) {
						TDEmVectorResponse<RT> P = A.forward_model_primary_field(GStore[si].pfr);
						if (S.CompInfo[XCOMP].Use) S.CompInfo[XCOMP].data[si].P = P[XCOMP];
						if (S.CompInfo[YCOMP].Use) S.CompInfo[YCOMP].data[si].P = P[YCOMP];
						if (S.CompInfo[ZCOMP].Use) S.CompInfo[ZCOMP].data[si].P = P[ZCOMP];
					}

					if (S.InvertXZAmplitude) {
						for (size_t wi = 0; wi < S.nwindows; wi++) {
							//XZ Amplitude
							int vi = vindex(si, sysi, XZAMP, wi);
							RT X = S.CompInfo[XCOMP].data[si].S[wi];
							RT Z = S.CompInfo[ZCOMP].data[si].S[wi];
							if (S.InvertTotalField) {
								X += S.CompInfo[XCOMP].data[si].P[wi];
								Z += S.CompInfo[ZCOMP].data[si].P[wi];
							}
							ov = AEM::hypot(X, Z);
							set_data_vector(obs, vi, ov);


							const RT& Xerr = S.CompInfo[XCOMP].data[si].E[wi];
							const RT& Zerr = S.CompInfo[ZCOMP].data[si].E[wi];

							// Todo fix this 
							if (ov == 0.0) {
								ev = AEM::hypot(Xerr, Zerr);
							}
							else {
								ev = AEM::ewise_div(AEM::hypot(AEM::ewise_mul(X, Xerr), AEM::ewise_mul(Z, Zerr)), ov);
							}
							set_data_vector(err, vi, ev);


							//Y Comp
							if (S.CompInfo[YCOMP].Use) {
								int vi = vindex(si, sysi, YCOMP, wi);
								ov = S.CompInfo[YCOMP].data[si].S[wi];
								if (S.InvertTotalField) {
									ov += S.CompInfo[YCOMP].data[si].P[wi];
								}
								ev = S.CompInfo[YCOMP].data[si].E[wi];
								set_data_vector(obs, vi, ov);
								set_data_vector(err, vi, ev);
							}

						}
					}
					else {
						for (size_t ci = 0; ci < NCOMP; ci++) {
							if (S.CompInfo[ci].Use) {
								const SoundingData<RT>& d = S.CompInfo[ci].data[si];

								for (size_t wi = 0; wi < S.nwindows; wi++) {
									int vi = vindex(si, sysi, ci, wi);

									if (S.InvertTotalField) ov = d.T[wi];
									else ov = d.S[wi];

									ev = d.E[wi];
									set_data_vector(obs, vi, ov);
									set_data_vector(err, vi, ev);
								}
							}
						}
					}
				}
			}

			if (ErrorAddition > 0.0) {
				for (size_t k = 0; k < err.size(); k++) {
					err[k] = err[k] + ErrorAddition;
				}
			}

			//Work out indices to be culled
			ActiveData.clear();
			for (size_t i = 0; i < nAllData; i++) {
				if (!isnull(obs[i]) && !isnull(err[i])) ActiveData.push_back(i);
			}
			nData = ActiveData.size();

			if (nData != nAllData) {
				size_t ncull = nAllData - nData;
				OutputMessage += strprint(", %d null data/noise were culled", (int)ncull);
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

		void initialise_bunch_parameters() {

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
					const std::string& gname = TDEmGeometry::element_name(gi);
					const int pi = gindex(si, gi);
					if (pi >= 0) {
						RefParam[pi] = g.refval[gname];
						RefParamStd[pi] = g.refvalstd[gname];
					}
				}

				//Scaling params
				if (solve_scalingfactors()) {
					for (size_t sysi = 0; sysi < SysInfo.size(); sysi++) {
						for (size_t ci = 0; ci < NCOMP; ci++) {
							const AEMComponentInversionInfo<RT>& C = SysInfo[sysi].CompInfo[ci];
							const int pi = sfindex(sysi, ci);
							if (pi >= 0) {
								RefParam[pi] = C.sfStore.refval;
								RefParamStd[pi] = C.sfStore.refvalstd;
							}
						}
					}
				}
			}
		}

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

		Vector elementwise_bound_restrict(const Vector& m_old, const Vector& m_new, Vector& dm) {
			if (fdC.bound()) {
				for (size_t si = 0; si < nSoundings; si++) {
					const EarthStore& e = EStore[si];
					for (size_t li = 0; li < nLayers; li++) {
						const int pindex = cindex(si, li);
						const double lmin = std::log10(e.minval.conductivity[li]);
						const double lmax = std::log10(e.maxval.conductivity[li]);
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
					const EarthStore& e = EStore[si];
					for (size_t li = 0; li < nLayers - 1; li++) {
						const int pindex = tindex(si, li);
						const double lmin = std::log10(e.minval.thickness[li]);
						const double lmax = std::log10(e.maxval.thickness[li]);
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
				GeometryStore& g = GStore[si];
				for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
					const std::string ename = TDEmGeometry::element_name(gi);
					const cInvertibleFieldDefinition& e = ifdGMap.at(ename);
					if (e.bound()) {
						const int pi = gindex(si, ename);
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

		Vector parameter_change(const double& lambda, const Vector& m_old, const Vector& pred) {
			const Vector m_new = solve_linear_system(lambda, m_old, pred);
			Vector dm = m_new - m_old;
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

		std::vector<TDEmGeometry> get_geometry(const Vector& parameters) {
			std::vector<TDEmGeometry> gv(nSoundings);
			for (size_t si = 0; si < nSoundings; si++) {
				gv[si] = GStore[si].input;
				for (int gi = 0; gi < TDEmGeometry::NELEM; gi++) {
					const std::string& gname = TDEmGeometry::element_name(gi);
					const int pi = gindex(si, gname);
					if (pi >= 0) {
						gv[si][gname] = parameters[pi];
					}
				}
			}
			return gv;
		}

		Vec3d get_scalefactors(const size_t sysi, const Vector& parameters) const {
			Vec3d sf;
			for (int ci = 0; ci < NCOMP; ci++) {
				sf[ci] = 1.0;
				const int pi = sfindex(sysi, ci);
				if (pi >= 0) sf[ci] = parameters[pi];
			}
			return sf;
		}

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
			Vector pred_all(nAllData);
			Matrix J_all;
			if (computederivatives) {
				J_all.resize(nAllData, nParam);
				J_all.setZero();
			}

			const size_t nsys = nSystems();
			std::vector<Earth1D> ev = get_earth(parameters);
			std::vector<TDEmGeometry> gv = get_geometry(parameters);
			for (size_t sysi = 0; sysi < nsys; sysi++) {
				const auto& S = SysInfo[sysi];
				AEMSystem<RT>& A = *S.System;
				const size_t& nw = A.nWindows();

				Vec3d scalefactors(1.0, 1.0, 1.0);
				const bool solvesf = solve_scalingfactors();
				if (solvesf) scalefactors = get_scalefactors(sysi, parameters);

				for (size_t si = 0; si < nSoundings; si++) {
					const Earth1D& e = ev[si];
					const TDEmGeometry& g = gv[si];

					// R is a reference to the Work Response struct
					const TDEmResponse<RT>& R = A.forward_model(e, g);
					
					TDEmVectorResponse<RT> FM;
					if (S.InvertTotalField) FM = R.totalfield();
					else FM = R.S;
					//std::cout << FM << std::endl;

					if (solvesf) FM.scale_components(scalefactors);

					TDEmScalarResponse<RT> XZFM;
					if (S.InvertXZAmplitude) {
						XZFM = FM.xzamp();
					}

					// Predicted
					if (S.InvertXZAmplitude) {
						for (size_t wi = 0; wi < nw; wi++) {
							int di = vindex(si, sysi, XZAMP, wi);
							//pred_all[di] = XZFM[wi];
							set_data_vector(pred_all, di, XZFM[wi]);
							if (S.CompInfo[YCOMP].Use) {
								di = vindex(si, sysi, YCOMP, wi);
								//pred_all[di] = FM[wi][YCOMP];
								set_data_vector(pred_all, di, FM[wi][YCOMP]);
							}
						}
					}
					else {
						for (size_t ci = 0; ci < NCOMP; ci++) {
							if (S.CompInfo[ci].Use) {
								for (size_t wi = 0; wi < nw; wi++) {
									int di = vindex(si, sysi, ci, wi);
									//pred_all[di] = FM[ci][wi];
									set_data_vector(pred_all, di, FM[ci][wi]);
								}
							}
						}
					}

					// Jacobian
					if (computederivatives) {
						TDEmVectorResponse<RT> DRV;

						// Scale factor derivatives
						if (solvesf) {
							for (size_t ci = 0; ci < NCOMP; ci++) {
								if (S.CompInfo[ci].Use) {
									const int pindex = sfindex(sysi, ci);
									if (pindex >= 0) {
										DRV = FM;
										// Here filling with the forward itself as derivative w.r.t scale factor param is the forward model itself
										// But zero for other components
										Vec3d f(0, 0, 0);
										f[ci] = 1.0;
										DRV.scale_components(f);
										fillMatrixColumn(J_all, si, sysi, pindex, FM, XZFM, DRV);
									}
								}
							}
						}

						if (solve_conductivity()) {
							for (size_t li = 0; li < nLayers; li++) {
								const int pindex = cindex(si, li);
								A.derivative(CalculationType(CMode::DC, li));
								if (S.InvertTotalField) DRV = R.totalfield();
								else DRV = R.S;
								//multiply by natural log(10) as parameters are in logbase10 units
								const double f = log(10.0) * e.conductivity[li];
								DRV *= f;
								if (solvesf) DRV.scale_components(scalefactors);
								//std::cout << FM << std::endl;
								//std::cout << DRV << std::endl;
								fillMatrixColumn(J_all, si, sysi, pindex, FM, XZFM, DRV);
							}
						}

						if (solve_thickness()) {
							for (size_t li = 0; li < nLayers - 1; li++) {
								const int pindex = tindex(si, li);
								A.derivative(CalculationType(CMode::DT, li));
								if (S.InvertTotalField) DRV = R.totalfield();
								else DRV = R.S;
								//multiply by natural log(10) as parameters are in logbase10 units
								double f = log(10.0) * e.thickness[li];
								DRV *= f;
								if (solvesf) DRV.scale_components(scalefactors);
								fillMatrixColumn(J_all, si, sysi, pindex, FM, XZFM, DRV);
							}
						}

						if (FreeGeometry) {
							if (solve_geometry_elementname("tx_height")) {
								const size_t pindex = gindex(si, "tx_height");
								A.derivative(CalculationType(CMode::DTX_HEIGHT));
								if (S.InvertTotalField) DRV = R.totalfield();
								else DRV = R.S;
								if (solvesf) DRV.scale_components(scalefactors);
								fillMatrixColumn(J_all, si, sysi, pindex, FM, XZFM, DRV);
							}

							if (solve_geometry_elementname("txrx_dx")) {
								const size_t pindex = gindex(si, "txrx_dx");
								A.derivative(CalculationType(CMode::DX));
								if (S.InvertTotalField) DRV = R.totalfield();
								else DRV = R.S;
								if (solvesf) DRV.scale_components(scalefactors);
								//std::cout << FM << std::endl;
								//std::cout << DRV << std::endl;
								fillMatrixColumn(J_all, si, sysi, pindex, FM, XZFM, DRV);
							}

							if (solve_geometry_elementname("txrx_dy")) {
								const size_t pindex = gindex(si, "txrx_dy");
								A.derivative(CalculationType(CMode::DY));
								if (S.InvertTotalField) DRV = R.totalfield();
								else DRV = R.S;
								if (solvesf) DRV.scale_components(scalefactors);
								fillMatrixColumn(J_all, si, sysi, pindex, FM, XZFM, DRV);
							}

							if (solve_geometry_elementname("txrx_dz")) {
								const size_t pindex = gindex(si, "txrx_dz");
								A.derivative(CalculationType(CMode::DZ));
								if (S.InvertTotalField) DRV = R.totalfield();
								else DRV = R.S;
								if (solvesf) DRV.scale_components(scalefactors);
								//std::cout << FM << std::endl;
								//std::cout << DRV << std::endl;
								fillMatrixColumn(J_all, si, sysi, pindex, FM, XZFM, DRV);
							}

							if (solve_geometry_elementname("rx_roll")) {
								const size_t pindex = gindex(si, "rx_roll");
								DRV = A.derivative(CalculationType(CMode::DRX_ROLL), g, FM);
								if (solvesf) DRV.scale_components(scalefactors);
								fillMatrixColumn(J_all, si, sysi, pindex, FM, XZFM, DRV);
							}

							if (solve_geometry_elementname("rx_pitch")) {
								const size_t pindex = gindex(si, "rx_pitch");
								DRV = A.derivative(CalculationType(CMode::DRX_PITCH), g, FM);
								if (solvesf) DRV.scale_components(scalefactors);
								fillMatrixColumn(J_all, si, sysi, pindex, FM, XZFM, DRV);
							}

							if (solve_geometry_elementname("rx_yaw")) {
								const size_t pindex = gindex(si, "rx_yaw");
								DRV = A.derivative(CalculationType(CMode::DRX_YAW), g, FM);
								if (solvesf) DRV.scale_components(scalefactors);
								fillMatrixColumn(J_all, si, sysi, pindex, FM, XZFM, DRV);
							}
						}
					}
				}
			}
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
				auto& S = SysInfo[sysi];
				S.predicted.resize(nSoundings);

				AEMSystem<RT>& A = *S.System;
				const size_t& nw = A.nWindows();
				for (size_t si = 0; si < nSoundings; si++) {
					const Earth1D& e = ev[si];
					const TDEmGeometry& g = gv[si];
					S.predicted[si] = A.forward_model(e, g);
				}
			}
		};

		void fillMatrixColumn(Matrix& M, const size_t& si, const size_t& sysi, const size_t& pindex,
			const TDEmVectorResponse<RT>& FM,
			const TDEmScalarResponse<RT>& XZFM,
			const TDEmVectorResponse<RT>& DRV) {
			const auto& S = SysInfo[sysi];
			const AEMSystem<RT>& A = *S.System;
			const size_t& nw = A.nWindows();
			if (S.InvertXZAmplitude) {
				// dr/dp = (x/r)dx/dp + (y/r)dy/dp
				for (size_t wi = 0; wi < nw; wi++) {
					const int di = vindex(si, sysi, XZAMP, wi);
					//M(di, pindex) = (FM[XCOMP][wi] * DRV[XCOMP][wi] + FM[ZCOMP][wi] * DRV[ZCOMP][wi]) / XZFM[wi];
					set_mat(M, di, pindex, (FM[XCOMP][wi] * DRV[XCOMP][wi] + FM[ZCOMP][wi] * DRV[ZCOMP][wi]) / XZFM[wi]);
				}

				if (S.CompInfo[YCOMP].Use) {
					for (size_t wi = 0; wi < nw; wi++) {
						const int di = vindex(si, sysi, YCOMP, wi);
						//M(di, pindex) = DRV[YCOMP][wi];
						set_mat(M, di, pindex, DRV[YCOMP][wi]);
					}
				}
			}
			else {
				for (size_t ci = 0; ci < NCOMP; ci++) {
					if (S.CompInfo[ci].Use) {
						for (size_t wi = 0; wi < nw; wi++) {
							const int di = vindex(si, sysi, ci, wi);
							//M(di, pindex) = DRV[ci][wi];
							set_mat(M, di, pindex, DRV[ci][wi]);
						}
					}
				}
			}
		}

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

			if (solve_scalingfactors()) {
				for (size_t sysi = 0; sysi < SysInfo.size(); sysi++) {
					for (size_t ci = 0; ci < NCOMP; ci++) {
						const cInvertibleFieldDefinition& ifd = get_ifd(sysi, ci, SCALEFACTOR);
						if (ifd.solve) {
							auto& C = SysInfo[sysi].CompInfo[ci];
							IM->read(ifd.refval, C.sfStore.refval);
							IM->read(ifd.refvalstd, C.sfStore.refvalstd);
							IM->read(ifd.minval, C.sfStore.minval);
							IM->read(ifd.maxval, C.sfStore.maxval);
						}
					}
				}
			}

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
			status = IM->read(fdC.input, e.refval.conductivity, nLayers); if (status == false) readstatus = false;
			if (solve_conductivity()) {
				status = IM->read(fdC.refval, e.refval.conductivity, nLayers); if (status == false) readstatus = false;
				status = IM->read(fdC.refvalstd, e.refvalstd.conductivity, nLayers); if (status == false) readstatus = false;
				status = IM->read(fdC.minval, e.minval.conductivity, nLayers); if (status == false) readstatus = false;
				status = IM->read(fdC.maxval, e.maxval.conductivity, nLayers); if (status == false) readstatus = false;
			}

			status = IM->read(fdT.input, e.refval.thickness, nLayers - 1); if (status == false) readstatus = false;
			if (solve_thickness()) {
				status = IM->read(fdT.refval, e.refval.thickness, nLayers - 1); if (status == false) readstatus = false;
				status = IM->read(fdT.refvalstd, e.refvalstd.thickness, nLayers - 1); if (status == false) readstatus = false;
				status = IM->read(fdT.minval, e.minval.thickness, nLayers - 1); if (status == false) readstatus = false;
				status = IM->read(fdT.maxval, e.maxval.thickness, nLayers - 1); if (status == false) readstatus = false;
			}
			e.sanity_check();

			status = IM->read(fdT.maxval, e.maxval.thickness, nLayers - 1); if (status == false) readstatus = false;

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
				const cInvertibleFieldDefinition ge = map.at(ename);
				bool inpstatus = IM->read(ge.input, gstore.input[gi]);
				bool refstatus = IM->read(ge.refval, gstore.refval[gi]);

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

				if (ge.solve) {
					bool stdstatus = IM->read(ge.refvalstd, gstore.refvalstd[gi]);
					if (stdstatus == false) {
						std::ostringstream msg;
						msg << "No 'Std' defined for " << ename << std::endl;
						glog.errormsg(_SRC_, msg.str());
					}

					bool minstatus = IM->read(ge.minval, gstore.minval[gi]);
					bool maxstatus = IM->read(ge.maxval, gstore.maxval[gi]);
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
		}

		void read_system_data(size_t& sysindex, const size_t& soundingindex)
		{
			AEMSystemInversionInfo<AEMSystemClass, RT>& S = SysInfo[sysindex];
			S.CompInfo[XCOMP].readdata(IM, soundingindex);
			S.CompInfo[YCOMP].readdata(IM, soundingindex);
			S.CompInfo[ZCOMP].readdata(IM, soundingindex);
		}

		void dump_first_iteration() {

			const std::string dp = dumppath().string();
			makedirectory(dumppath());

			const size_t si = Bunch.master_index();
			GeometryStore& g = GStore[si];
			EarthStore& e = EStore[si];
			SampleId& id = Id[si];

			write(Obs, dp + "observed.dat");
			write(Err, dp + "observed_std.dat");

			g.refval.write(dp + "geometry_start.dat");
			e.refval.write(dp + "earth_start.dat");

			g.refval.write(dp + "geometry_ref.dat");
			e.refval.write(dp + "earth_ref.dat");

			g.refvalstd.write(dp + "geometry_std.dat");
			e.refvalstd.write(dp + "earth_std.dat");

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

					Vector g;
					forwardmodel_and_jacobian(CIS.param, g, J);
					if (CIS.iteration == 0) {
						CIS.lambda = 1e8;
						if (OutputOpt.Dump) {
							dump_first_iteration();
							dump_iteration(CIS);
						}
					}

					const double targetphid = std::max(CIS.phid * 0.7, MinimumPhiD);
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
			bool readstatus = true;
			int paralleljob = 0;
			do {
				int record = ((int)StartRecord - 1) + paralleljob * (int)IM->subsamplerate();
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
						glog.logmsg(0, s.str());
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

		void write_result(const int& pointindex) {
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
				if (fi == 0) c.width++;//First column after unique_id needs to have extra space in case it fills the width to enfore space between columns
				OM->writevrnt(pi, fdv.vnt, c);
			}

			//Geometry Input
			bool invertedfieldsonly = false;
			for (size_t i = 0; i < TDEmGeometry::NELEM; i++) {
				if (invertedfieldsonly && solve_geometry_index(i) == false)continue;
				OM->writefield(pi, GStore[si].input[i], "input_" + GStore[si].input.element_name(i), "Input " + GStore[si].input.description(i), GStore[si].input.units(i), 1, ST_FLOAT, DN_NONE, 'F', 9, 2);
			}

			//Geometry Modelled		
			const TDEmGeometry& g = GStore[si].invmodel;
			invertedfieldsonly = true;
			for (size_t gi = 0; gi < TDEmGeometry::NELEM; gi++) {
				if (invertedfieldsonly && solve_geometry_index(gi) == false)continue;
				OM->writefield(pi, g[gi], "inverted_" + g.element_name(gi), "Inverted " + g.description(gi), g.units(gi), 1, ST_FLOAT, DN_NONE, 'F', 9, 2);
			}

			//ndata
			OM->writefield(pi, nData, "ndata", "Number of data in inversion", UNITLESS, 1, ST_UINT, DN_NONE, 'I', 4, 0);

			//Scaling factors
			if (solve_scalingfactors()) {
				const size_t nsys = nSystems();
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					Vec3d sf = get_scalefactors(sysi, m);
					for (size_t ci = 0; ci < NCOMP; ci++) {
						const AEMComponentInversionInfo<RT>& C = SysInfo[sysi].CompInfo[ci];
						if (C.Use && solve_scalingfactor(sysi, ci)) {
							std::string comp = C.Name;
							std::string fname = "scalingfactor" + strprint("_EMSystem_%d_", (int)sysi + 1) + comp;
							std::string fdesc = "Scaling factor" + strprint(" EMSystem %d ", (int)sysi + 1) + comp + "-component";
							OM->writefield(pi, sf[ci], fname, fdesc, UNITLESS, 1, ST_FLOAT, DN_NONE, 'F', 6, 3);
						}
					}
				}
			}

			//Earth	
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
						const std::string& gname = g.element_name(gi);
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
					const std::string& gname = g.element_name(gi);
					std::string name = "inverted_" + gname + "_uncertainty";
					std::string desc = g.description(gi) + " parameter uncertainty";
					OM->writefield(pi, pu[gindex(si, gname)], name, desc, g.units(gi), 1, ST_FLOAT, DN_NONE, 'E', 15, 6);
				}
			}

			//ObservedData
			const size_t nsys = nSystems();
			cAsciiColumnFormat emfmt('E', 15, 6);
			if (OutputOpt.ObservedData) {
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					const auto& S = SysInfo[sysi];
					bool reconstructed_primaryfield_flag = S.ReconstructPrimary;//Only do this for the observed data but not for predicted data or noise
					for (size_t ci = 0; ci < NCOMP; ci++) {
						const SoundingData<RT>& d = S.CompInfo[ci].data[si];
						if (S.CompInfo[ci].Use) writeresult_emdata(pi, sysi, ci, "Observed", "Observed", S.units, emfmt, d.P, d.S, d.T);
					}
				}
			}


			//Noise Estimates
			if (OutputOpt.NoiseEstimates) {
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					const auto& S = SysInfo[sysi];
					for (size_t ci = 0; ci < NCOMP; ci++) {
						const SoundingData<RT>& d = S.CompInfo[ci].data[si];
						if (S.CompInfo[ci].Use) writeresult_emdata(pi, sysi, ci, "Noise", "Estimated noise", S.units, emfmt, d.P, d.E, d.T);
					}
				}
			}

			//PredictedData
			if (OutputOpt.PredictedData) {
				for (size_t sysi = 0; sysi < nsys; sysi++) {
					const auto& S = SysInfo[sysi];
					for (size_t ci = 0; ci < NCOMP; ci++) {
						const auto p = S.predicted[si].primary(ci);
						const auto s = S.predicted[si].secondary(ci);
						const auto t = p + s;
						if (S.CompInfo[ci].Use) writeresult_emdata(pi, sysi, ci, "Predicted", "Predicted", S.units, emfmt, p, s, t);
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
		std::vector<T> real(const std::vector<std::complex<T>>& cv) const {
			const size_t n = cv.size();
			std::vector<T> v(n);
			for (size_t i = 0; i < n; i++) {
				v[i] = cv[i].real();
			}
			return v;
		};

		template <typename T>
		std::vector<T> imaginary(const std::vector<std::complex<T>>& cv) const {
			const size_t n = cv.size();
			std::vector<T> v(n);
			for (size_t i = 0; i < n; i++) {
				v[i] = cv[i].imag();
			}
			return v;
		};

		template <typename T>
		void write_plain(const int& pointindex, const cOutputField& of, const std::vector<T>& values)
		{
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

		void writeresult_emdata(const int& pointindex, const size_t& sysindex, const size_t& compindex, const std::string& nameprefix, const std::string& descprefix, const std::string& units, const cAsciiColumnFormat& fmt, const std::vector<RT>& Primary, const std::vector<RT>& Secondary, const std::vector<RT>& Total) {
			const AEMSystemInversionInfo<AEMSystemClass, RT>& S = SysInfo[sysindex];
			const std::string compname = S.CompInfo[compindex].Name;

			const BinaryStorageType btype = ST_FLOAT;
			std::string dimensionname = "em_window";
			std::string sysname = nameprefix + strprint("_EMSystem_%d_", (int)sysindex + 1);
			std::string sysdesc = descprefix + strprint(" EMSystem %d ", (int)sysindex + 1);

			const int nbands = Secondary.size();

			//Primary field
			if (S.InvertTotalField) {
				std::string name = sysname + compname + "P";
				std::string desc = sysdesc + compname + "-component primary field";
				if (S.ReconstructPrimary) desc += " reconstructed from input geometry";
				cOutputField of(name, desc, units, nbands, btype, dimensionname, fmt);
				write_plain_or_complex(pointindex, of, Primary);
			}

			// Secondary field
			{
				std::string name = sysname + compname + "S";
				std::string desc = sysdesc + compname + "-component secondary field";
				cOutputField of(name, desc, units, nbands, btype, dimensionname, fmt);
				write_plain_or_complex(pointindex, of, Secondary);
			}

			// Toatal field
			if (S.InvertTotalField) {
				std::string name = sysname + compname + "T";
				std::string desc = sysdesc + compname + "-component total field";
				cOutputField of(name, desc, units, nbands, btype, dimensionname, fmt);
				write_plain_or_complex(pointindex, of, Total);
			}
		};
	};
};//End namespace


