/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include <stdexcept>
#include <complex>

#include "eigen_utils.hpp"
#include "fftw3.h"
#include "file_utils.hpp"
#include "vector_utils.hpp"
#include "general_utils.hpp"
#include "blocklanguage.hpp"

#include "earth1d.hpp"
#include "fixed_point_spline.hpp"
#include "aem_coredefs.hpp"
#include "tdemgeometry.hpp"

namespace AEM {
	
	class FFTWPlanWrapper {

	private:
		fftw_plan Plan = nullptr;

	public:

		// Default constructor
		FFTWPlanWrapper() {
			Plan = nullptr;
		}

		FFTWPlanWrapper(const fftw_plan plan) {
			setplan(plan);
		}

		// Move constructor
		FFTWPlanWrapper(FFTWPlanWrapper&& other) noexcept
			: Plan(other.Plan)
		{
			other.Plan = nullptr;
		};

		// Copy assignment operator
		FFTWPlanWrapper& operator=(FFTWPlanWrapper& other) noexcept {
			setplan(other.Plan);
			other.Plan = nullptr;
			return *this;
		};

		~FFTWPlanWrapper() {
			destroy();
		};

		void setplan(const fftw_plan plan) {
			Plan = plan;
		}

		void destroy() const {
			if (Plan) {
				fftw_destroy_plan(Plan);
			}
		}

		void execute() const {
			fftw_execute(Plan);
		}

		void print() const {
			fftw_print_plan(Plan);
		}

	};

	class TDEmScalarResponse {

	private:
		std::vector<double> v;

	public:

		TDEmScalarResponse() {};

		TDEmScalarResponse(const size_t& nwindows) {
			resize(nwindows);
		}

		inline const size_t size() const { return v.size(); }

		double& operator[](const size_t& i) {
			return v[i];
		}

		double operator[](const size_t& i) const {
			return v[i];
		}

		TDEmScalarResponse& operator+=(const TDEmScalarResponse& rhs) {
			v += rhs.v;
			return *this;
		}

		TDEmScalarResponse& operator*=(const double& rhs) {
			v *= rhs;
			return *this;
		}

	private:

		void  resize(const size_t nwindows) {
			v.resize(nwindows);
		}
	};
	
	class TDEmVectorResponse {

		size_t nwindows=0;
		std::array<std::vector<double>, 3> v;

	public:

		TDEmVectorResponse(const size_t _nwindows = 0) { nwindows = _nwindows; };

		inline const size_t size() const { return nwindows; }

		void resize(const size_t _nwindows) {
			nwindows = _nwindows;
			v[0].resize(nwindows);
			v[1].resize(nwindows);
			v[2].resize(nwindows);
		};

		Vec3d get_vec3d(const size_t window) const {
			return Vec3d(v[0][window], v[1][window], v[2][window]);
		};

		void set_vec3d(const size_t window, const Vec3d& vec) {
			v[0][window] = vec[0];
			v[1][window] = vec[1];
			v[2][window] = vec[2];
		};

		std::vector<double>& operator[](const size_t& component) {
			return v[component];
		}

		const std::vector<double>& operator[](const size_t& component) const {
			return v[component];
		}

		TDEmVectorResponse& operator+=(const TDEmVectorResponse& rhs) {
			v[0] += rhs.v[0];
			v[1] += rhs.v[1];
			v[2] += rhs.v[2];
			return *this;
		}

		TDEmVectorResponse& operator*=(const double& s) {
			v[0] *= s;
			v[1] *= s;
			v[2] *= s;
			return *this;
		}

		double& operator()(const size_t& component, const size_t& window) {
			assert(component < NCOMP);
			assert(window < size());
			return v[component][window];
		}

		void scale_components(const Vec3d& scalefactors) {
			const size_t nw = v.size();
			v[0] *= scalefactors[0];
			v[1] *= scalefactors[1];
			v[2] *= scalefactors[2];
		};

		TDEmScalarResponse xzamp() {
			TDEmScalarResponse r(nwindows);
			for (size_t i = 0; i < nwindows; i++) {
				r[i] = std::hypot(v[XCOMP][i], v[ZCOMP][i]);
			}
			return r;
		};

	private:

	};

	class TDEmResponse {
	
	public:
		TDEmVectorResponse P;
		TDEmVectorResponse S;

		TDEmResponse() {};

		TDEmResponse(const size_t& _nwindows) {
			resize(_nwindows);
		}

		void resize(const size_t& _nwindows) {
			P.resize(_nwindows);
			S.resize(_nwindows);
		}

		const size_t& size() const {
			return S.size();
		}

		const double primary(const size_t& component) const {
			assert(component < NCOMP);
			return P[component][0];
		}

		const double secondary(const size_t& component, const size_t& window) const {
			assert(component < NCOMP);
			assert(window < size());
			return S[component][window];
		};

		const std::vector<double> secondary(const size_t& component) const {
			assert(component < NCOMP);
			return S[component];
		};

		TDEmVectorResponse totalfield() const {
			TDEmVectorResponse T = S;
			T += P;
			return T;
		};
	};

	class ComponentWorkStore {

	public:

		std::vector<double> IR_discrete_real;// Real impulse response discrete frequency nodes
		std::vector<double> IR_discrete_imag;// Imaginary impulse response discrete frequency nodes
		std::vector<cdouble> IR_splined;// Complex splines impulse response

		ComponentWorkStore() {};

		void resize(const size_t nnodes, const size_t nfftfreq, const size_t nwindows) {
			IR_discrete_real.resize(nnodes);
			IR_discrete_imag.resize(nnodes);
			IR_splined.resize(nfftfreq);
		};
	};

	class ModellingOptions {

	public:
		enum class OutputType { BFIELD, DBDT };
		enum class NormalizationType { NONE, PPM, PPM_PEAKTOPEAK };

		size_t FrequenciesPerDecade = 6;
		size_t NumAbscissa = 17;
		OutputType OutputType = OutputType::DBDT;
		NormalizationType NormalisationType = NormalizationType::NONE;
		double XOutputScaling = 1.0;
		double YOutputScaling = 1.0;
		double ZOutputScaling = 1.0;
		double ModellingLoopRadius = 0.0;

		bool SaveDiagnosticFiles = false;

		ModellingOptions() {};

		ModellingOptions(const cBlock& b) {
			read_modelling_options(b);
		};

		void read_modelling_options(const cBlock& b) {
			ModellingLoopRadius = b.getdoublevalue("ModellingLoopRadius");
			if (!isdefined(ModellingLoopRadius)) {
				ModellingLoopRadius = 0.0;
			}

			std::string ot = b.getstringvalue("OutputType");
			if (strcasecmp(ot, "B") == 0) {
				OutputType = OutputType::BFIELD;
			}
			else if (strcasecmp(ot, "dB/dt") == 0) {
				OutputType = OutputType::DBDT;
			}
			else {
				glog.errormsg(_SRC_, "OutputType %s unknown (must be one of \"B\" or \"dB/dt\")\n", ot.c_str());
			}

			FrequenciesPerDecade = b.getsizetvalue("FrequenciesPerDecade");
			if (FrequenciesPerDecade < 5) {
				glog.warningmsg(_SRC_, "It is wise to use at least 5 frequencies per decade\n");
			}

			NumAbscissa = b.getsizetvalue("NumberOfAbsiccaInHankelTransformEvaluation");
			if (NumAbscissa < 17) {
				glog.warningmsg(_SRC_, "It is wise to use at least 17 Absicca for integrating the Hankel Transforms");
			}

			std::string n = b.getstringvalue("SecondaryFieldNormalisation");
			if (strcasecmp(n, "None") == 0) {
				NormalisationType = NormalizationType::NONE;
			}
			else if (strcasecmp(n, "PPM") == 0) {
				NormalisationType = NormalizationType::PPM;
			}
			else if (strcasecmp(n, "PPMPEAKTOPEAK") == 0) {
				NormalisationType = NormalizationType::PPM_PEAKTOPEAK;
			}
			else {
				glog.errormsg(_SRC_, "Normalisation %s unknown (must be one of \"None,PPM,PPMPEAKTOPEAK\")\n", n.c_str());
			}

			XOutputScaling = b.getdoublevalue("XOutputScaling");
			YOutputScaling = b.getdoublevalue("YOutputScaling");
			ZOutputScaling = b.getdoublevalue("ZOutputScaling");


			SaveDiagnosticFiles = b.getboolvalue("SaveDiagnosticFiles");

		}

	};

	class AEMSystem {

	protected:
		std::string SystemName;
		std::string SystemType;
		cBlock STM;
		LEModeller LEM;

	public:

		const cBlock& system_descriptor_block() const { return STM; };
		LEModeller& lem() { return LEM; };

	};

	class TDEmSystem : public AEMSystem {

	private:
		FixedPointSpline<double> FrequencySpliner;
		FFTWPlanWrapper InverseFFTPlan;

		std::vector<ComponentWorkStore> Component;
		TDEmResponse WR; // Work response class

		double FrequencyLog10Spacing = 0.0;
		double DiscreteFrequencyLow = 0.0;
		double DiscreteFrequencyHigh = 0.0;
		size_t NumberOfDiscreteFrequencies = 0;
		std::vector<double> DiscreteFrequencies;
		std::vector<double> DiscreteFrequenciesLog10;

		size_t NumberOfSplinedFrequencies = 0;
		std::vector<double> SplinedFrequencieslog10;

		std::vector<LowPassFilter> Filters;

		WindowingScheme WindScheme;
		Waveform WvForm;
		Transmitter Tx;

		Vec3d Scale;
		Vec3d RefGeomPrimary; // For PPM systems
		Mat3d RotMatrixToRxFrame;

	public:

		const WindowSpecification& window(const size_t w) const { return WindScheme.Windows[w]; }
		const Waveform& waveform() const { return WvForm; }
		const Transmitter& transmitter() const { return Tx; }

		ModellingOptions MO;

		TDEmSystem(const fs::path& descriptorpath) {
			read_system_descriptor_file(descriptorpath);
		};

		TDEmSystem() {};
		
		const size_t& nWindows() const {
			return WindScheme.nWindows();
		}
				
		double primary(const size_t component) const {
			assert(component < NCOMP);
			return WR.P[component][0];
		}

		double secondary(const size_t component, const size_t window) {
			assert(component < NCOMP);
			return WR.S[component][window];
		}

		void read_system_descriptor_file(const fs::path& systemdescriptorfile) {
			if (!fs::exists(systemdescriptorfile)) {
				std::string msg = strprint("\n\tD'Oh! the specified system descriptor file (%s) does not exist\n", systemdescriptorfile.string().c_str());
				glog.errormsg(_SRC_, msg);
			}

			STM = cBlock(systemdescriptorfile);
			SystemName = STM.getstringvalue("Name");
			SystemType = STM.getstringvalue("Type");

			if (strcasecmp(SystemType, "Time Domain") != 0) {
				glog.errormsg(_SRC_, "System Type is not Time Domain\n");
			}

			cBlock b = STM.findblock("Transmitter");
			Tx.NumberOfTurns = b.getdoublevalue("NumberOfTurns");
			Tx.PeakCurrent = b.getdoublevalue("PeakCurrent");
			Tx.LoopArea = b.getdoublevalue("LoopArea");

			WvForm.initialise(b, systemdescriptorfile);

			cBlock rxblock = STM.findblock("Receiver");
			WindScheme = WindowingScheme(rxblock, WvForm);
			setup_nwindows();

			if (WvForm.Time.size() <= 2 || WvForm.Time.size() != WvForm.TD_Waveform.size()) {
				glog.errormsg(_SRC_, "The number of WaveformTime values must match number of WaveformCurrent/WaveformReceived values and also be more than two\n");
			}

			//Load low pass filters
			auto v1 = STM.getdoublevector("Receiver.LowPassFilter.Order");
			auto v2 = STM.getdoublevector("Receiver.LowPassFilter.CutOffFrequency");
			if (v1.size() != v2.size()) {
				glog.errormsg(_SRC_, "Filter CutOffFrequency and Order sizes must be equal.");
			}
			for (size_t i = 0; i < v1.size(); i++) {
				Filters.push_back(LowPassFilter(v1[i], v2[i]));
			}

			b = STM.findblock("ForwardModelling");
			MO = ModellingOptions(b);

			setup_discrete_frequencies();
			setup_transforms();
			setup_splines();
			setup_scaling();
		};

		// Modelling

		const TDEmVectorResponse& forward_model_primary_field(const TDEmGeometry& G) {
			set_geometry(G);
			set_calculationtype(CMode::FM);
			setprimaryfields();
			return WR.P;
		};

		const TDEmResponse& forward_model(const Earth1D& E, const TDEmGeometry& G) {
			set_earth(E);
			set_geometry(G);
			setup_computations();
			set_calculationtype(CMode::FM);
			setprimaryfields();
			setsecondaryfields();
			return WR;
		};

		const TDEmResponse& derivative(const CalculationType& calc) {
			set_calculationtype(calc);
			setprimaryfields();
			setsecondaryfields();
			return WR;
		}

		void set_earth(const Earth1D& earth) {
			lem().set_earth(earth);
		};

		void set_geometry(const TDEmGeometry& G) {
			// Set geometry inside the LE Modeller
			Vec3d tx_reference_orientation = Vec3d::UnitZ();
			const Vec3d sep = G.txrx_separation();
			const double& h = G.tx_height;
			const double& x = sep.x();
			const double& y = sep.y();
			const double& z = h + sep.z();
			const Vec3d tx_orientation = G.tx_orientation(Tx.Reference_Orientation);
			lem().set_geometry(tx_orientation, h, x, y, z);
			// Set the rotation matrix for rotating vector fields to Rx frame of reference
			RotMatrixToRxFrame = G.inertial_to_rx_frame_rotation_matrix();
		};

		void set_calculationtype(const CalculationType& _calculationtype) {
			LEM.set_calculationtype(_calculationtype);
		};

		void setup_computations() {
			lem().setup_computations();
		}

		void drx_pitch(double xb, double zb, double p, double& dxbdp, double& dzbdp) {
			//xi = (  xb*cosp  + zb*sinp);Inertial
			//zi = ( -xb*sinp  + zb*cosp);
			//xb = (  xi*cosp  - zi*sinp);As bird sees it
			//zb = (  xi*sinp  + zi*cosp);						

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM || 
				MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				//Must work with true field vector directions, not the PPM scaled versinn
				xb *= RefGeomPrimary[XCOMP];
				zb *= RefGeomPrimary[ZCOMP];
			}

			double cosp = cos(D2R<double> *p);
			double sinp = sin(D2R<double> *p);

			double xi = (xb * cosp + zb * sinp);//convert back to real coordinate system
			double zi = (-xb * sinp + zb * cosp);

			dxbdp = D2R<double> *(-xi * sinp - zi * cosp);
			dzbdp = D2R<double> *(+xi * cosp - zi * sinp);

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM || 
				MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dxbdp /= RefGeomPrimary[XCOMP];
				dzbdp /= RefGeomPrimary[ZCOMP];
			}
		}

		void drx_pitch(std::vector<double> xb, std::vector<double> zb, double p, std::vector<double>& dxbdp, std::vector<double>& dzbdp) {
			//xi = (  xb*cosp  + zb*sinp);Inertial
			//zi = ( -xb*sinp  + zb*cosp);
			//xb = (  xi*cosp  - zi*sinp);As bird sees it
			//zb = (  xi*sinp  + zi*cosp);						

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM ||
				MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				//Must work with true field vector directions, not the PPM scaled versinn
				xb *= RefGeomPrimary[XCOMP];
				zb *= RefGeomPrimary[ZCOMP];
			}


			double cosp = cos(D2R<double> *p);
			double sinp = sin(D2R<double> *p);

			//convert back to real coordinate system
			std::vector<double> xi = (xb * cosp + zb * sinp);
			std::vector<double> zi = (xb * -sinp + zb * cosp);

			dxbdp = (xi * -sinp - zi * cosp) * D2R<double>;
			dzbdp = (xi * cosp - zi * sinp) * D2R<double>;

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM ||
				MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dxbdp /= RefGeomPrimary[XCOMP];
				dzbdp /= RefGeomPrimary[ZCOMP];
			}
		}

		void drx_roll(double yb, double zb, double r, double& dybdr, double& dzbdr) {
			//yi = (  yb*cosr  - zb*sinr);Inertial
			//zi = (  yb*sinr  + zb*cosr);
			//yb = (  yi*cosr  + zi*sinr);As bird sees it
			//zb = ( -yi*sinr  + zi*cosr);						

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM ||
				MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				//Must work with true field vector directions, not the PPM scaled versinn
				yb *= RefGeomPrimary[YCOMP];
				zb *= RefGeomPrimary[ZCOMP];
			}

			double cosr = cos(D2R<double> *r);
			double sinr = sin(D2R<double> *r);

			double yi = (yb * cosr - zb * sinr);//convert back to real coordinate system
			double zi = (yb * sinr + zb * cosr);

			dybdr = D2R<double> *(-yi * sinr + zi * cosr);
			dzbdr = D2R<double> *(-yi * cosr - zi * sinr);

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM ||
				MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dybdr /= RefGeomPrimary[YCOMP];
				dzbdr /= RefGeomPrimary[ZCOMP];
			}
		}

		void  drx_roll(std::vector<double> yb, std::vector<double> zb, double r, std::vector<double>& dybdr, std::vector<double>& dzbdr) {
			//yi = (  yb*cosr  - zb*sinr);Inertial
			//zi = (  yb*sinr  + zb*cosr);
			//yb = (  yi*cosr  + zi*sinr);As bird sees it
			//zb = ( -yi*sinr  + zi*cosr);						

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM ||
				MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				//Must work with true field vector directions, not the PPM scaled versinn
				yb *= RefGeomPrimary[YCOMP];
				zb *= RefGeomPrimary[ZCOMP];
			}

			double cosr = cos(D2R<double> *r);
			double sinr = sin(D2R<double> *r);

			//convert back to real coordinate system
			std::vector<double> yi = (yb * cosr - zb * sinr);
			std::vector<double> zi = (yb * sinr + zb * cosr);

			dybdr = (yi * -sinr + zi * cosr) * D2R<double>;
			dzbdr = (yi * -cosr - zi * sinr) * D2R<double>;

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM ||
				MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dybdr /= RefGeomPrimary[YCOMP];
				dzbdr /= RefGeomPrimary[ZCOMP];
			}
		}

		void drx_roll_new(const TDEmGeometry& g, const TDEmVectorResponse& fields, TDEmVectorResponse& derivatives) const {
			Mat3d dM = g.rx_roll_derivative_matrix();
			apply_rx_derivative_matrix(dM, fields, derivatives);
		};

		void drx_pitch_new(const TDEmGeometry& g, const TDEmVectorResponse& fields, TDEmVectorResponse& derivatives) const {
			Mat3d dM = g.rx_pitch_derivative_matrix();
			apply_rx_derivative_matrix(dM, fields, derivatives);
		};

		void drx_yaw_new(const TDEmGeometry& g, const TDEmVectorResponse& fields, TDEmVectorResponse& derivatives) const {
			Mat3d dM = g.rx_yaw_derivative_matrix();
			apply_rx_derivative_matrix(dM, fields, derivatives);
		};

	private:

		void setprimaryfields() {
			Vec3d v = lem().primaryfield_inertial();
			//std::cout << v << std::endl;

			// Rotate field to Rx frame
			v = RotMatrixToRxFrame * v;

			if (lem().cmode() == CMode::DH) {
				//This is because when H changes Z also changes and DZ == DH ... //but they should be all zero anyway
				v *= 2.0;
			}

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				v *= 2.0;
			}

			if (MO.OutputType == ModellingOptions::OutputType::DBDT) {
				//Must convert to dB/dt. This happens implicitly for the secondary via the waveform.
				v *= Tx.PeakdIdT;
			}

			WR.P[XCOMP][0] = v.x() * Scale[XCOMP];
			WR.P[YCOMP][0] = v.y() * Scale[YCOMP];
			WR.P[ZCOMP][0] = v.z() * Scale[ZCOMP];
		};

		void setsecondaryfields() {
			//Computation for discrete frequencies 	
			for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++) {
				Vec3cd v = lem().secondaryfield_inertial(fi);
				//std::cout << v << std::endl;

				// Rotate field to Rx frame
				v = RotMatrixToRxFrame * v;

				//std::cout << RotMatrixToRxFrame << std::endl;

				if (lem().cmode() == CMode::DH) {
					//This is because when H changes Z also changes and DZ == DH
					v *= 2.0;
				}

				for (size_t ci = 0; ci < NCOMP; ci++) {
					Component[ci].IR_discrete_real[fi] = v[ci].real();
					Component[ci].IR_discrete_imag[fi] = v[ci].imag();
				}
			};

			//Spline discreet frequencies
			for (size_t ci = 0; ci < NCOMP; ci++) {
				if (Scale[ci] == 0.0) return;
				spline_component(ci);
				inverse_fft_window_scale_component(ci);
			}

			if (MO.SaveDiagnosticFiles) {
				write_discretefrequencies("diag_discretefrequencies.txt");
				write_splinedfrequencies("diag_splinedfrequencies.txt");
				WvForm.write_frequencydomainwaveform("diag_frequencydomainwaveform.txt");
			}

			if (MO.SaveDiagnosticFiles) {
				WindScheme.write_windows("diag_windows.txt", WR.secondary(XCOMP), WR.secondary(YCOMP), WR.secondary(ZCOMP));
			}
		}

		void setup_nwindows() {
			WR.resize(nWindows());
		}

		void apply_rx_derivative_matrix(const Mat3d& dM, const TDEmVectorResponse& fields, TDEmVectorResponse& derivatives) const {
			const size_t n = fields.size();
			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM || MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				for (size_t i = 0; i < n; i++) {
					Vec3d ftrue = fields.get_vec3d(i);
					//Must work with true field vector directions (not the PPM scaled versinn)
					ftrue[XCOMP] *= RefGeomPrimary[XCOMP];
					ftrue[YCOMP] *= RefGeomPrimary[YCOMP];
					ftrue[ZCOMP] *= RefGeomPrimary[ZCOMP];
					derivatives.set_vec3d(i, dM * ftrue);
					//Convert back to PPMS
					derivatives[XCOMP][i] /= RefGeomPrimary[XCOMP];
					derivatives[YCOMP][i] /= RefGeomPrimary[YCOMP];
					derivatives[ZCOMP][i] /= RefGeomPrimary[ZCOMP];
				}
			}
			else {
				for (size_t i = 0; i < n; i++) {
					derivatives.set_vec3d(i,dM * fields.get_vec3d(i));
				}
			}
		}

		void setup_transforms() {
			WvForm.NumFrequencies = WvForm.NumSamples / 2 + 1;
			WvForm.FFT_Frequency.resize(WvForm.NumFrequencies);
			for (size_t k = 0; k < WvForm.NumFrequencies; k++) {
				WvForm.FFT_Frequency[k] = WvForm.calculate_fft_frequency(k);
			}

			int N = WvForm.NumSamples;
			size_t NComplex = N;
			size_t NReal = 2 * (N / 2 + 1);

			// Forward transform
			WvForm.FD_Waveform.resize(NComplex);
			WvForm.FFT_WorkArray.resize(NReal);//Inverse transform work array	
			WvForm.TransferFunction.resize(WvForm.NumFrequencies);

			FFTWPlanWrapper ForwardFFTPlan(fftw_plan_dft_r2c_1d(N, (double*)WvForm.TD_Waveform.data(), (fftw_complex*)WvForm.FD_Waveform.data(), FFTW_ESTIMATE));
			ForwardFFTPlan.execute();
			const double scale = 1.0 / (double)WvForm.NumSamples;
			WvForm.FD_Waveform *= scale; // Scale the spectrum
			
			bool convert_B_2_dBdT = false;
			bool convert_dBdT_2_B = false;
			if (WvForm.Type == Waveform::Type::TX) {
				if (MO.OutputType == ModellingOptions::OutputType::DBDT) {
					convert_B_2_dBdT = true;
				}
			}

			if (WvForm.Type == Waveform::Type::RX) {
				if (MO.OutputType == ModellingOptions::OutputType::BFIELD) {
					convert_dBdT_2_B = true;
				}
			}

			WvForm.TransferFunction = WvForm.FD_Waveform;
			for (size_t k = 0; k < WvForm.NumFrequencies; k++) {
				const double& frequency = WvForm.FFT_Frequency[k];
				if (convert_B_2_dBdT == true) {
					WvForm.TransferFunction[k] *= cdouble(0.0, -TWOPI<double> *frequency);
				}
				if (convert_dBdT_2_B == true) {
					WvForm.TransferFunction[k] *= cdouble(0.0, -1.0 / (TWOPI<double> *frequency));
				}

				for (size_t fi = 0; fi < Filters.size(); fi++) {
					const cdouble w = Filters[fi].weight(frequency);
					WvForm.TransferFunction[k] *= w;
				}
			}

			// FM to be splined
			NumberOfSplinedFrequencies = WvForm.NumFrequencies / 2;
			SplinedFrequencieslog10.resize(NumberOfSplinedFrequencies);
			for (size_t k = 0; k < NumberOfSplinedFrequencies; k++) {
				SplinedFrequencieslog10[k] = log10(fabs(WvForm.FFT_Frequency[k * 2 + 1]));
			}

			// FFTW_MEASURE does not seem to be thread safe
			#if defined MULTITHREADED
				unsigned int FFTW_FLAGS = FFTW_ESTIMATE;
			#else
				unsigned int FFTW_FLAGS = FFTW_MEASURE;
			#endif

			InverseFFTPlan.setplan(fftw_plan_dft_c2r_1d(N, (fftw_complex*)WvForm.FFT_WorkArray.data(), (double*)WvForm.FFT_WorkArray.data(), FFTW_FLAGS));
		}

		void setup_scaling() {
			Tx.PeakdIdT = WvForm.compute_peak_didt();
			double tx_scale = MUZERO<double> * Tx.LoopArea * Tx.NumberOfTurns * Tx.PeakCurrent;
			
			//ModellingOptions
			Scale[XCOMP] = tx_scale * MO.XOutputScaling;
			Scale[YCOMP] = tx_scale * MO.YOutputScaling;
			Scale[ZCOMP] = tx_scale * MO.ZOutputScaling;

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM || MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				cBlock b = STM.findblock("ReferenceGeometry");
				if (b.Entries.size() == 0) {
					glog.errormsg(_SRC_, "Must define a ReferenceGeometry for PPM or PPMPEAKTOPEAK normalisation\n");
				}

				//Todo check this is working okay
				TDEmGeometry NormalizationGeometry(b);
				set_geometry(NormalizationGeometry);
				setprimaryfields();

				double s = 1.0;
				if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM) {
					s *= 1.0e6;
				}
				else if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
					s *= 1.0e6;
				}

				for (size_t ci = 0; ci < NCOMP; ci++) {
					RefGeomPrimary[ci] = WR.P[ci][0];
					if (RefGeomPrimary[ci] == 0.0) Scale[ci] = 0.0;
					else Scale[ci] *= (s / RefGeomPrimary[ci]);
				}
			}
		}

		void setup_discrete_frequencies() {
			double lf1 = log10(WvForm.BaseFrequency);
			double lf2 = log10(WvForm.SampleFrequency / 2);
			double dlf = 1.0 / MO.FrequenciesPerDecade;

			lf1 = lf1 - 2.0 * dlf;
			lf2 = lf2 + 2.0 * dlf;

			size_t nf = (size_t)ceil((lf2 - lf1) * MO.FrequenciesPerDecade);
			dlf = (lf2 - lf1) / double(nf - 1);

			NumberOfDiscreteFrequencies = nf;
			FrequencyLog10Spacing = dlf;
			DiscreteFrequencyLow = std::pow(10.0, lf1);
			DiscreteFrequencyHigh = std::pow(10.0, lf2);

			DiscreteFrequenciesLog10 = std::vector<double>(NumberOfDiscreteFrequencies);
			DiscreteFrequencies = std::vector<double>(NumberOfDiscreteFrequencies);
			for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++) {
				DiscreteFrequenciesLog10[fi] = log10(DiscreteFrequencyLow) + FrequencyLog10Spacing * (double)fi;
				DiscreteFrequencies[fi] = pow(10.0, DiscreteFrequenciesLog10[fi]);
			}
			lem().initialise(DiscreteFrequencies, MO.NumAbscissa, MO.ModellingLoopRadius);
		}

		void setup_splines() {
			Component.resize(NCOMP);
			Component[XCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies, nWindows());
			Component[YCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies, nWindows());
			Component[ZCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies, nWindows());
			FrequencySpliner.initialise(DiscreteFrequenciesLog10, SplinedFrequencieslog10);
			//lem().initialise_frequencies(DiscreteFrequencies);
		};
		
		void spline_component(const size_t& component) {
			ComponentWorkStore& C = Component[component];
			const std::vector<double>& v = FrequencySpliner.interpolated_values();

			const size_t n = v.size();
			double* a = (double*)(C.IR_splined.data());

			FrequencySpliner.compute_interpolation(C.IR_discrete_real);
			for (size_t i = 0; i < n; i++) {
				a[i * 2] = v[i];
			}

			FrequencySpliner.compute_interpolation(C.IR_discrete_imag);
			for (size_t i = 0; i < n; i++) {
				a[i * 2 + 1] = v[i];
			}
		};

		void inverse_fft_window_scale_component(const size_t& component) {
			ComponentWorkStore& C = Component[component];
			// Reset to the stored transfer function
			WvForm.FFT_WorkArray = WvForm.TransferFunction;
			// Apply transfer function
			size_t n = 0;
			for (size_t k = 1; k < WvForm.NumFrequencies; k += 2) {
				WvForm.FFT_WorkArray[k] *= C.IR_splined[n];
				n++;
			}

			// Inverse FFT

			InverseFFTPlan.execute();

			// Window
			WindScheme.computewindow((double*)WvForm.FFT_WorkArray.data(), WR.S[component]);
			if (MO.SaveDiagnosticFiles) {
				write_timesseries("diag_xtimeseries.txt");
			}
			// Scale
			WR.S[component] *= Scale[component];
		}

		void write_discretefrequencies(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumberOfDiscreteFrequencies; i++) {
				ofs << strprint("%15le\t%15le\t%15le\t%15le\t%15le\t%15le\t%15le\n", DiscreteFrequencies[i],
					Component[XCOMP].IR_discrete_real[i],
					Component[XCOMP].IR_discrete_imag[i],
					Component[YCOMP].IR_discrete_real[i],
					Component[YCOMP].IR_discrete_imag[i],
					Component[ZCOMP].IR_discrete_real[i],
					Component[ZCOMP].IR_discrete_imag[i]);
			}
		}

		void write_splinedfrequencies(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumberOfSplinedFrequencies; i++) {
				double f = pow10(SplinedFrequencieslog10[i]);
				ofs << strprint("%15le\t%15le\t%15le\t%15le\t%15le\t%15le\t%15le\n",
					f,
					Component[XCOMP].IR_splined[i].real(),
					Component[XCOMP].IR_splined[i].imag(),
					Component[YCOMP].IR_splined[i].real(),
					Component[YCOMP].IR_splined[i].imag(),
					Component[ZCOMP].IR_splined[i].real(),
					Component[ZCOMP].IR_splined[i].imag());
			}
		}

		void write_timesseries(const std::string& path) const {
			std::ofstream ofs = ofstream_ex(path);
			double* ts = (double*)(WvForm.FFT_WorkArray.data());
			for (size_t i = 0; i < WvForm.NumSamples; i++) {
				ofs << strprint("%20.10le\t%20.10le\n", WvForm.Time[i], ts[i]);
			}
		}
	};
};


