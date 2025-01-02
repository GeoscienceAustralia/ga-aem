/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include <cassert>
#include <stdexcept>
#include <complex>

#include "aemsystem.hpp"
#include "fftwplanwrapper.hpp"
#include "vector_utils.hpp"
#include "file_utils.hpp"
#include "general_utils.hpp"
#include "blocklanguage.hpp"
#include "eigen_utils.hpp"
#include "earth1d.hpp"
#include "fixed_point_spline.hpp"
#include "tdemgeometry.hpp"

namespace AEM {

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
		};

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
			set_nwindows();

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
			set_primaryfields();
			return WR.P;
		};

		const TDEmResponse& forward_model(const Earth1D& E, const TDEmGeometry& G) {
			set_earth(E);
			set_geometry(G);
			setup_computations();
			set_calculationtype(CMode::FM);
			set_primaryfields();
			set_secondaryfields();
			return WR;
		};

		TDEmVectorResponse derivative(const CalculationType& calc, const TDEmGeometry& G, const TDEmVectorResponse& forward_model) {
			TDEmVectorResponse derivative(nWindows());
			if (calc.get_mode() == CMode::DRX_ROLL) {
				drx_roll(G, forward_model, derivative);
			}
			else if (calc.get_mode() == CMode::DRX_PITCH) {
				drx_pitch(G, forward_model, derivative);
			}
			else if (calc.get_mode() == CMode::DRX_YAW) {
				drx_yaw(G, forward_model, derivative);
			}
			else {
				glog.errormsg(_SRC_, "Invalid derivative operation.");
			}
			return derivative;
		};

		const TDEmResponse& derivative(const CalculationType& calc, const TDEmGeometry& G, const TDEmResponse& forward_model) {
			TDEmResponse& derivative = WR;
			if (calc.get_mode() == CMode::DRX_ROLL) {
				drx_roll(G, forward_model, derivative);
			}
			else if (calc.get_mode() == CMode::DRX_PITCH) {
				drx_pitch(G, forward_model, derivative);
			}
			else if (calc.get_mode() == CMode::DRX_YAW) {
				drx_yaw(G, forward_model, derivative);
			}
			else{
				glog.errormsg(_SRC_,"Invalid derivative operation.");
			}
			return derivative;
		};

		const TDEmResponse& derivative(const CalculationType& calc) {
			if (calc.get_mode() == CMode::DTX_HEIGHT) {
				// This is because when H changes Z also changes
				set_calculationtype(CMode::DZ);
				set_primaryfields();
				set_secondaryfields();
				TDEmResponse DZ = WR;

				set_calculationtype(CMode::DH);
				set_primaryfields();
				set_secondaryfields();
				WR += DZ;
			}
			else {
				set_calculationtype(calc);
				set_primaryfields();
				set_secondaryfields();
			}
			return WR;
		};

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
		};

	private:

		void drx_roll(const TDEmGeometry& G, const TDEmVectorResponse& forward_model, TDEmVectorResponse& derivatives) const {
			const Mat3d dM = G.rx_roll_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model, derivatives);
		};

		void drx_pitch(const TDEmGeometry& G, const TDEmVectorResponse& forward_model, TDEmVectorResponse& derivatives) const {
			const Mat3d dM = G.rx_pitch_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model, derivatives);
		};

		void drx_yaw(const TDEmGeometry& G, const TDEmVectorResponse& forward_model, TDEmVectorResponse& derivatives) const {
			const Mat3d dM = G.rx_yaw_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model, derivatives);
		};

		void drx_roll(const TDEmGeometry& G, const TDEmResponse& forward_model, TDEmResponse& derivatives) const {
			const Mat3d dM = G.rx_roll_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model.P, derivatives.P);
			apply_rx_derivative_matrix(dM, forward_model.S, derivatives.S);
		};

		void drx_pitch(const TDEmGeometry& G, const TDEmResponse& forward_model, TDEmResponse& derivatives) const {
			const Mat3d dM = G.rx_pitch_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model.P, derivatives.P);
			apply_rx_derivative_matrix(dM, forward_model.S, derivatives.S);
		};

		void drx_yaw(const TDEmGeometry& G, const TDEmResponse& forward_model, TDEmResponse& derivatives) const {
			const Mat3d dM = G.rx_yaw_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model.P, derivatives.P);
			apply_rx_derivative_matrix(dM, forward_model.S, derivatives.S);
		};

		void apply_rx_derivative_matrix(const Mat3d& dM, const TDEmVectorResponse& fields, TDEmVectorResponse& derivatives) const {
			const size_t n = fields.nWindows();
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
				for (size_t wi = 0; wi < n; wi++) {
					derivatives.set_vec3d(wi, dM * fields.get_vec3d(wi));
				}
			}
		};

		void set_primaryfields() {
			Vec3d v = lem().primaryfield_inertial();
			//std::cout << v << std::endl;

			// Rotate field to Rx frame
			v = RotMatrixToRxFrame * v;
			//std::cout << v << std::endl;

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				v *= 2.0;
			}

			if (MO.OutputType == ModellingOptions::OutputType::DBDT) {
				//Must convert to dB/dt. This happens implicitly for the secondary via the waveform.
				v *= Tx.PeakdIdT;
			}

			v[XCOMP] *= Scale[XCOMP];
			v[YCOMP] *= Scale[YCOMP];
			v[ZCOMP] *= Scale[ZCOMP];

			//WR.P[XCOMP][0] = v.x() * Scale[XCOMP];
			//WR.P[YCOMP][0] = v.y() * Scale[YCOMP];
			//WR.P[ZCOMP][0] = v.z() * Scale[ZCOMP];
			//for (size_t wi = 1; wi < nWindows(); wi++) {
			//	WR.P[XCOMP][wi] = WR.P[XCOMP][0];
			//	WR.P[YCOMP][wi] = WR.P[XCOMP][0];
			//	WR.P[ZCOMP][wi] = WR.P[XCOMP][0];
			//}

			std::fill(WR.P[XCOMP].begin(), WR.P[XCOMP].end(), v[XCOMP]);
			std::fill(WR.P[YCOMP].begin(), WR.P[YCOMP].end(), v[YCOMP]);
			std::fill(WR.P[ZCOMP].begin(), WR.P[ZCOMP].end(), v[ZCOMP]);
		};

		void set_secondaryfields() {
			//Computation for discrete frequencies 	
			for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++) {
				Vec3cd v = lem().secondaryfield_inertial(fi);
				//std::cout << v << std::endl;

				// Rotate field to Rx frame
				v = RotMatrixToRxFrame * v;

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

		void set_nwindows() {
			WR.set_nWindows(nWindows());
		};

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
		};

		void setup_scaling() {
			Tx.PeakdIdT = WvForm.compute_peak_didt();
			double tx_scale = MUZERO<double> *Tx.LoopArea * Tx.NumberOfTurns * Tx.PeakCurrent;

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
				set_primaryfields();

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
		};

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
		};

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
		};

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
		};

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
		};

		void write_timesseries(const std::string& path) const {
			std::ofstream ofs = ofstream_ex(path);
			double* ts = (double*)(WvForm.FFT_WorkArray.data());
			for (size_t i = 0; i < WvForm.NumSamples; i++) {
				ofs << strprint("%20.10le\t%20.10le\n", WvForm.Time[i], ts[i]);
			}
		};
	};
};


