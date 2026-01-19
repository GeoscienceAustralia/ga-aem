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
#include <iostream>

#include <fftw3.h>
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

	class TDEmSystem : public AEMSystem<double> {

	private:

		using Response = TDEmResponse<double>;
		using VectorResponse = TDEmVectorResponse<double>;
		using ScalarResponse = TDEmScalarResponse<double>;


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

		FixedPointSpline<double> FrequencySpliner;
		FFTWPlanWrapper InverseFFTPlan;

		std::vector<ComponentWorkStore> Component;

		double FrequencyLog10Spacing = 0.0;
		double DiscreteFrequencyLow = 0.0;
		double DiscreteFrequencyHigh = 0.0;
		size_t NumberOfDiscreteFrequencies = 0;
		std::vector<double> DiscreteFrequencies;
		std::vector<double> DiscreteFrequenciesLog10;

		size_t NumberOfSplinedFrequencies = 0;
		std::vector<double> SplinedFrequencieslog10;

		std::vector<LowPassFilter> Filters;
		Waveform WvForm;

	public:

		const WindowSpecification& window(const size_t w) const { return WindScheme.Windows[w]; }
		const Waveform& waveform() const { return WvForm; }
		const Transmitter& transmitter() const { return Tx; }

		ModellingOptions MO;

		TDEmSystem(const std::filesystem::path& descriptorpath) {
			read_system_descriptor_file(descriptorpath);
		};

		TDEmSystem() {};

		static AEM::SystemType get_type() { return AEM::SystemType::TimeDomain; };

		static std::unique_ptr<AEMSystem<double>> unique_ptr(const std::filesystem::path stmfile) {
			return std::make_unique<TDEmSystem>(stmfile);
		};


		SystemType type() const {
			return AEM::SystemType::TimeDomain;
		};

		std::string type_string() const {
			return "TimeDomain";
		};

		void read_system_descriptor_file(const std::filesystem::path& systemdescriptorfile) {
			if (!std::filesystem::exists(systemdescriptorfile)) {
				std::string msg = strprint("\n\tD'Oh! the specified system descriptor file (%s) does not exist\n", systemdescriptorfile.string().c_str());
				glog.errormsg(_SRC_, msg);
			}

			STM = cBlock(systemdescriptorfile);
			SystemName = STM.getstringvalue("Name");

			std::string typestr;
			if (STM.getvalue("Type", typestr)) {
				AEM::SystemType systype = aem_system_type(typestr);
				if (systype != SystemType::TimeDomain) {
					glog.errormsg(_SRC_, "System Type must be 'Time Domain'.\n");
				}
			}
			else glog.errormsg(_SRC_,"The AEM System 'Type' is not specified.\n");
			
			

			cBlock txblock = STM.findblock("Transmitter");
			Tx = Transmitter(txblock);
			
			WvForm.initialise(txblock, systemdescriptorfile);

			cBlock rxblock = STM.findblock("Receiver");
			WindScheme = WindowingScheme(rxblock, WvForm.Time);
			//set_nwindows();

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

			MO = ModellingOptions(STM.findblock("ForwardModelling"));

			setup_discrete_frequencies();
			setup_transforms();
			setup_splines();
			setup_scaling();
		};

	private:

		bool is_ppm_system() const {
			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM || MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				return true;
			}
			return false;
		};

		void set_primaryfields(VectorResponse& P, const Vec3d& txvec, const Mat3d& rxmat) {
			Vec3d v = lem().primaryfield_inertial(txvec);

			// Rotate field to Rx frame
			v = rxmat * v;

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				v *= 2.0;
			}

			if (MO.OutputType == ModellingOptions::OutputType::DBDT) {
				//Must convert to dB/dt. This happens implicitly for the secondary via the waveform.
				v *= Tx.PeakdIdT;
			}

			for (size_t ci = 0; ci < NCOMP; ci++) {
				v[ci] *= Scale[ci];
				std::fill(P[ci].begin(), P[ci].end(), v[ci]);
			}
		};

		void set_secondaryfields(VectorResponse& S, const Vec3d& txvec, const Mat3d& rxmat) {
			//Computation for discrete frequencies 	
			for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++) {
				Vec3cd v = lem().secondaryfield_inertial(fi,txvec);
				// Rotate field to Rx frame
				v = rxmat * v;

				for (size_t ci = 0; ci < NCOMP; ci++) {
					Component[ci].IR_discrete_real[fi] = v[ci].real();
					Component[ci].IR_discrete_imag[fi] = v[ci].imag();
				}
			};

			//Spline discreet frequencies
			for (size_t ci = 0; ci < NCOMP; ci++) {
				if (Scale[ci] == 0.0) return;
				spline_component(ci);
				inverse_fft_window_scale_component(S,ci);
			}

			if (MO.SaveDiagnosticFiles) {
				write_discretefrequencies("diag_discretefrequencies.txt");
				write_splinedfrequencies("diag_splinedfrequencies.txt");
				WvForm.write_frequencydomainwaveform("diag_frequencydomainwaveform.txt");
			}

			if (MO.SaveDiagnosticFiles) {
				WindScheme.write_windows<double,std::vector>("diag_windows.txt", S[XCOMP], S[YCOMP], S[ZCOMP]);
			}
		}

		void setup_transforms() {
			WvForm.NumFrequencies = WvForm.NumSamples / 2 + 1;
			WvForm.FFT_Frequency.resize(WvForm.NumFrequencies);
			for (size_t k = 0; k < WvForm.NumFrequencies; k++) {
				WvForm.FFT_Frequency[k] = WvForm.calculate_fft_frequency(k);
			}

			const size_t N = WvForm.NumSamples;
			const size_t NComplex = N;
			const size_t NReal = 2 * (N / 2 + 1);

			// Forward transform
			WvForm.FD_Waveform.resize(NComplex);
			WvForm.FFT_WorkArray.resize(NReal);//Inverse transform work array	
			WvForm.TransferFunction.resize(WvForm.NumFrequencies);

			FFTWPlanWrapper ForwardFFTPlan(fftw_plan_dft_r2c_1d((int)N, (double*)WvForm.TD_Waveform.data(), (fftw_complex*)WvForm.FD_Waveform.data(), FFTW_ESTIMATE));
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

			InverseFFTPlan.setplan(fftw_plan_dft_c2r_1d((int)N, (fftw_complex*)WvForm.FFT_WorkArray.data(), (double*)WvForm.FFT_WorkArray.data(), FFTW_FLAGS));
		};

		void setup_scaling() {
			Tx.PeakdIdT = WvForm.compute_peak_didt();
			double tx_scale = MUZERO<double> * Tx.LoopArea * Tx.nTurns * Tx.PeakCurrent;

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
				VectorResponse P(nWindows());
				set_primaryfields(P, Tx.Orientation,InertialToRxFrame);

				double s = 1.0;
				if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM) {
					s *= 1.0e6;
				}
				else if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
					s *= 1.0e6;
				}

				for (size_t ci = 0; ci < NCOMP; ci++) {
					RefGeomPrimary[ci] = P[ci][0];
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

		void inverse_fft_window_scale_component(VectorResponse& S, const size_t& component) {
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
			WindScheme.computewindow((double*)WvForm.FFT_WorkArray.data(), S[component]);
			if (MO.SaveDiagnosticFiles) {
				write_timesseries("diag_xtimeseries.txt");
			}
			// Scale
			S[component] *= Scale[component];
		};

		void write_discretefrequencies(const std::filesystem::path& path) const {
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

		void write_splinedfrequencies(const std::filesystem::path& path) const {
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


