/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Brodie Geophysics
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
		size_t NumAbscissa = 21;
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
			if (FrequenciesPerDecade < 6) {
				glog.warningmsg(_SRC_, "It is wise to use at least 6 frequencies per decade\n");
			}

			NumAbscissa = b.getsizetvalue("NumberOfAbsiccaInHankelTransformEvaluation");
			if (NumAbscissa < 21) {
				glog.warningmsg(_SRC_, "It is wise to use at least 21 Absicca for integrating the Hankel Transforms");
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

	class SpectralAEMSystem : public AEMSystem<cdouble> {

	private:

		using Response = TDEmResponse<cdouble>;
		using VectorResponse = TDEmVectorResponse<cdouble>;
		using ScalarResponse = TDEmScalarResponse<cdouble>;

		std::vector<ComponentWorkStore> Component;
		
		// Discrete frequency knots to be splined
		size_t NumberOfKnots = 0;
		double KnotLog10Spacing = 0.0;
		double KnotLow = 0.0;
		double KnotHigh = 0.0;
		std::vector<double> Knots;
		std::vector<double> KnotsLog10;
		FixedPointSpline<double> FrequencySpliner;

		// Odd harmonic frequencies where splines are to be evaluated
		std::vector<double> FrequencySeries;
		std::vector<double> FrequencySeriesLog10;

	public:

		const WindowSpecification& window(const size_t w) const { return WindScheme.Windows[w]; }
		const Transmitter& transmitter() const { return Tx; }

		ModellingOptions MO;

		SpectralAEMSystem(const fs::path& descriptorpath) {
			read_system_descriptor_file(descriptorpath);
		};

		SpectralAEMSystem() {};

		static AEM::SystemType get_type() { return AEM::SystemType::SpectralTimeDomain; };
		
		static std::unique_ptr<AEMSystem<cdouble>> unique_ptr(const fs::path stmfile) {
			return std::make_unique<SpectralAEMSystem>(stmfile);
		};

		SystemType type() const {
			return AEM::SystemType::SpectralTimeDomain;
		};

		std::string type_string() const {
			return "SpectralTimeDomain";
		};

		void read_system_descriptor_file(const fs::path& systemdescriptorfile) {
			if (!fs::exists(systemdescriptorfile)) {
				std::string msg = strprint("\n\tD'Oh! the specified system descriptor file (%s) does not exist\n", systemdescriptorfile.string().c_str());
				glog.errormsg(_SRC_, msg);
			}

			STM = cBlock(systemdescriptorfile);
			SystemName = STM.getstringvalue("Name");

			std::string typestr;
			if (STM.getvalue("Type", typestr)) {
				AEM::SystemType systype = aem_system_type(typestr);
				if (systype != SystemType::SpectralTimeDomain) {
					glog.errormsg(_SRC_, "System Type must be 'Spectral Time Domain'.\n");
				}
			}
			else glog.errormsg(_SRC_, "The AEM System 'Type' is not specified.\n");

			cBlock txb = STM.findblock("Transmitter");
			Tx = Transmitter(txb);

			cBlock rxb = STM.findblock("Receiver");
			Rx = Receiver(rxb);

			MO = ModellingOptions(STM.findblock("ForwardModelling"));
			setup_frequencies();
			WindScheme = WindowingScheme(rxb, FrequencySeries);
			//set_nwindows();
			
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

			v[XCOMP] *= Scale[XCOMP];
			v[YCOMP] *= Scale[YCOMP];
			v[ZCOMP] *= Scale[ZCOMP];
			std::fill(P[XCOMP].begin(), P[XCOMP].end(), v[XCOMP]);
			std::fill(P[YCOMP].begin(), P[YCOMP].end(), v[YCOMP]);
			std::fill(P[ZCOMP].begin(), P[ZCOMP].end(), v[ZCOMP]);
		};

		void set_secondaryfields(VectorResponse& S, const Vec3d& txvec, const Mat3d& rxmat) {
			//Computation for discrete frequencies 	
			for (size_t fi = 0; fi < NumberOfKnots; fi++) {
				Vec3cd v = lem().secondaryfield_inertial(fi,txvec);

				//std::cout << v << std::endl;

				// Rotate field to Rx frame
				v = rxmat * -v;
				
				for (size_t ci = 0; ci < NCOMP; ci++) {
					Component[ci].IR_discrete_real[fi] = v[ci].real();
					Component[ci].IR_discrete_imag[fi] = v[ci].imag();
				}
			};

			//Spline discreet frequencies
			for (size_t ci = 0; ci < NCOMP; ci++) {
				if (Scale[ci] == 0.0) return;
				spline_component(ci);
				window_component(S,ci);
				scale_component(S,ci);
			}

			if (MO.SaveDiagnosticFiles) {
				write_discretefrequencies("diag_discretefrequencies.txt");
				write_splinedfrequencies("diag_splinedfrequencies.txt");
			}

			if (MO.SaveDiagnosticFiles) {
				WindScheme.write_windows<cdouble,std::vector>("diag_windows.txt", S[XCOMP], S[YCOMP], S[ZCOMP]);
			}
		}

		void setup_scaling() {
			//Tx.PeakdIdT = WvForm.compute_peak_didt();
			double tx_scale = MUZERO<double> * Tx.LoopArea * Tx.nTurns * Tx.PeakCurrent;
			double rx_scale = Rx.Area * Rx.nTurns;

			//ModellingOptions
			Scale[XCOMP] = tx_scale * rx_scale * MO.XOutputScaling;
			Scale[YCOMP] = tx_scale * rx_scale * MO.YOutputScaling;
			Scale[ZCOMP] = tx_scale * rx_scale * MO.ZOutputScaling;

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM || MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				cBlock b = STM.findblock("ReferenceGeometry");
				if (b.Entries.size() == 0) {
					glog.errormsg(_SRC_, "Must define a ReferenceGeometry for PPM or PPMPEAKTOPEAK normalisation\n");
				}

				//Todo check this is working okay
				TDEmGeometry NormalizationGeometry(b);
				set_geometry(NormalizationGeometry);
				
				VectorResponse P(nWindows());
				set_primaryfields(P,Tx.Orientation,InertialToRxFrame);

				double s = 1.0;
				if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM) {
					s *= 1.0e6;
				}
				else if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
					s *= 1.0e6;
				}

				for (size_t ci = 0; ci < NCOMP; ci++) {
					RefGeomPrimary[ci] = P[ci][0].real();
					if (RefGeomPrimary[ci] == 0.0) Scale[ci] = 0.0;
					else Scale[ci] *= (s / RefGeomPrimary[ci]);
				}
			}
		};

		void setup_frequencies() {
			size_t N = (size_t) ((Rx.SamplingFrequency / Tx.BaseFrequency) / 2.0 / 2.0);
			FrequencySeries = increment(N, Tx.BaseFrequency, Tx.BaseFrequency * 2.0);
			FrequencySeriesLog10 = log10(FrequencySeries);

			double lf1 = log10(Tx.BaseFrequency); // Lowest harmonic
			double lf2 = log10(Rx.SamplingFrequency / 2.0); // Nyquist
			double dlf = 1.0 / MO.FrequenciesPerDecade;

			lf1 = lf1 - 2.0 * dlf;
			lf2 = lf2 + 2.0 * dlf;

			size_t nf = (size_t)ceil((lf2 - lf1) * MO.FrequenciesPerDecade);
			dlf = (lf2 - lf1) / double(nf - 1);

			NumberOfKnots = nf;
			KnotLog10Spacing = dlf;
			KnotLow = std::pow(10.0, lf1);
			KnotHigh = std::pow(10.0, lf2);

			KnotsLog10 = std::vector<double>(NumberOfKnots);
			Knots = std::vector<double>(NumberOfKnots);
			for (size_t fi = 0; fi < NumberOfKnots; fi++) {
				KnotsLog10[fi] = log10(KnotLow) + KnotLog10Spacing * (double)fi;
				Knots[fi] = pow(10.0, KnotsLog10[fi]);
			}
			lem().initialise(Knots, MO.NumAbscissa, MO.ModellingLoopRadius);
		};

		void setup_splines() {
			const size_t ndf = Knots.size();
			const size_t nsf = FrequencySeriesLog10.size();
			Component.resize(NCOMP);
			Component[XCOMP].resize(ndf, nsf, nWindows());
			Component[YCOMP].resize(ndf, nsf, nWindows());
			Component[ZCOMP].resize(ndf, nsf, nWindows());
			FrequencySpliner.initialise(KnotsLog10, FrequencySeriesLog10);
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

		void window_component(VectorResponse& V, const size_t& component) {
			ComponentWorkStore& C = Component[component];
			// Window
			WindScheme.computewindow(C.IR_splined.data(), V[component]);
			//if (MO.SaveDiagnosticFiles) {
			//	write_frequencyseries("diag_xtimeseries.txt");
			//}
		};

		void scale_component(VectorResponse& V, const size_t& component) {
			//ComponentWorkStore& C = Component[component];
			V[component] *= Scale[component];
		};

		void write_discretefrequencies(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumberOfKnots; i++) {
				ofs << strprint("%15le\t%15le\t%15le\t%15le\t%15le\t%15le\t%15le\n", Knots[i],
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
			const size_t nsf = FrequencySeriesLog10.size();
			for (size_t i = 0; i < nsf; i++) {
				double f = pow10(FrequencySeriesLog10[i]);
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

		//void write_timesseries(const std::string& path) const {
		//	std::ofstream ofs = ofstream_ex(path);
		//	double* ts = (double*)(WvForm.FFT_WorkArray.data());
		//	for (size_t i = 0; i < WvForm.NumSamples; i++) {
		//		ofs << strprint("%20.10le\t%20.10le\n", WvForm.Time[i], ts[i]);
		//	}
		//};
	};
};


