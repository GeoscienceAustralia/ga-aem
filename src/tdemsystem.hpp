/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include <stdexcept>
#include <complex>
#include "fftw3.h"
#include "file_utils.hpp"
#include "vector_utils.hpp"
#include "general_utils.hpp"
#include "geometry3d.hpp"
#include "blocklanguage.hpp"
#include "earth1d.hpp"
#include "lem.hpp"
#include "fixed_point_spline.hpp"

namespace AEM {

	struct sTDEmNoiseModelComponent {
		double MultiplicativeNoise;
		std::vector<double> AdditiveNoise;
	};

	struct sTDEmNoiseModel {
		sTDEmNoiseModelComponent xcomponent;
		sTDEmNoiseModelComponent ycomponent;
		sTDEmNoiseModelComponent zcomponent;
	};

	struct cTDEmComponent {
		double Primary = 0.0;
		std::vector<double> Secondary;
	};

	class cTDEmData {

	private:
		cTDEmComponent data[3];

	public:

		cTDEmData() {};

		cTDEmComponent& component(const size_t i) { return data[i]; }
		cTDEmComponent& xcomponent() { return data[0]; }
		cTDEmComponent& ycomponent() { return data[1]; }
		cTDEmComponent& zcomponent() { return data[2]; }

	};

	class cTDEmResponse {

	public:
		double PX;
		double PY;
		double PZ;
		std::vector<double> SX;
		std::vector<double> SY;
		std::vector<double> SZ;
	};

	class cTDEmGeometry {

	public:
		enum class ElementType {
			tx_height,
			tx_roll, tx_pitch, tx_yaw,
			txrx_dx, txrx_dy, txrx_dz,
			rx_roll, rx_pitch, rx_yaw,
			unknown
		};


		double tx_height = 0.0;
		double tx_roll = 0.0;
		double tx_pitch = 0.0;
		double tx_yaw = 0.0;
		double txrx_dx = 0.0;
		double txrx_dy = 0.0;
		double txrx_dz = 0.0;
		double rx_roll = 0.0;
		double rx_pitch = 0.0;
		double rx_yaw = 0.0;

		cTDEmGeometry() {};

		cTDEmGeometry(const double& _tx_height, const double& _tx_roll, const double& _tx_pitch, const double& _tx_yaw, const double& _txrx_dx, const double& _txrx_dy, const double& _txrx_dz, const double& _rx_roll, const double& _rx_pitch, const double& _rx_yaw) {
			tx_height = _tx_height;
			tx_roll = _tx_roll; tx_pitch = _tx_pitch; tx_yaw = _tx_yaw;
			txrx_dx = _txrx_dx; txrx_dy = _txrx_dy; txrx_dz = _txrx_dz;
			rx_roll = _rx_roll; rx_pitch = _rx_pitch; rx_yaw = _rx_yaw;
		}

		cTDEmGeometry(const std::vector<double> gvector) {
			for (size_t i = 0; i < size(); i++) {
				(*this)[i] = gvector[i];
			}
		}

		cTDEmGeometry(const cBlock& b) {
			set_zero();
			b.getvalue("tx_height", tx_height);
			b.getvalue("tx_roll", tx_roll);
			b.getvalue("tx_pitch", tx_pitch);
			b.getvalue("tx_yaw", tx_yaw);
			b.getvalue("txrx_dx", txrx_dx);
			b.getvalue("txrx_dy", txrx_dy);
			b.getvalue("txrx_dz", txrx_dz);
			b.getvalue("rx_roll", rx_roll);
			b.getvalue("rx_pitch", rx_pitch);
			b.getvalue("rx_yaw", rx_yaw);
		}

		inline static size_t size() {
			return 10;
		}

		double& operator[](const size_t& index)
		{
			switch (index) {
			case 0: return tx_height; break;
			case 1: return tx_roll; break;
			case 2: return tx_pitch; break;
			case 3: return tx_yaw; break;
			case 4: return txrx_dx; break;
			case 5: return txrx_dy; break;
			case 6: return txrx_dz; break;
			case 7: return rx_roll; break;
			case 8: return rx_pitch; break;
			case 9: return rx_yaw; break;
			default:
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
			}
			return tx_height;//this will never be reached
		}


		double operator[](const size_t& index) const
		{
			//Remove implied constness using const_cast
			return (*(const_cast<cTDEmGeometry*>(this)))[index]; // Correctly calls the function above.		
		};


		double& operator[](const std::string& gname)
		{
			const size_t& i = eindex(gname);
			return (*this)[i]; // Correctly calls the function above.
		};

		double operator[](const std::string& gname) const
		{
			//Remove implied constness using const_cast
			const size_t& i = eindex(gname);
			return (*(const_cast<cTDEmGeometry*>(this)))[i]; // Correctly calls the function above.
		};


		void set_zero() {
			for (size_t i = 0; i < size(); i++) {
				(*this)[i] = 0.0;
			}
		}

		void fillundefined(const cTDEmGeometry& g)
		{
			for (size_t i = 0; i < size(); i++) {
				if ((*this)[i] == undefinedvalue<double>()) {
					(*this)[i] = g[i];
				}
			}
		}

		static std::string element_name(const size_t& index) {

			switch (index) {
			case 0: return "tx_height"; break;
			case 1: return "tx_roll"; break;
			case 2: return "tx_pitch"; break;
			case 3: return "tx_yaw"; break;
			case 4: return "txrx_dx"; break;
			case 5: return "txrx_dy"; break;
			case 6: return "txrx_dz"; break;
			case 7: return "rx_roll"; break;
			case 8: return "rx_pitch"; break;
			case 9: return "rx_yaw"; break;
			default:
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
			}
			return "unknown";
		};

		static size_t eindex(const std::string& name) {

			for (size_t i = 0; i < size(); i++) {
				if (strcasecmp(name, element_name(i)) == 0) return i;
			}
			glog.errormsg(_SRC_, "Geometry field name %s is bad\n", name.c_str());
			return 0;
		};

		static std::string units(const size_t& index) {

			switch (index) {
			case 0: return "m"; break;
			case 1: return "degrees"; break;
			case 2: return "degrees"; break;
			case 3: return "degrees"; break;
			case 4: return "m"; break;
			case 5: return "m"; break;
			case 6: return "m"; break;
			case 7: return "degrees"; break;
			case 8: return "degrees"; break;
			case 9: return "degrees"; break;
			default:
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
				break;
			}
			return "unknown";
		};

		static std::string description(const size_t& index) {

			switch (index) {
			case 0: return "Tx height above ground level"; break;
			case 1: return "Tx roll - left side up + ve";   break;
			case 2: return "Tx pitch - nose down + ve";  break;
			case 3: return "Tx yaw - turn left + ve";    break;
			case 4: return "Tx - Rx horizonatl inline separation";   break;
			case 5: return "Tx - Rx horizonatl transverse separation";   break;
			case 6: return "Tx - Rx vertical separation";   break;
			case 7: return "Rx roll - left side up + ve";   break;
			case 8: return "Rx pitch - nose down + ve";  break;
			case 9: return "Rx yaw - turn left + ve";    break;
			default:
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
			}
			return "Error unknown geometry parameter";
		};

		static ElementType elementtype(const size_t& index) {
			switch (index) {
			case 0: return ElementType::tx_height; break;
			case 1: return ElementType::tx_roll;   break;
			case 2: return ElementType::tx_pitch;  break;
			case 3: return ElementType::tx_yaw;    break;
			case 4: return ElementType::txrx_dx;   break;
			case 5: return ElementType::txrx_dy;   break;
			case 6: return ElementType::txrx_dz;   break;
			case 7: return ElementType::rx_roll;   break;
			case 8: return ElementType::rx_pitch;  break;
			case 9: return ElementType::rx_yaw;    break;
			default:
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
				break;
			}
			return ElementType::unknown;
		}

		static cLEM::CalculationType derivativetype(const size_t& index) {
			switch (index) {
			case 0: return cLEM::CalculationType::HDERIVATIVE; break;
			case 1: return cLEM::CalculationType::NONE; break;
			case 2: return cLEM::CalculationType::NONE; break;
			case 3: return cLEM::CalculationType::NONE; break;
			case 4: return cLEM::CalculationType::XDERIVATIVE; break;
			case 5: return cLEM::CalculationType::YDERIVATIVE; break;
			case 6: return cLEM::CalculationType::ZDERIVATIVE; break;
			case 7: return cLEM::CalculationType::NONE; break;
			case 8: return cLEM::CalculationType::NONE; break;
			case 9: return cLEM::CalculationType::NONE; break;
			default:
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
				break;
			}
			return cLEM::CalculationType::NONE;
		}

		void write(std::string path) const
		{
			std::ofstream ofs(path);
			for (size_t i = 0; i < size(); i++) {
				ofs << element_name(i) << "\t" << (*this)[i] << std::endl;
			}
		}

		double txrx_dh() const {
			return std::hypot(txrx_dx, txrx_dy);
		};

		double txrx_dv() const {
			return std::hypot(txrx_dx, txrx_dz);
		};

		double txrx_dr() const {
			return std::sqrt(txrx_dx * txrx_dx + txrx_dy * txrx_dy + txrx_dz * txrx_dz);
		};

		cVec tx_orientation(const cVec& tx_reference_orientation) const {
			cVec v = tx_reference_orientation;
			v.rotate_inplace(tx_yaw, Geometry3D::zaxis);
			v.rotate_inplace(tx_pitch, Geometry3D::yaxis);
			v.rotate_inplace(tx_roll, Geometry3D::xaxis);
			return v;
		}

		cVec txrx_separation() const {
			cVec v = cVec(txrx_dx, txrx_dy, txrx_dz);
			return v;
		}

	};

	class LowPassFilter {

	private:
		double Order;
		double CutoffFrequency;

	public:

		LowPassFilter(const double order, const double cutoff_frequency) {
			Order = order;
			CutoffFrequency = cutoff_frequency;
		}

		cdouble weight(const double& frequency) const {
			cdouble a = 1.0 / cdouble(1.0, frequency / CutoffFrequency);
			return std::pow(a, Order);
		}

	};

	class Transmitter {

	public:
		double LoopArea = 0.0;
		double NumberOfTurns = 0.0;
		double PeakCurrent = 0.0;
		double PeakdIdT = 0.0;
		cVec Reference_Orientation = Geometry3D::zaxis;
	};

	class Waveform {

	public:
		enum class Type { TX, RX };

		double BaseFrequency = 0.0;
		double BasePeriod = 0.0;
		double SampleFrequency = 0.0;
		double SampleInterval = 0.0;
		size_t NumSamples = 0;
		size_t NumFrequencies = 0;

		Waveform::Type Type;
		std::vector<double> Time;
		//std::vector<double>  Current;
		//std::vector<double>  Received;
		std::vector<double> Value;
		std::vector<double>  T_Waveform;
		std::vector<cdouble> F_Waveform;
		std::vector<cdouble> Transfer;
		std::vector<cdouble> FFTWork;
		std::vector<double>  fft_frequency;

		void initialise(const cBlock& b, const fs::path& systemdescriptorfile) {
			BaseFrequency = b.getdoublevalue("BaseFrequency");
			BasePeriod = 1.0 / BaseFrequency;
			SampleFrequency = b.getdoublevalue("WaveformDigitisingFrequency");
			bool wavformdefined = false;

			if (wavformdefined == false) {
				std::string path = b.getstringvalue("WaveformReceived.File");
				if (isdefined(path)) {
					FilePathParts fpp(systemdescriptorfile);
					std::vector<std::vector<double>> wp = readwaveformfile(fpp.directory + path);
					if (wp.size() > 0) {
						digitisewaveform(wp, Time, Value);
						Type = Waveform::Type::RX;
						T_Waveform = Value;
						wavformdefined = true;
					}
				}
			}

			if (wavformdefined == false) {
				std::string path = b.getstringvalue("WaveformCurrent.File");
				if (isdefined(path)) {
					FilePathParts fpp(systemdescriptorfile);
					std::vector<std::vector<double>> wp = readwaveformfile(fpp.directory + path);
					if (wp.size() > 0) {
						digitisewaveform(wp, Time, Value);
						Type = Waveform::Type::TX;
						T_Waveform = Value;
						wavformdefined = true;
					}
				}
			}

			if (wavformdefined == false) {
				std::vector<std::vector<double>> wp = b.getdoublematrix("WaveformCurrent");
				if (wp.size() > 0) {
					digitisewaveform(wp, Time, Value);
					Type = Waveform::Type::TX;
					T_Waveform = Value;
					wavformdefined = true;
				}
			}

			if (wavformdefined == false) {
				std::vector<std::vector<double>> wp = b.getdoublematrix("WaveformReceived");
				if (wp.size() > 0) {
					digitisewaveform(wp, Time, Value);
					Type = Waveform::Type::RX;
					T_Waveform = Value;
					wavformdefined = true;
				}
			}

			if (wavformdefined == false) {
				glog.errormsg(_SRC_, "The waveform is not defined\n");
			}

		};

		double calculate_fft_frequency(const size_t index) const
		{
			double deltaF = 1.0 / ((double)NumSamples * SampleInterval);
			double s = (double)index;
			if (index > (NumSamples / 2)) {
				s = (double)index - (double)NumSamples;
			}
			return s * deltaF;
		}

		static std::vector<std::vector<double>> readwaveformfile(const std::string& filename)
		{
			std::vector<std::vector<double>> w;
			if (!fs::exists(filename)) {
				glog.errormsg(_SRC_, "\n\tD'Oh! the specified waveform file (%s) does not exist\n", filename.c_str());
			}

			std::ifstream ifs = ifstream_ex(filename);
			if (ifs.fail()) {
				glog.errormsg(_SRC_, "Unable to open waveformfile %s\n", filename.c_str());
			}

			std::string s;
			while (filegetline_ifs(ifs, s)) {
				trim_inplace(s);
				if (s.size() > 0) {
					std::vector<double> v(2);
					std::istringstream iss(s);
					iss >> v[0];
					iss >> v[1];
					w.push_back(v);
				}
			}
			return w;
		}

		void digitisewaveform(const std::vector<std::vector<double>>& wp, std::vector<double>& t, std::vector<double>& v) {
			double hp = 0.5 / BaseFrequency;
			SampleInterval = 1.0 / SampleFrequency;
			NumSamples = (size_t)(SampleFrequency / BaseFrequency);

			t.resize(NumSamples);
			v.resize(NumSamples);

			size_t np = wp.size();

			if (wp[np - 1][0] - wp[0][0] < hp) {
				glog.errormsg(_SRC_, "One complete halfcycle of the waveform has not been specified\n");

				glog.errormsg(_SRC_, "One complete halfcycle of the waveform has not been specified\n \
                                      Last waveform time - first waveform time must be >= 0.5/BaseFrequency\n");
			}

			std::vector<double> x(np);
			std::vector<double> y(np);
			for (size_t i = 0; i < np; i++) {
				x[i] = wp[i][0];
				y[i] = wp[i][1];
			}


			for (size_t i = 0; i < NumSamples / 2; i++) {
				bool set = false;

				double time = (double)i * SampleInterval;
				t[i] = time;

				if (time >= x[0] && time <= x[np - 1]) {
					v[i] = linearinterp(x, y, time);
					set = true;
				}
				else if ((time - hp) >= x[0] && (time - hp) <= x[np - 1]) {
					v[i] = -linearinterp(x, y, time - hp);
					set = true;
				}
				else if ((time - 2 * hp) >= x[0] && (time - 2 * hp) <= x[np - 1]) {
					v[i] = linearinterp(x, y, time - 2.0 * hp);
					set = true;
				}
				else if ((time + hp) >= x[0] && (time + hp) <= x[np - 1]) {
					v[i] = -linearinterp(x, y, time + hp);
					set = true;
				}
				else if ((time + 2 * hp) >= x[0] && (time + 2 * hp) <= x[np - 1]) {
					v[i] = linearinterp(x, y, time + 2.0 * hp);
					set = true;
				}

				if (set == false) {
					glog.errormsg(_SRC_, "Error in waveform - not all defined\n");
				}
			}

			for (size_t i = NumSamples / 2; i < NumSamples; i++) {
				t[i] = hp + t[i - NumSamples / 2];
				v[i] = -v[i - NumSamples / 2];
			}
		}

		double compute_peak_didt() const {
			double maxdidt = 0.0;
			for (size_t i = 1; i < Value.size(); i++) {
				const double di = Value[i] - Value[i - 1];
				const double dt = Time[i] - Time[i - 1];
				double didt = std::fabs(di / dt);
				if (didt > maxdidt) maxdidt = didt;
			}
			return maxdidt;
		}

		void write_timedomainwaveform(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumSamples; i++) {
				ofs << strprint("%20le\t%20le\n", Time[i], T_Waveform[i]);
			}
		}

		void write_frequencydomainwaveform(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumFrequencies; i++) {
				ofs << strprint("%15le\t%15le\t%15le\t%15le\t%15le\n", fft_frequency[i], F_Waveform[i].real(), F_Waveform[i].imag(), Transfer[i].real(), Transfer[i].imag());
			}
		}

		void write_frequencyseries(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumFrequencies; i++) {
				ofs << strprint("%15le\t%15le\t%15le\n", fft_frequency[i], FFTWork[i].real(), FFTWork[i].imag());
			}
		}
	};

	struct WindowSpecification {
		size_t SampleLow = 0;
		size_t SampleHigh = 0;
		size_t NumberOfSamples = 0;

		double TimeLow = 0.0;
		double TimeHigh = 0.0;
		double TimeWidth = 0.0;

		std::vector<size_t> Sample;
		std::vector<double> Weight;
	};

	struct WindowingScheme {

	private:


	public:

		size_t NumberOfWindows = 0;
		double WindowsTimeShift = 0.0;
		std::string WindowWeightingScheme;
		std::vector<WindowSpecification> WinSpec;

		WindowingScheme() {};

		WindowingScheme(const cBlock& receiverblock, const Waveform& WFM) {
			const cBlock& b = receiverblock;

			if (!b.getvalue("NumberOfWindows", NumberOfWindows)) {
				glog.errormsg(_SRC_, "NumberOfWindows is not specified");
			}

			if (!b.getvalue("TimeShift", WindowsTimeShift)) {
				WindowsTimeShift = 0.0;
			}

			WinSpec.resize(NumberOfWindows);
			//X.resize(NumberOfWindows);
			//Y.resize(NumberOfWindows);
			//Z.resize(NumberOfWindows);

			//Read window times
			std::vector<std::vector<double>> wt;
			if (!b.getvalue("WindowTimes", wt)) {
				glog.errormsg(_SRC_, "The WindowTimes have not been specified.");
			}
			size_t nw = wt.size();
			if (nw != NumberOfWindows) {
				glog.errormsg(_SRC_, "The number of WindowTimes does not match the NumberOfWindows\n");
			}

			for (size_t i = 0; i < NumberOfWindows; i++) {
				if (wt[i].size() != 2) {
					glog.errormsg(_SRC_, "The number of WindowTimes must have exactly 2 columns (error in window %lu)\n", i + 1);
				}
				WinSpec[i].TimeLow = wt[i][0] + WindowsTimeShift;
				WinSpec[i].TimeHigh = wt[i][1] + WindowsTimeShift;
			}

			WindowWeightingScheme = b.getstringvalue("WindowWeightingScheme");
			if (strcasecmp(WindowWeightingScheme, "AreaUnderCurve") == 0) {
				initialise_area(WFM);
			}
			else if (strcasecmp(WindowWeightingScheme, "Boxcar") == 0) {
				initialise_boxcar(WFM);
			}
			else if (strcasecmp(WindowWeightingScheme, "LinearTaper") == 0) {
				initialise_lineartaper(WFM);
			}
			else {
				glog.errormsg(_SRC_, "WindowWeightingScheme %s unknown (must be \"AreaUnderCurve\" or  \"Boxcar\" or \"LinearTaper\")\n", WindowWeightingScheme.c_str());
			}
		}

		const size_t& nwindows() const { return NumberOfWindows; }

		void initialise_area(const Waveform& WFM)
		{
			double tlow, thigh, t, tp, tn, tleft, tright;

			double dwt = WFM.Time[1] - WFM.Time[0];
			double eps = 1.0e-7;

			for (size_t w = 0; w < NumberOfWindows; w++) {
				WinSpec[w].TimeWidth = WinSpec[w].TimeHigh - WinSpec[w].TimeLow;
				for (size_t s = 0; s < WFM.NumSamples; s++) {
					t = WFM.Time[s];
					if (t + eps >= WinSpec[w].TimeLow) {
						WinSpec[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = WFM.NumSamples; s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					t = WFM.Time[s];
					if (t - eps <= WinSpec[w].TimeHigh) {
						WinSpec[w].SampleHigh = s;
						break;
					}
				}

				WinSpec[w].NumberOfSamples = WinSpec[w].SampleHigh - WinSpec[w].SampleLow + 1;
				WinSpec[w].Sample.resize(WinSpec[w].NumberOfSamples);
				WinSpec[w].Weight.resize(WinSpec[w].NumberOfSamples);
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Sample[k] = WinSpec[w].SampleLow + k;
					WinSpec[w].Weight[k] = 0.0;
				}

				tlow = WinSpec[w].TimeLow;
				thigh = WinSpec[w].TimeHigh;
				double wsum = 0;
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					size_t s = WinSpec[w].Sample[k];
					t = WFM.Time[s];
					tp = WFM.Time[s] - dwt;
					tn = WFM.Time[s] + dwt;
					tleft = std::max(tp, tlow);
					tright = std::min(tn, thigh);
					WinSpec[w].Weight[k] = 0.5 * (t - tleft) + 0.5 * (tright - t);
					WinSpec[w].Weight[k] /= WinSpec[w].TimeWidth;
					wsum += WinSpec[w].Weight[k];
				}
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Weight[k] /= wsum;
				}
			}
		}

		void initialise_boxcar(const Waveform& WFM)
		{
			double eps = 1.0e-7;
			for (size_t w = 0; w < NumberOfWindows; w++) {
				WinSpec[w].TimeWidth = WinSpec[w].TimeHigh - WinSpec[w].TimeLow;
				for (size_t s = 0; s < WFM.NumSamples; s++) {
					double t = WFM.Time[s];
					if (t + eps >= WinSpec[w].TimeLow) {
						WinSpec[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = WFM.NumSamples; s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					double t = WFM.Time[s];
					if (t - eps <= WinSpec[w].TimeHigh) {
						WinSpec[w].SampleHigh = s;
						break;
					}
				}


				WinSpec[w].NumberOfSamples = WinSpec[w].SampleHigh - WinSpec[w].SampleLow + 1;
				WinSpec[w].Sample.resize(WinSpec[w].NumberOfSamples);
				WinSpec[w].Weight.resize(WinSpec[w].NumberOfSamples);
				double weightsum = 0.0;
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Sample[k] = WinSpec[w].SampleLow + k;
					WinSpec[w].Weight[k] = 1.0;
					weightsum += WinSpec[w].Weight[k];
				}
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Weight[k] /= weightsum;
				}
				//printf("%lu %lu\n", w, WinSpec[w].NumberOfSamples);
			}
		}

		void initialise_lineartaper(const Waveform& WFM)
		{
			double eps = 1.0e-7;
			for (size_t w = 0; w < NumberOfWindows; w++) {
				WinSpec[w].TimeWidth = WinSpec[w].TimeHigh - WinSpec[w].TimeLow;
				for (size_t s = 0; s < WFM.NumSamples; s++) {
					double t = WFM.Time[s];
					if (t + eps >= WinSpec[w].TimeLow) {
						WinSpec[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = WFM.NumSamples; s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					double t = WFM.Time[s];
					if (t - eps <= WinSpec[w].TimeHigh) {
						WinSpec[w].SampleHigh = s;
						break;
					}
				}

				size_t ns = WinSpec[w].SampleHigh - WinSpec[w].SampleLow + 1;
				WinSpec[w].NumberOfSamples = ns * 3;
				WinSpec[w].SampleLow -= ns;
				WinSpec[w].SampleHigh += ns;
				WinSpec[w].Sample.resize(WinSpec[w].NumberOfSamples);
				WinSpec[w].Weight.resize(WinSpec[w].NumberOfSamples);

				double weightsum = 0.0;
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Sample[k] = WinSpec[w].SampleLow + k;
					if (k < ns) {
						WinSpec[w].Weight[k] = (double)(k + 1) / (double)(ns + 1);
					}
					else if (k >= 2 * ns) {
						WinSpec[w].Weight[k] = 1.0 - (double)((k + 1) - 2 * ns) / (double)(ns + 1);
					}
					else {
						WinSpec[w].Weight[k] = 1.0;
					}
					weightsum += WinSpec[w].Weight[k];
				}
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Weight[k] /= weightsum;
				}
			}
		}

		std::vector<double> computewindow(const double* timeseries)
		{
			std::vector<double> W(NumberOfWindows, 0.0);
			for (size_t w = 0; w < NumberOfWindows; w++) {
				for (size_t k = 0; k < WinSpec[w].Sample.size(); k++) {
					W[w] += timeseries[WinSpec[w].Sample[k]] * WinSpec[w].Weight[k];
				}
			}
			return W;
		}

		void printwindows(const double& PX, const double& PY, const double& PZ, const std::vector<double>& SX, const std::vector<double>& SY, const std::vector<double>& SZ) const {
			printf("Primary   %15.8lf%15.8lf%15.8lf\n\n", PX, PY, PZ);
			printf("Window#             X               Y               Z\n");
			for (size_t w = 0; w < NumberOfWindows; w++) {
				printf("%2zu        %15.8lf%15.8lf%15.8lf\n", w + 1, SX[w], SY[w], SZ[w]);
			}
		};

		void write_windows(const fs::path& path, const std::vector<double>& SX, const std::vector<double>& SY, const std::vector<double>& SZ) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t w = 0; w < NumberOfWindows; w++) {
				ofs << strprint("%2zu\t%20e\t%20e\t%15e%15e%15e\n", w + 1, WinSpec[w].TimeLow, WinSpec[w].TimeHigh, SX[w], SY[w], SZ[w]);
			}
		};
	};

	class AEMSystem {

	protected:
		std::string SystemName;
		std::string SystemType;
		cLEM LEM;
		cBlock STM;

		bool SaveDiagnosticFiles = false;

	public:

		inline static const size_t XCOMP = 0;
		inline static const size_t YCOMP = 1;
		inline static const size_t ZCOMP = 2;
		inline static const size_t NCOMP = 3;

		const cBlock& stm() const { return STM; };
		cLEM& lem() { return LEM; };

	};

	class ComponentWorkStore {

	public:

		std::vector<double> R;
		std::vector<double> I;
		std::vector<cdouble> complex_splined_values;
		std::vector<double> windows;

		ComponentWorkStore() {};
		void resize(const size_t nnodes, const size_t nfftfreq) {
			R.resize(nnodes);
			I.resize(nnodes);
			complex_splined_values.resize(nfftfreq);
		};
	};

	class cTDEmSystem : public AEMSystem {

	private:
		std::vector<ComponentWorkStore> Comp;
		std::vector<double> HxR;
		std::vector<double> HxI;
		std::vector<double> HyR;
		std::vector<double> HyI;
		std::vector<double> HzR;
		std::vector<double> HzI;

		FixedPointSpline<double> FSPLINE;
		std::vector<double> HxR_spline;
		std::vector<double> HxI_spline;
		std::vector<double> HyR_spline;
		std::vector<double> HyI_spline;
		std::vector<double> HzR_spline;
		std::vector<double> HzI_spline;



		std::vector<double> h2_spline;
		std::vector<double> a_spline;
		std::vector<double> b_spline;
		std::vector<size_t> klo_spline;
		std::vector<size_t> khi_spline;
		std::vector<double> a3ma_spline;
		std::vector<double> b3mb_spline;

		std::vector<cdouble> X_splined;
		std::vector<cdouble> Y_splined;
		std::vector<cdouble> Z_splined;

	public:

		enum class OutputType { BFIELD, DBDT };
		enum class NormalizationType { NONE, PPM, PPM_PEAKTOPEAK };

		OutputType OutputType;
		NormalizationType Normalisation;

		Waveform WFM;
		fftw_plan inverse_fftplan;
		size_t FrequenciesPerDecade = 0;
		size_t NumberOfSplinedFrequencies = 0;

		size_t NumberOfDiscreteFrequencies = 0;
		std::vector<double> DiscreteFrequencies;
		std::vector<double> DiscreteFrequenciesLog10;
		double DiscreteFrequencyLow = 0.0;
		double DiscreteFrequencyHigh = 0.0;
		double FrequencyLog10Spacing = 0.0;
		std::vector<double> SplinedFrequencieslog10;;

		std::vector<LowPassFilter> Filters;

		Transmitter   Tx;
		cTDEmGeometry Geometry;
		cTDEmGeometry NormalizationGeometry;

		std::vector<double> Scale;
		std::vector<double> X; //Secondary X field
		std::vector<double> Y; //Secondary Y field 
		std::vector<double> Z; //Secondary Z field 
		double PrimaryX = 0.0;  //Primary X field
		double PrimaryY = 0.0;  //Primary Y field
		double PrimaryZ = 0.0;  //Primary Z field
		double RefGeomPrimaryX = 0.0;  //Primary X ref field for PPM normalisation
		double RefGeomPrimaryY = 0.0;  //Primary Y ref field for PPM normalisation
		double RefGeomPrimaryZ = 0.0;  //Primary Z ref field for PPM normalisation

		WindowingScheme Win;

		cTDEmSystem() { initialise(); };

		cTDEmSystem(std::string systemdescriptorfile) {
			initialise();
			read_system_descriptor_file(systemdescriptorfile);
		};

		~cTDEmSystem()
		{
			if (inverse_fftplan) {
				fftw_destroy_plan(inverse_fftplan);
			}
		};

		const size_t& nwindows() const {
			return Win.nwindows();
		}

		double primary(const size_t component) const {
			if (component == 0) return PrimaryX;
			else if (component == 1) return PrimaryY;
			else if (component == 2) return PrimaryZ;
			else return 0;
		}

		double secondary(const size_t component, const size_t window) {
			if (component == 0) return X[window];
			else if (component == 1) return Y[window];
			else if (component == 2) return Z[window];
			else return 0;
		}

		std::vector<double> secondary(const size_t component) const {
			if (component == 0) return X;
			else if (component == 1) return Y;
			else if (component == 2) return Z;
			else return std::vector<double>(0);
		}

		void initialise() {
			LEM.calculation_type = cLEM::CalculationType::FORWARDMODEL;
			LEM.rzerotype = cLEM::RZeroMethod::PROPOGATIONMATRIX;
			inverse_fftplan = 0;
		}

		void read_system_descriptor_file(const std::string& systemdescriptorfile)
		{
			if (!fs::exists(systemdescriptorfile)) {
				std::string msg = strprint("\n\tD'Oh! the specified system descriptor file (%s) does not exist\n", systemdescriptorfile.c_str());
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

			WFM.initialise(b, systemdescriptorfile);

			/*
			BaseFrequency = b.getdoublevalue("BaseFrequency");
			BasePeriod = 1.0 / BaseFrequency;
			SampleFrequency = b.getdoublevalue("WaveformDigitisingFrequency");
			bool wavformdefined = false;
			if (wavformdefined == false) {
				std::string path = b.getstringvalue("WaveformReceived.File");
				if (isdefined(path)) {
					FilePathParts fpp(systemdescriptorfile);
					std::vector<std::vector<double>> wp = readwaveformfile(fpp.directory + path);
					if (wp.size() > 0) {
						digitisewaveform(wp, WaveformTime, WaveformReceived);
						WaveformType = Waveform::Type::RX;
						T_Waveform = WaveformReceived;
						wavformdefined = true;
					}
				}
			}

			if (wavformdefined == false) {
				std::string path = b.getstringvalue("WaveformCurrent.File");
				if (isdefined(path)) {
					FilePathParts fpp(systemdescriptorfile);
					std::vector<std::vector<double>> wp = readwaveformfile(fpp.directory + path);
					if (wp.size() > 0) {
						digitisewaveform(wp, WaveformTime, WaveformCurrent);
						WaveformType = Waveform::Type::TX;
						T_Waveform = WaveformCurrent;
						wavformdefined = true;
					}
				}
			}

			if (wavformdefined == false) {
				std::vector<std::vector<double>> wp = b.getdoublematrix("WaveformCurrent");
				if (wp.size() > 0) {
					digitisewaveform(wp, WaveformTime, WaveformCurrent);
					WaveformType = Waveform::Type::TX;
					T_Waveform = WaveformCurrent;
					wavformdefined = true;
				}
			}

			if (wavformdefined == false) {
				std::vector<std::vector<double>> wp = b.getdoublematrix("WaveformReceived");
				if (wp.size() > 0) {
					digitisewaveform(wp, WaveformTime, WaveformReceived);
					WaveformType = Waveform::Type::RX;
					T_Waveform = WaveformReceived;
					wavformdefined = true;
				}
			}

			if (wavformdefined == false) {
				glog.errormsg(_SRC_, "The waveform is not defined\n");
			}
			*/

			cBlock rxblock = STM.findblock("Receiver");
			initialise_windows(rxblock);

			LEM.ModellingLoopRadius = STM.getdoublevalue("ForwardModelling.ModellingLoopRadius");
			if (!isdefined(LEM.ModellingLoopRadius)) {
				LEM.ModellingLoopRadius = 0.0;
			}

			std::string ot = STM.getstringvalue("ForwardModelling.OutputType");
			if (strcasecmp(ot, "B") == 0) {
				OutputType = OutputType::BFIELD;
			}
			else if (strcasecmp(ot, "dB/dt") == 0) {
				OutputType = OutputType::DBDT;
			}
			else {
				glog.errormsg(_SRC_, "OutputType %s unknown (must be one of \"B\" or \"dB/dt\")\n", ot.c_str());
			}

			FrequenciesPerDecade = (size_t)STM.getintvalue("ForwardModelling.FrequenciesPerDecade");
			if (FrequenciesPerDecade < 5) {
				glog.warningmsg(_SRC_, "It is wise to use at least 5 frequencies per decade\n");
			}

			LEM.NumAbscissa = (size_t)STM.getintvalue("ForwardModelling.NumberOfAbsiccaInHankelTransformEvaluation");

			std::string n = STM.getstringvalue("ForwardModelling.SecondaryFieldNormalisation");
			if (strcasecmp(n, "None") == 0) {
				Normalisation = NormalizationType::NONE;
			}
			else if (strcasecmp(n, "PPM") == 0) {
				Normalisation = NormalizationType::PPM;
			}
			else if (strcasecmp(n, "PPMPEAKTOPEAK") == 0) {
				Normalisation = NormalizationType::PPM_PEAKTOPEAK;
			}
			else {
				glog.errormsg(_SRC_, "Normalisation %s unknown (must be one of \"None,PPM,PPMPEAKTOPEAK\")\n", n.c_str());
			}

			SaveDiagnosticFiles = STM.getboolvalue("ForwardModelling.SaveDiagnosticFiles");

			if (WFM.Time.size() <= 2 || WFM.Time.size() != WFM.T_Waveform.size()) {
				glog.errormsg(_SRC_, "The number of WaveformTime values must match number of WaveformCurrent/WaveformReceived values and also be more than two\n");
			}

			if (LEM.NumAbscissa < 17) {
				glog.warningmsg(_SRC_, "It is wise to use at least 17 Absicca for integrating the Hankel Transforms");
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


			systeminitialise();
		};

		void systeminitialise() {
			create_transforms();
			setupdiscretefrequencies();
			setup_splines();
			setup_scaling();
		}

		void create_transforms() {
			WFM.NumFrequencies = WFM.NumSamples / 2 + 1;
			WFM.fft_frequency.resize(WFM.NumFrequencies);

			size_t N = WFM.NumSamples;
			size_t NC = N;
			size_t NR = 2 * (N / 2 + 1);

			//Forward transform
			WFM.F_Waveform.resize(NC);
			double* in = (double*)&(WFM.T_Waveform[0]);
			fftw_complex* out = (fftw_complex*)&(WFM.F_Waveform[0]);
			fftw_plan fftwplan_forward = fftw_plan_dft_r2c_1d((int)N, in, out, FFTW_ESTIMATE);
			for (size_t k = 0; k < WFM.NumSamples; k++) {
				WFM.T_Waveform[k] /= (double)WFM.NumSamples;
			}
			fftw_execute(fftwplan_forward);
			fftw_destroy_plan(fftwplan_forward);

			for (size_t k = 0; k < WFM.NumFrequencies; k++) {
				WFM.fft_frequency[k] = WFM.calculate_fft_frequency(k);
			}

			bool convert_B_2_dBdT = false;
			bool convert_dBdT_2_B = false;
			if (WFM.Type == Waveform::Type::TX) {
				if (OutputType == OutputType::DBDT) {
					convert_B_2_dBdT = true;
				}
			}

			if (WFM.Type == Waveform::Type::RX) {
				if (OutputType == OutputType::BFIELD) {
					convert_dBdT_2_B = true;
				}
			}

			for (size_t k = 0; k < WFM.NumFrequencies; k += 2) {
				WFM.F_Waveform[k] = 0.0;
			}

			WFM.Transfer.resize(WFM.NumFrequencies);
			WFM.Transfer = WFM.F_Waveform;

			for (size_t k = 0; k < WFM.NumFrequencies; k++) {
				const double& frequency = WFM.fft_frequency[k];
				if (convert_B_2_dBdT == true) {
					WFM.Transfer[k] *= cdouble(0.0, -TWOPI<double> *frequency);
				}
				if (convert_dBdT_2_B == true) {
					WFM.Transfer[k] *= cdouble(0.0, -1.0 / (TWOPI<double> *frequency));
				}

				for (size_t fi = 0; fi < Filters.size(); fi++) {
					const cdouble w = Filters[fi].weight(frequency);
					WFM.Transfer[k] *= w;
				}
			}

			//Frequencies to be splined
			NumberOfSplinedFrequencies = WFM.NumFrequencies / 2;
			SplinedFrequencieslog10.resize(NumberOfSplinedFrequencies);
			for (size_t k = 0; k < NumberOfSplinedFrequencies; k++) {
				SplinedFrequencieslog10[k] = log10(fabs(WFM.fft_frequency[k * 2 + 1]));
			}

			//Setup inverse transform work array	
			WFM.FFTWork.resize(NR);

#if defined MULTITHREADED
			//FFTW_MEASURE does not seem to be thread safe
			unsigned int FLAGS = FFTW_ESTIMATE;
#else
			unsigned int FLAGS = FFTW_MEASURE;
#endif

			fftw_complex* invin = (fftw_complex*)(WFM.FFTWork.data());
			double* invout = (double*)(WFM.FFTWork.data());
			inverse_fftplan = fftw_plan_dft_c2r_1d((int)N, invin, invout, FLAGS);
		}

		void setup_splines() {
			Comp.resize(3);
			Comp[XCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies);
			Comp[YCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies);
			Comp[ZCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies);

			HxR = std::vector<double>(NumberOfDiscreteFrequencies);
			HxI = std::vector<double>(NumberOfDiscreteFrequencies);
			HyR = std::vector<double>(NumberOfDiscreteFrequencies);
			HyI = std::vector<double>(NumberOfDiscreteFrequencies);
			HzR = std::vector<double>(NumberOfDiscreteFrequencies);
			HzI = std::vector<double>(NumberOfDiscreteFrequencies);

			HxR_spline = std::vector<double>(NumberOfDiscreteFrequencies);
			HxI_spline = std::vector<double>(NumberOfDiscreteFrequencies);
			HyR_spline = std::vector<double>(NumberOfDiscreteFrequencies);
			HyI_spline = std::vector<double>(NumberOfDiscreteFrequencies);
			HzR_spline = std::vector<double>(NumberOfDiscreteFrequencies);
			HzI_spline = std::vector<double>(NumberOfDiscreteFrequencies);

			size_t ns = NumberOfSplinedFrequencies;
			h2_spline = std::vector<double>(ns);
			a_spline = std::vector<double>(ns);
			b_spline = std::vector<double>(ns);
			klo_spline = std::vector<size_t>(ns);
			khi_spline = std::vector<size_t>(ns);
			a3ma_spline = std::vector<double>(ns);
			b3mb_spline = std::vector<double>(ns);

			X_splined = std::vector<cdouble>(ns);
			Y_splined = std::vector<cdouble>(ns);
			Z_splined = std::vector<cdouble>(ns);

			FSPLINE.initialise(DiscreteFrequenciesLog10, SplinedFrequencieslog10);

			setup_splineinterp(DiscreteFrequenciesLog10, SplinedFrequencieslog10);

			LEM.init_frequencies(DiscreteFrequencies);
		}

		void spline(const std::vector<double>& x, const std::vector<double>& y, double yp1, double ypn, std::vector<double>& y2) {
			double p, qn, sig, un;

			size_t n = y2.size();
			std::vector<double> u(n - 1);
			if (yp1 > 0.99e30)
				y2[0] = u[0] = 0.0;
			else {
				y2[0] = -0.5;
				u[0] = (3.0 / (x[1] - x[0])) * ((y[1] - y[0]) / (x[1] - x[0]) - yp1);
			}
			for (size_t i = 1; i < n - 1; i++) {
				sig = (x[i] - x[i - 1]) / (x[i + 1] - x[i - 1]);
				p = sig * y2[i - 1] + 2.0;
				y2[i] = (sig - 1.0) / p;
				u[i] = (y[i + 1] - y[i]) / (x[i + 1] - x[i]) - (y[i] - y[i - 1]) / (x[i] - x[i - 1]);
				u[i] = (6.0 * u[i] / (x[i + 1] - x[i - 1]) - sig * u[i - 1]) / p;
			}
			if (ypn > 0.99e30)
				qn = un = 0.0;
			else {
				qn = 0.5;
				un = (3.0 / (x[n - 1] - x[n - 2])) * (ypn - (y[n - 1] - y[n - 2]) / (x[n - 1] - x[n - 2]));
			}
			y2[n - 1] = (un - qn * u[n - 2]) / (qn * y2[n - 2] + 1.0);

			for (size_t k = n - 1; k-- > 0;) {
				//Note the unusual syntax for decrement of unsigned variable
				y2[k] = y2[k] * y2[k + 1] + u[k];
			}

			//for (size_t k = n - 2; k >= 0; k--){
			//	y2[k] = y2[k] * y2[k + 1] + u[k];
			//	if (k == 0)break;
			//}
		}

		void setup_splineinterp(const std::vector<double>& xn, const std::vector<double>& xi) {
			const size_t n = xn.size();
			for (size_t fi = 0; fi < NumberOfSplinedFrequencies; fi++) {
				size_t klo = 0;
				size_t khi = n - 1;

				while (khi - klo > 1) {
					size_t k = (khi + klo) >> 1;
					if (xn[k] > xi[fi]) {
						khi = k;
					}
					else klo = k;
				}
				double h = xn[khi] - xn[klo];
				if (h == 0.0) {
					glog.errormsg(_SRC_, "setup_splineinterp(): Bad xa input to routine splintsetup\n");
				}
				double a = (xn[khi] - xi[fi]) / h;
				double b = (xi[fi] - xn[klo]) / h;

				a_spline[fi] = a;
				b_spline[fi] = b;
				h2_spline[fi] = h * h;
				klo_spline[fi] = klo;
				khi_spline[fi] = khi;
				a3ma_spline[fi] = a * a * a - a;
				b3mb_spline[fi] = b * b * b - b;
			}

		}

		void spline_interp() {
			//y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;			
			for (size_t k = 0; k < NumberOfSplinedFrequencies; k++) {
				const double& a = a_spline[k];
				const double& b = b_spline[k];
				const double& a3 = a3ma_spline[k];
				const double& b3 = b3mb_spline[k];
				const size_t& klo = klo_spline[k];
				const size_t& khi = khi_spline[k];
				const double& scale = (h2_spline[k]) / 6.0;

				if (Scale[XCOMP] != 0.0) {
					double rv, iv;
					rv = (a * HxR[klo] + b * HxR[khi] + scale * (a3 * HxR_spline[klo] + b3 * HxR_spline[khi])),
						iv = (a * HxI[klo] + b * HxI[khi] + scale * (a3 * HxI_spline[klo] + b3 * HxI_spline[khi]));
					X_splined[k] = cdouble(rv, iv);
				}

				if (Scale[YCOMP] != 0.0) {
					double rv, iv;
					rv = (a * HyR[klo] + b * HyR[khi] + scale * (a3 * HyR_spline[klo] + b3 * HyR_spline[khi]));
					iv = (a * HyI[klo] + b * HyI[khi] + scale * (a3 * HyI_spline[klo] + b3 * HyI_spline[khi]));
					Y_splined[k] = cdouble(rv, iv);
				}

				if (Scale[ZCOMP] != 0.0) {
					double rv, iv;
					rv = (a * HzR[klo] + b * HzR[khi] + scale * (a3 * HzR_spline[klo] + b3 * HzR_spline[khi]));
					iv = (a * HzI[klo] + b * HzI[khi] + scale * (a3 * HzI_spline[klo] + b3 * HzI_spline[khi]));
					Z_splined[k] = cdouble(rv, iv);
				}
			}
		}

		void inversefft() const {
			fftw_execute(inverse_fftplan);
		}

		void setearthproperties(const cEarth1D& E)
		{
			LEM.setproperties(E);
		}

		void setconductivitythickness(const size_t nlayers, const double* conductivity, const double* thickness)
		{
			LEM.setconductivitythickness(nlayers, conductivity, thickness);
		}

		void setconductivitythickness(const std::vector<double>& conductivity, const std::vector<double>& thickness)
		{
			LEM.setconductivitythickness(conductivity, thickness);
		}

		void setgeometry(const cTDEmGeometry& G)
		{
			//X = +ve in flight direction
			//Y = +ve on left wing
			//Z = +ve vertical up
			//ie different to Fugro convention

			//Steer left is positive yaw     X->Y axis
			//Left wing up is positive roll  Y->Z axis
			//Nose down is positive pitch	 Z->X axis

			Geometry = G;
			cVec tx_reference_orientation = Geometry3D::zaxis;

			const cVec sep = Geometry.txrx_separation();
			const double& h = Geometry.tx_height;
			const double& x = sep.x;
			const double& y = sep.y;
			const double& z = h + sep.z;
			const cVec tx_orientation = Geometry.tx_orientation(Tx.Reference_Orientation);
			LEM.setgeometry(tx_orientation, h, x, y, z);

			//RX_height = TX_height + G.txrx_dz;
			//RX_roll = G.rx_roll;
			//RX_pitch = G.rx_pitch;
			//RX_yaw = G.rx_yaw;
		};

		void setgeometry(const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw)
		{
			cTDEmGeometry G(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
			setgeometry(G);
		}

		void setgeometry(const double* g)
		{
			//const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw)
			cTDEmGeometry G(g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9]);
			setgeometry(G);
		}

		void setupcomputations()
		{
			for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++) {
				LEM.init_frequency(fi);
			}
		}

		void setprimaryfields()
		{
			LEM.setprimaryfields();
			PrimaryX = LEM.Fields.t.p.x;
			PrimaryY = LEM.Fields.t.p.y;
			PrimaryZ = LEM.Fields.t.p.z;

			if (LEM.calculation_type == cLEM::CalculationType::HDERIVATIVE) {
				//This is because when H changes Z also changes
				//and DZ = DH
				//but they should be all zero anyway
				PrimaryX *= 2.0;
				PrimaryY *= 2.0;
				PrimaryZ *= 2.0;
			}

			if (Normalisation == NormalizationType::PPM_PEAKTOPEAK) {
				PrimaryX *= 2.0;
				PrimaryY *= 2.0;
				PrimaryZ *= 2.0;
			}

			if (OutputType == OutputType::DBDT) {
				//Must convert to dB/dt. This happens implicitly for the secondary via the waveform.
				PrimaryX *= Tx.PeakdIdT;
				PrimaryY *= Tx.PeakdIdT;
				PrimaryZ *= Tx.PeakdIdT;
			}

			cVec field = cVec(PrimaryX, PrimaryY, PrimaryZ);
			rotate_field_to_receiver_frame(field);
			PrimaryX = field.x;
			PrimaryY = field.y;
			PrimaryZ = field.z;

			PrimaryX *= Scale[XCOMP];
			PrimaryY *= Scale[YCOMP];
			PrimaryZ *= Scale[ZCOMP];

		}

		void rotate_field_to_receiver_frame(cVec& fvec)
		{
			//Rotating in opposite sense because we are rotating the axes and doing it in the reverse order
			fvec.rotate_inplace(-Geometry.rx_roll, Geometry3D::xaxis);
			fvec.rotate_inplace(-Geometry.rx_pitch, Geometry3D::yaxis);
			fvec.rotate_inplace(-Geometry.rx_yaw, Geometry3D::zaxis);
		}

		void spline_component(const size_t& component) {
			ComponentWorkStore& C = Comp[component];
			const std::vector<double>& v = FSPLINE.interpolated_values();

			const size_t n = v.size();
			double* a = (double*)(C.complex_splined_values.data());

			FSPLINE.compute_interpolation(C.R);
			for (size_t i = 0; i < n; i++) {
				a[i * 2] = v[i];
			}

			FSPLINE.compute_interpolation(C.I);
			for (size_t i = 0; i < n; i++) {
				a[i * 2 + 1] = v[i];
			}
		};

		void inverse_fft_window_scale_component(const size_t& component) {

			ComponentWorkStore& C = Comp[component];
			// Reset to the stored transfer function
			WFM.FFTWork = WFM.Transfer;
			// Apply transfer function
			size_t n = 0;
			for (size_t k = 1; k < WFM.NumFrequencies; k += 2) {
				WFM.FFTWork[k] *= C.complex_splined_values[n];
				n++;
			}

			// Inverse FFT
			fftw_execute(inverse_fftplan);

			// Window
			C.windows = Win.computewindow((double*)WFM.FFTWork.data());
			if (SaveDiagnosticFiles) {
				write_timesseries("diag_xtimeseries.txt");
			}

			// Scale
			C.windows *= Scale[component];
		}

		void setsecondaryfields() {
			//Computation for discrete frequencies 	
			for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++) {
				LEM.dointegrals(fi);
				LEM.setsecondaryfields(fi);
				const cdouble& x = LEM.Fields.t.s.x;
				const cdouble& y = LEM.Fields.t.s.y;
				const cdouble& z = LEM.Fields.t.s.z;
				cVec vr = cVec(x.real(), y.real(), z.real());
				cVec vi = cVec(x.imag(), y.imag(), z.imag());

				rotate_field_to_receiver_frame(vr);
				rotate_field_to_receiver_frame(vi);

				if (LEM.calculation_type == cLEM::CalculationType::HDERIVATIVE) {
					//This is because when H changes Z also changes
					//and DZ = DH
					vr *= 2.0;
					vi *= 2.0;
				}

				Comp[XCOMP].R[fi] = vr.x;
				Comp[XCOMP].I[fi] = vi.x;
				Comp[YCOMP].R[fi] = vr.y;
				Comp[YCOMP].I[fi] = vi.y;
				Comp[ZCOMP].R[fi] = vr.z;
				Comp[ZCOMP].I[fi] = vi.z;

				HxR[fi] = vr.x;
				HxI[fi] = vi.x;
				HyR[fi] = vr.y;
				HyI[fi] = vi.y;
				HzR[fi] = vr.z;
				HzI[fi] = vi.z;
			}

			for (size_t i = 0; i < NCOMP; i++) {
				if (Scale[i] == 0.0) return;
				spline_component(i);
			}

			for (size_t i = 0; i < NCOMP; i++) {
				if (Scale[i] == 0.0) return;
				inverse_fft_window_scale_component(i);
			}

			//Spline discreet frequencies		
			if (Scale[XCOMP] != 0.0) {
				spline(DiscreteFrequenciesLog10, HxR, 1e-30, 1e-30, HxR_spline);
				spline(DiscreteFrequenciesLog10, HxI, 1e-30, 1e-30, HxI_spline);
			}

			if (Scale[YCOMP] != 0.0) {
				spline(DiscreteFrequenciesLog10, HyR, 1e-30, 1e-30, HyR_spline);
				spline(DiscreteFrequenciesLog10, HyI, 1e-30, 1e-30, HyI_spline);
			}

			if (Scale[ZCOMP] != 0.0) {
				spline(DiscreteFrequenciesLog10, HzR, 1e-30, 1e-30, HzR_spline);
				spline(DiscreteFrequenciesLog10, HzI, 1e-30, 1e-30, HzI_spline);
			}

			//Interpolate
			spline_interp();

			if (SaveDiagnosticFiles) {
				write_discretefrequencies("diag_discretefrequencies.txt");
				write_splinedfrequencies("diag_splinedfrequencies.txt");
				WFM.write_frequencydomainwaveform("diag_frequencydomainwaveform.txt");
			}

			//Inverse FFT
			if (Scale[XCOMP] != 0.0) {
				size_t n = 0;
				WFM.FFTWork = WFM.Transfer;
				for (size_t k = 1; k < WFM.NumFrequencies; k += 2) {
					WFM.FFTWork[k] *= X_splined[n];
					n++;
				}
				//Inverse FFT
				fftw_execute(inverse_fftplan);
				X = Win.computewindow((double*)WFM.FFTWork.data());
				if (SaveDiagnosticFiles) {
					write_timesseries("diag_xtimeseries.txt");
				}
				X *= Scale[XCOMP];
			}

			if (Scale[YCOMP] != 0.0) {
				size_t n = 0;
				WFM.FFTWork = WFM.Transfer;
				for (size_t k = 1; k < WFM.NumFrequencies; k += 2) {
					WFM.FFTWork[k] = Y_splined[n];
					n++;
				}
				//Inverse FFT		
				fftw_execute(inverse_fftplan);
				Y = Win.computewindow((double*)WFM.FFTWork.data());
				if (SaveDiagnosticFiles) {
					write_timesseries("diag_ytimeseries.txt");
				}
				Y *= Scale[YCOMP];
			}

			if (Scale[ZCOMP] != 0.0) {
				size_t n = 0;
				WFM.FFTWork = WFM.Transfer;
				for (size_t k = 1; k < WFM.NumFrequencies; k += 2) {
					WFM.FFTWork[k] *= Z_splined[n];
					n++;
				}

				//Inverse FFT
				fftw_execute(inverse_fftplan);
				Z = Win.computewindow((double*)WFM.FFTWork.data());
				if (SaveDiagnosticFiles) {
					write_timesseries("diag_ztimeseries.txt");
				}
				Z *= Scale[ZCOMP];
			}

			for (size_t w = 0; w < Win.NumberOfWindows; w++) {
				//std::cout << X[w] << "\t" << Comp[0].windows[w] << std::endl;
			}
			for (size_t w = 0; w < Win.NumberOfWindows; w++) {
				//std::cout << Y[w] << "\t" << Comp[1].windows[w] << std::endl;
			}
			for (size_t w = 0; w < Win.NumberOfWindows; w++) {
				std::cout << Z[w] << "\t" << Comp[2].windows[w] << std::endl;
			}

			//bookmark
			X = Comp[XCOMP].windows;
			Y = Comp[YCOMP].windows;
			Z = Comp[ZCOMP].windows;

			if (SaveDiagnosticFiles) {
				Win.write_windows("diag_windows.txt", X, Y, Z);
			}
		}

		void initialise_windows(const cBlock& rxblock) {
			Win = WindowingScheme(rxblock, WFM);
		}

		/*
		void initialise_windows()
		{
			NumberOfWindows = (size_t)STM.getintvalue("Receiver.NumberOfWindows");
			if (STM.getvalue("Receiver.TimeShift", WindowsTimeShift) == false) {
				WindowsTimeShift = 0.0;
			}

			WinSpec.resize(NumberOfWindows);
			X.resize(NumberOfWindows);
			Y.resize(NumberOfWindows);
			Z.resize(NumberOfWindows);

			//Read window times
			std::vector<std::vector<double>> wt = STM.getdoublematrix("Receiver.WindowTimes");
			size_t nw = wt.size();
			if (nw != NumberOfWindows) {
				glog.errormsg(_SRC_, "The number of WindowTimes does not match the NumberOfWindows\n");
			}
			for (size_t i = 0; i < NumberOfWindows; i++) {
				if (wt[i].size() != 2) {
					glog.errormsg(_SRC_, "The number of WindowTimes must have exactly 2 columns (error in window %lu)\n", i + 1);
				}
				WinSpec[i].TimeLow = wt[i][0] + WindowsTimeShift;
				WinSpec[i].TimeHigh = wt[i][1] + WindowsTimeShift;
			}

			WindowWeightingScheme = STM.getstringvalue("Receiver.WindowWeightingScheme");
			if (strcasecmp(WindowWeightingScheme, "AreaUnderCurve") == 0) {
				initialise_windows_area();
			}
			else if (strcasecmp(WindowWeightingScheme, "Boxcar") == 0) {
				initialise_windows_boxcar();
			}
			else if (strcasecmp(WindowWeightingScheme, "LinearTaper") == 0) {
				initialise_windows_lineartaper();
			}
			else {
				glog.errormsg(_SRC_, "WindowWeightingScheme %s unknown (must be \"AreaUnderCurve\" or  \"Boxcar\" or \"LinearTaper\")\n", WindowWeightingScheme.c_str());
			}
		}

		void initialise_windows_area()
		{
			double tlow, thigh, t, tp, tn, tleft, tright;

			double dwt = WFM.Time[1] - WFM.Time[0];
			double eps = 1.0e-7;

			for (size_t w = 0; w < NumberOfWindows; w++) {
				WinSpec[w].TimeWidth = WinSpec[w].TimeHigh - WinSpec[w].TimeLow;
				for (size_t s = 0; s < WFM.NumSamples; s++) {
					t = WFM.Time[s];
					if (t + eps >= WinSpec[w].TimeLow) {
						WinSpec[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = WFM.NumSamples; s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					t = WFM.Time[s];
					if (t - eps <= WinSpec[w].TimeHigh) {
						WinSpec[w].SampleHigh = s;
						break;
					}
				}

				WinSpec[w].NumberOfSamples = WinSpec[w].SampleHigh - WinSpec[w].SampleLow + 1;
				WinSpec[w].Sample.resize(WinSpec[w].NumberOfSamples);
				WinSpec[w].Weight.resize(WinSpec[w].NumberOfSamples);
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Sample[k] = WinSpec[w].SampleLow + k;
					WinSpec[w].Weight[k] = 0.0;
				}

				tlow = WinSpec[w].TimeLow;
				thigh = WinSpec[w].TimeHigh;
				double wsum = 0;
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					size_t s = WinSpec[w].Sample[k];
					t = WFM.Time[s];
					tp = WFM.Time[s] - dwt;
					tn = WFM.Time[s] + dwt;
					tleft = std::max(tp, tlow);
					tright = std::min(tn, thigh);
					WinSpec[w].Weight[k] = 0.5 * (t - tleft) + 0.5 * (tright - t);
					WinSpec[w].Weight[k] /= WinSpec[w].TimeWidth;
					wsum += WinSpec[w].Weight[k];
				}
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Weight[k] /= wsum;
				}
			}
		}

		void initialise_windows_boxcar()
		{
			double eps = 1.0e-7;
			for (size_t w = 0; w < NumberOfWindows; w++) {
				WinSpec[w].TimeWidth = WinSpec[w].TimeHigh - WinSpec[w].TimeLow;
				for (size_t s = 0; s < WFM.NumSamples; s++) {
					double t = WFM.Time[s];
					if (t + eps >= WinSpec[w].TimeLow) {
						WinSpec[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = WFM.NumSamples; s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					double t = WFM.Time[s];
					if (t - eps <= WinSpec[w].TimeHigh) {
						WinSpec[w].SampleHigh = s;
						break;
					}
				}


				WinSpec[w].NumberOfSamples = WinSpec[w].SampleHigh - WinSpec[w].SampleLow + 1;
				WinSpec[w].Sample.resize(WinSpec[w].NumberOfSamples);
				WinSpec[w].Weight.resize(WinSpec[w].NumberOfSamples);
				double weightsum = 0.0;
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Sample[k] = WinSpec[w].SampleLow + k;
					WinSpec[w].Weight[k] = 1.0;
					weightsum += WinSpec[w].Weight[k];
				}
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Weight[k] /= weightsum;
				}
				//printf("%lu %lu\n", w, WinSpec[w].NumberOfSamples);
			}
		}

		void initialise_windows_lineartaper()
		{
			double eps = 1.0e-7;
			for (size_t w = 0; w < NumberOfWindows; w++) {
				WinSpec[w].TimeWidth = WinSpec[w].TimeHigh - WinSpec[w].TimeLow;
				for (size_t s = 0; s < WFM.NumSamples; s++) {
					double t = WFM.Time[s];
					if (t + eps >= WinSpec[w].TimeLow) {
						WinSpec[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = WFM.NumSamples; s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					double t = WFM.Time[s];
					if (t - eps <= WinSpec[w].TimeHigh) {
						WinSpec[w].SampleHigh = s;
						break;
					}
				}

				size_t ns = WinSpec[w].SampleHigh - WinSpec[w].SampleLow + 1;
				WinSpec[w].NumberOfSamples = ns * 3;
				WinSpec[w].SampleLow -= ns;
				WinSpec[w].SampleHigh += ns;
				WinSpec[w].Sample.resize(WinSpec[w].NumberOfSamples);
				WinSpec[w].Weight.resize(WinSpec[w].NumberOfSamples);

				double weightsum = 0.0;
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Sample[k] = WinSpec[w].SampleLow + k;
					if (k < ns) {
						WinSpec[w].Weight[k] = (double)(k + 1) / (double)(ns + 1);
					}
					else if (k >= 2 * ns) {
						WinSpec[w].Weight[k] = 1.0 - (double)((k + 1) - 2 * ns) / (double)(ns + 1);
					}
					else {
						WinSpec[w].Weight[k] = 1.0;
					}
					weightsum += WinSpec[w].Weight[k];
				}
				for (size_t k = 0; k < WinSpec[w].NumberOfSamples; k++) {
					WinSpec[w].Weight[k] /= weightsum;
				}
			}
		}

		void computewindow(const double* timeseries, std::vector<double>& W)
		{
			for (size_t w = 0; w < NumberOfWindows; w++) {
				W[w] = 0.0;
				for (size_t k = 0; k < WinSpec[w].Sample.size(); k++) {
					W[w] += timeseries[WinSpec[w].Sample[k]] * WinSpec[w].Weight[k];
				}
			}
		}
		*/

		void write_discretefrequencies(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumberOfDiscreteFrequencies; i++) {
				ofs << strprint("%15le\t%15le\t%15le\t%15le\t%15le\t%15le\t%15le\n", DiscreteFrequencies[i], HxR[i], HxI[i], HyR[i], HyI[i], HzR[i], HzI[i]);
			}
		}

		void write_splinedfrequencies(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumberOfSplinedFrequencies; i++) {
				double f = pow10(SplinedFrequencieslog10[i]);
				ofs << strprint("%15le\t%15le\t%15le\t%15le\t%15le\t%15le\t%15le\n", f, X_splined[i].real(), X_splined[i].imag(), Y_splined[i].real(), Y_splined[i].imag(), Z_splined[i].real(), Z_splined[i].imag());
			}
		}

		void write_timesseries(const std::string& path) const {
			std::ofstream ofs = ofstream_ex(path);
			double* ts = (double*)(WFM.FFTWork.data());
			for (size_t i = 0; i < WFM.NumSamples; i++) {
				ofs << strprint("%20.10le\t%20.10le\n", WFM.Time[i], ts[i]);
			}
		}

		void drx_pitch(double xb, double zb, double p, double& dxbdp, double& dzbdp)
		{
			//xi = (  xb*cosp  + zb*sinp);Inertial
			//zi = ( -xb*sinp  + zb*cosp);
			//xb = (  xi*cosp  - zi*sinp);As bird sees it
			//zb = (  xi*sinp  + zi*cosp);						

			if (Normalisation == NormalizationType::PPM || Normalisation == NormalizationType::PPM_PEAKTOPEAK) {
				//Must work with true field vector directions, not the PPM scaled versinn
				xb *= RefGeomPrimaryX;
				zb *= RefGeomPrimaryZ;
			}

			double cosp = cos(D2R<double> *p);
			double sinp = sin(D2R<double> *p);

			double xi = (xb * cosp + zb * sinp);//convert back to real coordinate system
			double zi = (-xb * sinp + zb * cosp);

			dxbdp = D2R<double> *(-xi * sinp - zi * cosp);
			dzbdp = D2R<double> *(+xi * cosp - zi * sinp);

			if (Normalisation == NormalizationType::PPM || Normalisation == NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dxbdp /= RefGeomPrimaryX;
				dzbdp /= RefGeomPrimaryZ;
			}
		}

		void drx_pitch(std::vector<double> xb, std::vector<double> zb, double p, std::vector<double>& dxbdp, std::vector<double>& dzbdp)
		{
			//xi = (  xb*cosp  + zb*sinp);Inertial
			//zi = ( -xb*sinp  + zb*cosp);
			//xb = (  xi*cosp  - zi*sinp);As bird sees it
			//zb = (  xi*sinp  + zi*cosp);						

			if (Normalisation == NormalizationType::PPM || Normalisation == NormalizationType::PPM_PEAKTOPEAK) {
				//Must work with true field vector directions, not the PPM scaled versinn
				xb *= RefGeomPrimaryX;
				zb *= RefGeomPrimaryZ;
			}


			double cosp = cos(D2R<double> *p);
			double sinp = sin(D2R<double> *p);

			//convert back to real coordinate system
			std::vector<double> xi = (xb * cosp + zb * sinp);
			std::vector<double> zi = (xb * -sinp + zb * cosp);

			dxbdp = (xi * -sinp - zi * cosp) * D2R<double>;
			dzbdp = (xi * cosp - zi * sinp) * D2R<double>;

			if (Normalisation == NormalizationType::PPM || Normalisation == NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dxbdp /= RefGeomPrimaryX;
				dzbdp /= RefGeomPrimaryZ;
			}
		}

		void drx_roll(double yb, double zb, double r, double& dybdr, double& dzbdr)
		{
			//yi = (  yb*cosr  - zb*sinr);Inertial
			//zi = (  yb*sinr  + zb*cosr);
			//yb = (  yi*cosr  + zi*sinr);As bird sees it
			//zb = ( -yi*sinr  + zi*cosr);						

			if (Normalisation == NormalizationType::PPM || Normalisation == NormalizationType::PPM_PEAKTOPEAK) {
				//Must work with true field vector directions, not the PPM scaled versinn
				yb *= RefGeomPrimaryY;
				zb *= RefGeomPrimaryZ;
			}

			double cosr = cos(D2R<double> *r);
			double sinr = sin(D2R<double> *r);

			double yi = (yb * cosr - zb * sinr);//convert back to real coordinate system
			double zi = (yb * sinr + zb * cosr);

			dybdr = D2R<double> *(-yi * sinr + zi * cosr);
			dzbdr = D2R<double> *(-yi * cosr - zi * sinr);

			if (Normalisation == NormalizationType::PPM || Normalisation == NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dybdr /= RefGeomPrimaryY;
				dzbdr /= RefGeomPrimaryZ;
			}
		}

		void  drx_roll(std::vector<double> yb, std::vector<double> zb, double r, std::vector<double>& dybdr, std::vector<double>& dzbdr)
		{
			//yi = (  yb*cosr  - zb*sinr);Inertial
			//zi = (  yb*sinr  + zb*cosr);
			//yb = (  yi*cosr  + zi*sinr);As bird sees it
			//zb = ( -yi*sinr  + zi*cosr);						

			if (Normalisation == NormalizationType::PPM || Normalisation == NormalizationType::PPM_PEAKTOPEAK) {
				//Must work with true field vector directions, not the PPM scaled versinn
				yb *= RefGeomPrimaryY;
				zb *= RefGeomPrimaryZ;
			}

			double cosr = cos(D2R<double> *r);
			double sinr = sin(D2R<double> *r);

			//convert back to real coordinate system
			std::vector<double> yi = (yb * cosr - zb * sinr);
			std::vector<double> zi = (yb * sinr + zb * cosr);

			dybdr = (yi * -sinr + zi * cosr) * D2R<double>;
			dzbdr = (yi * -cosr - zi * sinr) * D2R<double>;

			if (Normalisation == NormalizationType::PPM || Normalisation == NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dybdr /= RefGeomPrimaryY;
				dzbdr /= RefGeomPrimaryZ;
			}
		}

		void setup_scaling() {
			Tx.PeakdIdT = WFM.compute_peak_didt();
			double tx_scale = MUZERO<double> *Tx.LoopArea * Tx.NumberOfTurns * Tx.PeakCurrent;
			double xos = STM.getdoublevalue("ForwardModelling.XOutputScaling");
			double yos = STM.getdoublevalue("ForwardModelling.YOutputScaling");
			double zos = STM.getdoublevalue("ForwardModelling.ZOutputScaling");

			Scale.resize(NCOMP);
			Scale[XCOMP] = tx_scale * xos;
			Scale[YCOMP] = tx_scale * yos;
			Scale[ZCOMP] = tx_scale * zos;

			if (Normalisation == NormalizationType::PPM || Normalisation == NormalizationType::PPM_PEAKTOPEAK) {
				cBlock b = STM.findblock("ReferenceGeometry");
				if (b.Entries.size() == 0) {
					glog.errormsg(_SRC_, "Must define a ReferenceGeometry for PPM or PPMPEAKTOPEAK normalisation\n");
				}
				NormalizationGeometry = cTDEmGeometry(b);
				setgeometry(NormalizationGeometry);
				setprimaryfields();

				double s = 1.0;
				if (Normalisation == NormalizationType::PPM) {
					s *= 1.0e6;
				}
				else if (Normalisation == NormalizationType::PPM_PEAKTOPEAK) {
					s *= 1.0e6;
				}

				RefGeomPrimaryX = PrimaryX;
				RefGeomPrimaryY = PrimaryY;
				RefGeomPrimaryZ = PrimaryZ;

				if (RefGeomPrimaryX == 0.0) Scale[XCOMP] = 0.0;
				else Scale[XCOMP] *= (s / RefGeomPrimaryX);

				if (RefGeomPrimaryY == 0.0) Scale[YCOMP] = 0.0;
				else Scale[YCOMP] *= (s / RefGeomPrimaryY);

				if (RefGeomPrimaryZ == 0.0) Scale[ZCOMP] = 0.0;
				else Scale[ZCOMP] *= (s / RefGeomPrimaryZ);
			}

		}

		void setupdiscretefrequencies()
		{
			double lf1 = log10(WFM.BaseFrequency);
			double lf2 = log10(WFM.SampleFrequency / 2);
			double dlf = 1.0 / FrequenciesPerDecade;

			lf1 = lf1 - 2.0 * dlf;
			lf2 = lf2 + 2.0 * dlf;

			size_t nf = (size_t)ceil((lf2 - lf1) * FrequenciesPerDecade);
			dlf = (lf2 - lf1) / double(nf - 1);

			NumberOfDiscreteFrequencies = nf;
			FrequencyLog10Spacing = dlf;
			DiscreteFrequencyLow = pow(10.0, lf1);
			DiscreteFrequencyHigh = pow(10.0, lf2);

			DiscreteFrequenciesLog10 = std::vector<double>(NumberOfDiscreteFrequencies);
			DiscreteFrequencies = std::vector<double>(NumberOfDiscreteFrequencies);
			for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++) {
				DiscreteFrequenciesLog10[fi] = log10(DiscreteFrequencyLow) + FrequencyLog10Spacing * (double)fi;
				DiscreteFrequencies[fi] = pow(10.0, DiscreteFrequenciesLog10[fi]);
			}
		}

		void forwardmodel(const cTDEmGeometry& G, const cEarth1D& E, cTDEmResponse& R)
		{
			setgeometry(G);
			setearthproperties(E);
			setupcomputations();
			setprimaryfields();
			setsecondaryfields();

			R.PX = PrimaryX;
			R.PY = PrimaryY;
			R.PZ = PrimaryZ;
			R.SX = X;
			R.SY = Y;
			R.SZ = Z;
		}

		void forwardmodel(const std::vector<double>& conductivity, const std::vector<double>& thickness, const cTDEmGeometry& geometry)
		{
			setconductivitythickness(conductivity, thickness);
			setgeometry(geometry);
			LEM.calculation_type = cLEM::CalculationType::FORWARDMODEL;
			LEM.derivative_layer = undefinedvalue<size_t>();

			setupcomputations();
			setprimaryfields();
			setsecondaryfields();
		}

		void forwardmodel(const size_t nlayers, const double* conductivity, const double* thickness, const double* g, double& PX, double& PY, double& PZ, double* SX, double* SY, double* SZ)
		{
			//Order const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw)
			cTDEmGeometry G(g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9]);
			setgeometry(G);
			setconductivitythickness(nlayers, conductivity, thickness);
			setupcomputations();
			setprimaryfields();
			setsecondaryfields();
			getfields(PX, PY, PZ, SX, SY, SZ);
		}

		void getfields(double& PX, double& PY, double& PZ, double* SX, double* SY, double* SZ)
		{
			PX = PrimaryX;
			PY = PrimaryY;
			PZ = PrimaryZ;
			const size_t nw = nwindows();
			for (size_t i = 0; i < nw; i++) SX[i] = X[i];
			for (size_t i = 0; i < nw; i++) SY[i] = Y[i];
			for (size_t i = 0; i < nw; i++) SZ[i] = Z[i];
		}

	};
};


