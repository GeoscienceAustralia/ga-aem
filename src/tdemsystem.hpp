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
#include "lem.hpp"
#include "layeredearthmodeller.hpp"
#include "fixed_point_spline.hpp"
#include "rollpitchyaw.hpp"

namespace AEM {
	using namespace LEM2;
	using CalculationType = CT::CalculationType;
	using CMode = CT::CalculationType::Mode;
	
	inline constexpr size_t XCOMP = 0;
	inline constexpr size_t YCOMP = 1;
	inline constexpr size_t ZCOMP = 2;
	inline constexpr size_t NCOMP = 3;

	inline static Mat3d YPR(const double& roll_degrees, const double& pitch_degrees, const double& yaw_degrees) {
		const Mat3d Rot = yawpitchroll_matrix(roll_degrees * D2R<double>, pitch_degrees * D2R<double>, yaw_degrees * D2R<double>);
		return Rot;
	};

	inline static Mat3d invYPR(const double& roll_degrees, const double& pitch_degrees, const double& yaw_degrees) {
		const Mat3d Rot = yawpitchroll_matrix(roll_degrees * D2R<double>, pitch_degrees * D2R<double>, yaw_degrees * D2R<double>);
		Mat3d RotT = Rot.transpose();
		return RotT;
	};

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
		cTDEmComponent data[NCOMP];

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

	class TDEmGeometry {

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

		TDEmGeometry() {};

		TDEmGeometry(const double& _tx_height, const double& _tx_roll, const double& _tx_pitch, const double& _tx_yaw, const double& _txrx_dx, const double& _txrx_dy, const double& _txrx_dz, const double& _rx_roll, const double& _rx_pitch, const double& _rx_yaw) {
			initialise(_tx_height, _tx_roll, _tx_pitch, _tx_yaw, _txrx_dx, _txrx_dy, _txrx_dz, _rx_roll, _rx_pitch, _rx_yaw);
		}

		TDEmGeometry(const double* g) {
			//const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw)
			initialise(g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9]);
		}

		TDEmGeometry(const std::vector<double> gvector) {
			for (size_t i = 0; i < size(); i++) {
				(*this)[i] = gvector[i];
			}
		}

		TDEmGeometry(const cBlock& b) {
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

		void initialise(const double& _tx_height, const double& _tx_roll, const double& _tx_pitch, const double& _tx_yaw, const double& _txrx_dx, const double& _txrx_dy, const double& _txrx_dz, const double& _rx_roll, const double& _rx_pitch, const double& _rx_yaw) {
			tx_height = _tx_height;
			tx_roll = _tx_roll; tx_pitch = _tx_pitch; tx_yaw = _tx_yaw;
			txrx_dx = _txrx_dx; txrx_dy = _txrx_dy; txrx_dz = _txrx_dz;
			rx_roll = _rx_roll; rx_pitch = _rx_pitch; rx_yaw = _rx_yaw;
		};

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
			return (*(const_cast<TDEmGeometry*>(this)))[index]; // Correctly calls the function above.		
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
			return (*(const_cast<TDEmGeometry*>(this)))[i]; // Correctly calls the function above.
		};

		void set_zero() {
			for (size_t i = 0; i < size(); i++) {
				(*this)[i] = 0.0;
			}
		}

		void fillundefined(const TDEmGeometry& g)
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
			}
			return ElementType::unknown;
		}

		static CMode derivativetype(const size_t& index) {
			switch (index) {
			case 0: return CMode::DH; break;
			case 1: return CMode::NONE; break;
			case 2: return CMode::NONE; break;
			case 3: return CMode::NONE; break;
			case 4: return CMode::DX; break;
			case 5: return CMode::DY; break;
			case 6: return CMode::DZ; break;
			case 7: return CMode::NONE; break;
			case 8: return CMode::NONE; break;
			case 9: return CMode::NONE; break;
			default:
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
			}
			return CMode::NONE;
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

		Vec3d tx_orientation(const Vec3d& tx_reference_orientation) const {
			Mat3d Rot = YPR(tx_roll, tx_pitch, tx_yaw);
			return Rot * tx_reference_orientation;
		}

		Vec3d txrx_separation() const {
			Vec3d v = Vec3d(txrx_dx, txrx_dy, txrx_dz);
			return v;
		}

		inline Mat3d inertial_to_rx_frame_rotation_matrix() const {
			// Mat3d inertial_to_rx_frame_rotation_matrix() const {
			return invYPR(rx_roll, rx_pitch, rx_yaw);
		};

		inline Mat3d rx_roll_derivative_matrix() const {
			const Mat3d R = roll_matrix_degrees(rx_roll);
			const Mat3d P = pitch_matrix_degrees(rx_pitch);
			const Mat3d Y = yaw_matrix_degrees(rx_yaw);
			const Mat3d dR = roll_matrix_derivative_degrees(rx_roll);
			const Mat3d m = (dR * P * Y).transpose() * R * P * Y;
			return m;
		};

		inline Mat3d rx_pitch_derivative_matrix() const {
			// f(p) = [R(r) P(p) Y(y)]' * I     // fp is field vector in Rx reference frame
			// I    = [R(r) P(p) Y(y)]  * f(p)  //  I is field vector in the inertial reference frame
			// df(p)/dp = d([R(r) P(p) Y(y)]')/dp * I + [R(r) P(p) Y(y)]' * d(I)/dp
			// df(p)/dp = d([R(r) P(p) Y(y)]')/dp * I         since d(I)/dp = 0
			//          = d([R(r)    P(p)  Y(y)]')/dp * [R(r) P(p) Y(y)] * f(p)
			//          =   [R(r) dP(p)/dp Y(y)]'   * [R(r) P(p) Y(y)] * f(p)
			const Mat3d R = roll_matrix_degrees(rx_roll);
			const Mat3d P = pitch_matrix_degrees(rx_pitch);
			const Mat3d Y = yaw_matrix_degrees(rx_yaw);
			const Mat3d dP = pitch_matrix_derivative_degrees(rx_pitch);
			const Mat3d m = (R * dP * Y).transpose() * R * P * Y;
			return m;
		};

		inline Mat3d rx_yaw_derivative_matrix() const {
			const Mat3d R = roll_matrix_degrees(rx_roll);
			const Mat3d P = pitch_matrix_degrees(rx_pitch);
			const Mat3d Y = yaw_matrix_degrees(rx_yaw);
			const Mat3d dY = yaw_matrix_derivative_degrees(rx_yaw);
			const Mat3d m = (R * P * dY).transpose() * R * P * Y;
			return m;
		};
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
		Vec3d Reference_Orientation = Vec3d::UnitZ();
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

		Waveform::Type Type; // Is the time domain waveform specified as TX current or Rx voltage
		std::vector<double> Time; // Times in seconds
		std::vector<double>  TD_Waveform; // Time domain waveform
		std::vector<cdouble> FD_Waveform; // Pure frequency domain waveform
		std::vector<cdouble> TransferFunction; // FD_Waveform * RX Filters * (b->db/dt or db/dt->b) conversion
		std::vector<cdouble> FFT_WorkArray; // Work array for repeated inplace inverse FFTs
		std::vector<double>  FFT_Frequency; // Pre-computed FFT frequencies

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
						digitisewaveform(wp, Time, TD_Waveform);
						Type = Waveform::Type::RX;
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
						digitisewaveform(wp, Time, TD_Waveform);
						Type = Waveform::Type::TX;
						wavformdefined = true;
					}
				}
			}

			if (wavformdefined == false) {
				std::vector<std::vector<double>> wp = b.getdoublematrix("WaveformCurrent");
				if (wp.size() > 0) {
					digitisewaveform(wp, Time, TD_Waveform);
					Type = Waveform::Type::TX;
					wavformdefined = true;
				}
			}

			if (wavformdefined == false) {
				std::vector<std::vector<double>> wp = b.getdoublematrix("WaveformReceived");
				if (wp.size() > 0) {
					digitisewaveform(wp, Time, TD_Waveform);
					Type = Waveform::Type::RX;
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
			for (size_t i = 1; i < TD_Waveform.size(); i++) {
				const double di = TD_Waveform[i] - TD_Waveform[i - 1];
				const double dt = Time[i] - Time[i - 1];
				double didt = std::fabs(di / dt);
				if (didt > maxdidt) maxdidt = didt;
			}
			return maxdidt;
		}

		void write_timedomainwaveform(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumSamples; i++) {
				ofs << strprint("%20le\t%20le\n", Time[i], TD_Waveform[i]);
			}
		}

		void write_frequencydomainwaveform(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumFrequencies; i++) {
				ofs << strprint("%15le\t%15le\t%15le\t%15le\t%15le\n", FFT_Frequency[i], FD_Waveform[i].real(), FD_Waveform[i].imag(), TransferFunction[i].real(), TransferFunction[i].imag());
			}
		}

		void write_frequencyseries(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumFrequencies; i++) {
				ofs << strprint("%15le\t%15le\t%15le\n", FFT_Frequency[i], FFT_WorkArray[i].real(), FFT_WorkArray[i].imag());
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

		double centre_time() const {
			return (TimeLow + TimeHigh) / 2.0;
		}
	};

	class WindowingScheme {

	private:


	public:

		enum class WeightingMethod { BoxCar, AreaUnderCurve, LinearTaper };

		size_t nWindows = 0;
		double TimeShift = 0.0;
		WeightingMethod Method = WeightingMethod::BoxCar;
		std::vector<WindowSpecification> Windows;

		WindowingScheme() {};

		WindowingScheme(const cBlock& receiverblock, const Waveform& WFM) {
			const cBlock& b = receiverblock;

			if (!b.getvalue("NumberOfWindows", nWindows)) {
				glog.errormsg(_SRC_, "NumberOfWindows is not specified");
			}

			if (!b.getvalue("TimeShift", TimeShift)) {
				TimeShift = 0.0;
			}

			Windows.resize(nWindows);

			//Read window times
			std::vector<std::vector<double>> wt;
			if (!b.getvalue("WindowTimes", wt)) {
				glog.errormsg(_SRC_, "The WindowTimes have not been specified.");
			}
			size_t nw = wt.size();
			if (nw != nWindows) {
				glog.errormsg(_SRC_, "The number of WindowTimes does not match the NumberOfWindows\n");
			}

			for (size_t i = 0; i < nWindows; i++) {
				if (wt[i].size() != 2) {
					glog.errormsg(_SRC_, "The number of WindowTimes must have exactly 2 columns (error in window %lu)\n", i + 1);
				}
				Windows[i].TimeLow = wt[i][0] + TimeShift;
				Windows[i].TimeHigh = wt[i][1] + TimeShift;
			}

			std::string wmethod = b.getstringvalue("WindowWeightingScheme");
			if (strcasecmp(wmethod, "AreaUnderCurve") == 0) {
				initialise_area(WFM);
			}
			else if (strcasecmp(wmethod, "Boxcar") == 0) {
				initialise_boxcar(WFM);
			}
			else if (strcasecmp(wmethod, "LinearTaper") == 0) {
				initialise_lineartaper(WFM);
			}
			else glog.errormsg(_SRC_, "WindowWeightingScheme %s unknown (must be \"AreaUnderCurve\" or  \"Boxcar\" or \"LinearTaper\")\n", wmethod.c_str());
		}

		const size_t& nwindows() const { return nWindows; }

		void initialise_area(const Waveform& WFM) {
			double tlow, thigh, t, tp, tn, tleft, tright;

			double dwt = WFM.Time[1] - WFM.Time[0];
			double eps = 1.0e-7;

			for (size_t w = 0; w < nWindows; w++) {
				Windows[w].TimeWidth = Windows[w].TimeHigh - Windows[w].TimeLow;
				for (size_t s = 0; s < WFM.NumSamples; s++) {
					t = WFM.Time[s];
					if (t + eps >= Windows[w].TimeLow) {
						Windows[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = WFM.NumSamples; s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					t = WFM.Time[s];
					if (t - eps <= Windows[w].TimeHigh) {
						Windows[w].SampleHigh = s;
						break;
					}
				}

				Windows[w].NumberOfSamples = Windows[w].SampleHigh - Windows[w].SampleLow + 1;
				Windows[w].Sample.resize(Windows[w].NumberOfSamples);
				Windows[w].Weight.resize(Windows[w].NumberOfSamples);
				for (size_t k = 0; k < Windows[w].NumberOfSamples; k++) {
					Windows[w].Sample[k] = Windows[w].SampleLow + k;
					Windows[w].Weight[k] = 0.0;
				}

				tlow = Windows[w].TimeLow;
				thigh = Windows[w].TimeHigh;
				double wsum = 0;
				for (size_t k = 0; k < Windows[w].NumberOfSamples; k++) {
					size_t s = Windows[w].Sample[k];
					t = WFM.Time[s];
					tp = WFM.Time[s] - dwt;
					tn = WFM.Time[s] + dwt;
					tleft = std::max(tp, tlow);
					tright = std::min(tn, thigh);
					Windows[w].Weight[k] = 0.5 * (t - tleft) + 0.5 * (tright - t);
					Windows[w].Weight[k] /= Windows[w].TimeWidth;
					wsum += Windows[w].Weight[k];
				}
				for (size_t k = 0; k < Windows[w].NumberOfSamples; k++) {
					Windows[w].Weight[k] /= wsum;
				}
			}
		}

		void initialise_boxcar(const Waveform& WFM)
		{
			double eps = 1.0e-7;
			for (size_t w = 0; w < nWindows; w++) {
				Windows[w].TimeWidth = Windows[w].TimeHigh - Windows[w].TimeLow;
				for (size_t s = 0; s < WFM.NumSamples; s++) {
					double t = WFM.Time[s];
					if (t + eps >= Windows[w].TimeLow) {
						Windows[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = WFM.NumSamples; s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					double t = WFM.Time[s];
					if (t - eps <= Windows[w].TimeHigh) {
						Windows[w].SampleHigh = s;
						break;
					}
				}


				Windows[w].NumberOfSamples = Windows[w].SampleHigh - Windows[w].SampleLow + 1;
				Windows[w].Sample.resize(Windows[w].NumberOfSamples);
				Windows[w].Weight.resize(Windows[w].NumberOfSamples);
				double weightsum = 0.0;
				for (size_t k = 0; k < Windows[w].NumberOfSamples; k++) {
					Windows[w].Sample[k] = Windows[w].SampleLow + k;
					Windows[w].Weight[k] = 1.0;
					weightsum += Windows[w].Weight[k];
				}
				for (size_t k = 0; k < Windows[w].NumberOfSamples; k++) {
					Windows[w].Weight[k] /= weightsum;
				}
				//printf("%lu %lu\n", w, WinSpec[w].NumberOfSamples);
			}
		}

		void initialise_lineartaper(const Waveform& WFM)
		{
			double eps = 1.0e-7;
			for (size_t w = 0; w < nWindows; w++) {
				Windows[w].TimeWidth = Windows[w].TimeHigh - Windows[w].TimeLow;
				for (size_t s = 0; s < WFM.NumSamples; s++) {
					double t = WFM.Time[s];
					if (t + eps >= Windows[w].TimeLow) {
						Windows[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = WFM.NumSamples; s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					double t = WFM.Time[s];
					if (t - eps <= Windows[w].TimeHigh) {
						Windows[w].SampleHigh = s;
						break;
					}
				}

				size_t ns = Windows[w].SampleHigh - Windows[w].SampleLow + 1;
				Windows[w].NumberOfSamples = ns * 3;
				Windows[w].SampleLow -= ns;
				Windows[w].SampleHigh += ns;
				Windows[w].Sample.resize(Windows[w].NumberOfSamples);
				Windows[w].Weight.resize(Windows[w].NumberOfSamples);

				double weightsum = 0.0;
				for (size_t k = 0; k < Windows[w].NumberOfSamples; k++) {
					Windows[w].Sample[k] = Windows[w].SampleLow + k;
					if (k < ns) {
						Windows[w].Weight[k] = (double)(k + 1) / (double)(ns + 1);
					}
					else if (k >= 2 * ns) {
						Windows[w].Weight[k] = 1.0 - (double)((k + 1) - 2 * ns) / (double)(ns + 1);
					}
					else {
						Windows[w].Weight[k] = 1.0;
					}
					weightsum += Windows[w].Weight[k];
				}
				for (size_t k = 0; k < Windows[w].NumberOfSamples; k++) {
					Windows[w].Weight[k] /= weightsum;
				}
			}
		}

		void computewindow(const double* timeseries, std::vector<double>& windowed_values) {
			std::fill(windowed_values.begin(), windowed_values.end(), 0.0); // Reset to zero
			for (size_t w = 0; w < nWindows; w++) {
				for (size_t k = 0; k < Windows[w].Sample.size(); k++) {
					windowed_values[w] += timeseries[Windows[w].Sample[k]] * Windows[w].Weight[k];
				}
			}
		}

		void printwindows(const double& PX, const double& PY, const double& PZ, const std::vector<double>& SX, const std::vector<double>& SY, const std::vector<double>& SZ) const {
			printf("Primary   %15.8lf%15.8lf%15.8lf\n\n", PX, PY, PZ);
			printf("Window#             X               Y               Z\n");
			for (size_t w = 0; w < nWindows; w++) {
				printf("%2zu        %15.8lf%15.8lf%15.8lf\n", w + 1, SX[w], SY[w], SZ[w]);
			}
		};

		void write_windows(const fs::path& path, const std::vector<double>& SX, const std::vector<double>& SY, const std::vector<double>& SZ) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t w = 0; w < nWindows; w++) {
				ofs << strprint("%2zu\t%20e\t%20e\t%15e%15e%15e\n", w + 1, Windows[w].TimeLow, Windows[w].TimeHigh, SX[w], SY[w], SZ[w]);
			}
		};
	};

	class ComponentWorkStore {

	public:

		std::vector<double> IR_discrete_real;// Real impulse response discrete frequency nodes
		std::vector<double> IR_discrete_imag;// Imaginary impulse response discrete frequency nodes
		std::vector<cdouble> IR_splined;// Complex splines impulse response

		double Scale = 0.0;
		double Primary = 0.0;
		double RefGeomPrimary = 0.0;;
		std::vector<double> Secondary;
		//double Scale = 0.0;

		ComponentWorkStore() {};

		void resize(const size_t nnodes, const size_t nfftfreq, const size_t nwindows) {
			IR_discrete_real.resize(nnodes);
			IR_discrete_imag.resize(nnodes);
			IR_splined.resize(nfftfreq);
			Secondary.resize(nwindows);
		};
	};

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

	class cTDEmSystem : public AEMSystem {

	private:

		FFTWPlanWrapper InverseFFTPlan;

		std::vector<ComponentWorkStore> Component;
		FixedPointSpline<double> FrequencySpliner;
				
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

		TDEmGeometry Geometry;
		TDEmGeometry NormalizationGeometry;
		Mat3d RotMatrixToRxFrame;

	public:

		const WindowSpecification& window(const size_t w) const { return WindScheme.Windows[w]; }
		const Waveform& waveform() const { return WvForm; }
		const Transmitter& transmitter() const { return Tx; }

		ModellingOptions MO;

		cTDEmSystem(const fs::path& descriptorpath) {
			read_system_descriptor_file(descriptorpath);
		};

		cTDEmSystem() {};
		//cTDEmSystem(const cTDEmSystem& other) = delete;
		//cTDEmSystem& operator=(const cTDEmSystem& other) = delete;

		const size_t& nwindows() const {
			return WindScheme.nwindows();
		}

		const double& PX() const { return Component[XCOMP].Primary; };
		const double& PY() const { return Component[YCOMP].Primary; };
		const double& PZ() const { return Component[ZCOMP].Primary; };

		const std::vector<double>& XS() const { return Component[XCOMP].Secondary; }
		const std::vector<double>& YS() const { return Component[YCOMP].Secondary; }
		const std::vector<double>& ZS() const { return Component[ZCOMP].Secondary; }

		double primary(const size_t component) const {
			assert(component < NCOMP);
			return Component[component].Primary;
		}

		double secondary(const size_t component, const size_t window) {
			assert(component < NCOMP);
			return Component[component].Secondary[window];
		}

		std::vector<double> secondary(const size_t component) const {
			assert(component < NCOMP);
			return Component[component].Secondary;
		}

	public:

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
		void set_earth(const Earth1D& earth) {
			lem().set_earth(earth);
		};

		void set_geometry(const TDEmGeometry& G) {
			Geometry = G;

			// Set geometry inside the LE Modeller
			Vec3d tx_reference_orientation = Vec3d::UnitZ();
			const Vec3d sep = Geometry.txrx_separation();
			const double& h = Geometry.tx_height;
			const double& x = sep.x();
			const double& y = sep.y();
			const double& z = h + sep.z();
			const Vec3d tx_orientation = Geometry.tx_orientation(Tx.Reference_Orientation);
			lem().set_geometry(tx_orientation, h, x, y, z);
			// Set the rotation matrix for rotating vector fields to Rx frame of reference
			RotMatrixToRxFrame = Geometry.inertial_to_rx_frame_rotation_matrix();
		};

		void set_calculationtype(const CalculationType& _calculationtype) {
			LEM.set_calculationtype(_calculationtype);
		};

		void set_response(cTDEmResponse& Response) const {
			Response.PX = PX();
			Response.PY = PY();
			Response.PZ = PZ();
			Response.SX = XS();
			Response.SY = YS();
			Response.SZ = ZS();
		};

		void forwardmodel(const TDEmGeometry& G, const Earth1D& E, cTDEmResponse& R) {
			set_geometry(G);
			set_earth(E);
			setup_computations();
			setprimaryfields();
			setsecondaryfields();
			set_response(R);
		}

		void getfields(double& px, double& py, double& pz, double* sx, double* sy, double* sz) const {
			px = PX();
			py = PY();
			pz = PZ();
			const size_t nw = nwindows();
			for (size_t i = 0; i < nw; i++) sx[i] = XS()[i];
			for (size_t i = 0; i < nw; i++) sy[i] = YS()[i];
			for (size_t i = 0; i < nw; i++) sz[i] = ZS()[i];
		}

		void setup_computations() {
			lem().setup_computations();
		}

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

			Component[XCOMP].Primary = v.x() * Component[XCOMP].Scale;
			Component[YCOMP].Primary = v.y() * Component[YCOMP].Scale;
			Component[ZCOMP].Primary = v.z() * Component[ZCOMP].Scale;
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

				for(size_t ci = 0; ci < NCOMP; ci++){
					Component[ci].IR_discrete_real[fi] = v[ci].real();
					Component[ci].IR_discrete_imag[fi] = v[ci].imag();
				}
			};

			//Spline discreet frequencies		
			for (size_t i = 0; i < NCOMP; i++) {
				if (Component[i].Scale == 0.0) return;
				spline_component(i);
				inverse_fft_window_scale_component(i);
			}

			if (MO.SaveDiagnosticFiles) {
				write_discretefrequencies("diag_discretefrequencies.txt");
				write_splinedfrequencies("diag_splinedfrequencies.txt");
				WvForm.write_frequencydomainwaveform("diag_frequencydomainwaveform.txt");
			}

			if (MO.SaveDiagnosticFiles) {
				WindScheme.write_windows("diag_windows.txt", XS(), YS(), ZS());
			}
		}

		void drx_pitch(double xb, double zb, double p, double& dxbdp, double& dzbdp) {
			//xi = (  xb*cosp  + zb*sinp);Inertial
			//zi = ( -xb*sinp  + zb*cosp);
			//xb = (  xi*cosp  - zi*sinp);As bird sees it
			//zb = (  xi*sinp  + zi*cosp);						

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM || 
				MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				//Must work with true field vector directions, not the PPM scaled versinn
				xb *= Component[XCOMP].RefGeomPrimary;
				zb *= Component[ZCOMP].RefGeomPrimary;
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
				dxbdp /= Component[XCOMP].RefGeomPrimary;
				dzbdp /= Component[ZCOMP].RefGeomPrimary;
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
				xb *= Component[XCOMP].RefGeomPrimary;
				zb *= Component[ZCOMP].RefGeomPrimary;
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
				dxbdp /= Component[XCOMP].RefGeomPrimary;
				dzbdp /= Component[ZCOMP].RefGeomPrimary;
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
				yb *= Component[YCOMP].RefGeomPrimary;
				zb *= Component[ZCOMP].RefGeomPrimary;
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
				dybdr /= Component[YCOMP].RefGeomPrimary;
				dzbdr /= Component[ZCOMP].RefGeomPrimary;
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
				yb *= Component[YCOMP].RefGeomPrimary;
				zb *= Component[ZCOMP].RefGeomPrimary;
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
				dybdr /= Component[YCOMP].RefGeomPrimary;
				dzbdr /= Component[ZCOMP].RefGeomPrimary;
			}
		}

		void drx_roll_new(const TDEmGeometry& g, const std::vector<Vec3d>& fields, std::vector<Vec3d>& derivatives) const {
			Mat3d dM = g.rx_roll_derivative_matrix();
			apply_rx_derivative_matrix(dM, fields, derivatives);
		};

		void drx_pitch_new(const TDEmGeometry& g, const std::vector<Vec3d>& fields, std::vector<Vec3d>& derivatives) const {
			Mat3d dM = g.rx_pitch_derivative_matrix();
			apply_rx_derivative_matrix(dM, fields, derivatives);
		};

		void drx_yaw_new(const TDEmGeometry& g, const std::vector<Vec3d>& fields, std::vector<Vec3d>& derivatives) const {
			Mat3d dM = g.rx_yaw_derivative_matrix();
			apply_rx_derivative_matrix(dM, fields, derivatives);
		};

	private:

		void apply_rx_derivative_matrix(const Mat3d& dM, const std::vector<Vec3d>& fields, std::vector<Vec3d>& derivatives) const {
			const size_t n = fields.size();
			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM || MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				for (size_t i = 0; i < n; i++) {
					Vec3d ftrue = fields[i];
					//Must work with true field vector directions (not the PPM scaled versinn)
					ftrue[XCOMP] *= Component[XCOMP].RefGeomPrimary;
					ftrue[YCOMP] *= Component[YCOMP].RefGeomPrimary;
					ftrue[ZCOMP] *= Component[ZCOMP].RefGeomPrimary;
					derivatives[i] = dM * ftrue;
					//Convert back to PPMS
					derivatives[i][XCOMP] /= Component[XCOMP].RefGeomPrimary;
					derivatives[i][YCOMP] /= Component[YCOMP].RefGeomPrimary;
					derivatives[i][ZCOMP] /= Component[ZCOMP].RefGeomPrimary;
				}
			}
			else {
				for (size_t i = 0; i < n; i++) {
					derivatives[i] = dM * fields[i];
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
			// Make sure every odd harmonic to exactly zero
			//for (size_t k = 0; k < WvForm.NumFrequencies; k += 2) {
			//	WvForm.F_Waveform[k] = 0.0;
			//}

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
			Component[XCOMP].Scale = tx_scale * MO.XOutputScaling;
			Component[YCOMP].Scale = tx_scale * MO.YOutputScaling;
			Component[ZCOMP].Scale = tx_scale * MO.ZOutputScaling;

			if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM || MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
				cBlock b = STM.findblock("ReferenceGeometry");
				if (b.Entries.size() == 0) {
					glog.errormsg(_SRC_, "Must define a ReferenceGeometry for PPM or PPMPEAKTOPEAK normalisation\n");
				}
				NormalizationGeometry = TDEmGeometry(b);
				set_geometry(NormalizationGeometry);
				setprimaryfields();

				double s = 1.0;
				if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM) {
					s *= 1.0e6;
				}
				else if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
					s *= 1.0e6;
				}

				for (size_t i = 0; i < NCOMP; i++) {
					ComponentWorkStore& c = Component[i];
					c.RefGeomPrimary = c.Primary;
					if (c.RefGeomPrimary == 0.0) c.Scale = 0.0;
					else c.Scale *= (s / c.RefGeomPrimary);
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
			Component[XCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies, nwindows());
			Component[YCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies, nwindows());
			Component[ZCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies, nwindows());
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
			WindScheme.computewindow((double*)WvForm.FFT_WorkArray.data(), C.Secondary);
			if (MO.SaveDiagnosticFiles) {
				write_timesseries("diag_xtimeseries.txt");
			}

			// Scale
			C.Secondary *= Component[component].Scale;
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


