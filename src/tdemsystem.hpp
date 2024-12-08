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
#include "fixed_point_spline.hpp"
#include "rollpitchyaw.hpp"

namespace AEM {

	inline static Mat3 YPR(const double& roll_degrees, const double& pitch_degrees, const double& yaw_degrees) {
		const Mat3 Rot = yawpitchroll_matrix(roll_degrees * D2R<double>, pitch_degrees * D2R<double>, yaw_degrees * D2R<double>);
		return Rot;
	};

	inline static Mat3 invYPR(const double& roll_degrees, const double& pitch_degrees, const double& yaw_degrees) {
		const Mat3 Rot = yawpitchroll_matrix(roll_degrees * D2R<double>, pitch_degrees * D2R<double>, yaw_degrees * D2R<double>);
		//std::cout << Rot << std::endl;
		Mat3 RotT = Rot.transpose();
		//std::cout << RotT << std::endl;
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
			initialise(_tx_height, _tx_roll, _tx_pitch, _tx_yaw, _txrx_dx, _txrx_dy, _txrx_dz, _rx_roll, _rx_pitch, _rx_yaw);
		}

		cTDEmGeometry(const double* g) {
			//const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw)
			initialise(g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9]);
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

		Vec3 tx_orientation(const Vec3& tx_reference_orientation) const {
			//Vec3 v(1, 2, 3);
			//Mat3 Rot = YPR(tx_roll, tx_pitch, tx_yaw);
			//Mat3 invRot = invYPR(tx_roll, tx_pitch, tx_yaw);
			//Vec3 v1 = Rot * v;
			//Vec3 v2 = invRot * v1;
			//std::cout << v << std::endl << std::endl;
			//std::cout << v1 << std::endl << std::endl;
			//std::cout << v2 << std::endl << std::endl;

			//Vec3 v = tx_reference_orientation;
			//v.rotate_inplace(tx_yaw, Geometry3D::zaxis);
			//v.rotate_inplace(tx_pitch, Geometry3D::yaxis);
			//v.rotate_inplace(tx_roll, Geometry3D::xaxis);
			Mat3 Rot = YPR(tx_roll, tx_pitch, tx_yaw);
			return Rot * tx_reference_orientation;
		}

		Vec3 txrx_separation() const {
			Vec3 v = Vec3(txrx_dx, txrx_dy, txrx_dz);
			return v;
		}

		inline Mat3 inertial_to_rx_frame_rotation_matrix() const {
			// Mat3 inertial_to_rx_frame_rotation_matrix() const {
			return invYPR(rx_roll, rx_pitch, rx_yaw);
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
		Vec3 Reference_Orientation = Vec3::UnitZ();
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

		double centre_time() const {
			return (TimeLow + TimeHigh) / 2.0;
		}
	};

	class WindowingScheme {

	private:


	public:

		enum class WeightingMethod { BoxCar, AreaUnderCurve, LinearTaper};

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

		std::vector<double> computewindow(const double* timeseries)
		{
			std::vector<double> W(nWindows, 0.0);
			for (size_t w = 0; w < nWindows; w++) {
				for (size_t k = 0; k < Windows[w].Sample.size(); k++) {
					W[w] += timeseries[Windows[w].Sample[k]] * Windows[w].Weight[k];
				}
			}
			return W;
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

		double Primary = 0.0;
		std::vector<double> Secondary;

		ComponentWorkStore() {};
		void resize(const size_t nnodes, const size_t nfftfreq) {
			IR_discrete_real.resize(nnodes);
			IR_discrete_imag.resize(nnodes);
			IR_splined.resize(nfftfreq);
		};
	};

	class AEMSystem {

	protected:
		std::string SystemName;
		std::string SystemType;
		cBlock STM;
		cLEM LEM;

		bool SaveDiagnosticFiles = false;

	public:

		inline static const size_t XCOMP = 0;
		inline static const size_t YCOMP = 1;
		inline static const size_t ZCOMP = 2;
		inline static const size_t NCOMP = 3;

		const cBlock& system_descriptor_block() const { return STM; };
		cLEM& lem() { return LEM; };

	};
	
	class cTDEmSystem : public AEMSystem {

	private:
		enum class OutputType { BFIELD, DBDT };
		enum class NormalizationType { NONE, PPM, PPM_PEAKTOPEAK };
		OutputType OutputType;
		NormalizationType NormalisationType;

		fftw_plan inverse_fftplan;
		std::vector<ComponentWorkStore> Comp;
		FixedPointSpline<double> FrequencySpliner;

		size_t FrequenciesPerDecade = 0;
		double FrequencyLog10Spacing = 0.0;
		double DiscreteFrequencyLow = 0.0;
		double DiscreteFrequencyHigh = 0.0;
		size_t NumberOfDiscreteFrequencies = 0;
		std::vector<double> DiscreteFrequencies;
		std::vector<double> DiscreteFrequenciesLog10;
		
		size_t NumberOfSplinedFrequencies = 0;
		std::vector<double> SplinedFrequencieslog10;
		
		std::vector<LowPassFilter> Filters;

		cTDEmGeometry Geometry;
		cTDEmGeometry NormalizationGeometry;

		std::vector<double> Scale;
		double RefGeomPrimaryX = 0.0;  //Primary X ref field for PPM normalisation
		double RefGeomPrimaryY = 0.0;  //Primary Y ref field for PPM normalisation
		double RefGeomPrimaryZ = 0.0;  //Primary Z ref field for PPM normalisation
		
		WindowingScheme WindScheme;
		Waveform WvForm;
		Transmitter Tx;

	public:
		
		const WindowSpecification& window(const size_t w) const { return WindScheme.Windows[w]; }
		const Waveform& waveform() const { return WvForm; }
		const Transmitter& transmitter() const { return Tx; }
		
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
			return WindScheme.nwindows();
		}

		const double& PX() const { return Comp[XCOMP].Primary; };
		const double& PY() const { return Comp[YCOMP].Primary; };
		const double& PZ() const { return Comp[ZCOMP].Primary; };
		
		const std::vector<double>& XS() const { return Comp[XCOMP].Secondary; }
		const std::vector<double>& YS() const { return Comp[YCOMP].Secondary; }
		const std::vector<double>& ZS() const { return Comp[ZCOMP].Secondary; }

		double primary(const size_t component) const {
			assert(component < NCOMP);
			return Comp[component].Primary;
		}

		double secondary(const size_t component, const size_t window) {
			assert(component < NCOMP);
			return Comp[component].Secondary[window];
		}

		std::vector<double> secondary(const size_t component) const {
			assert(component < NCOMP);
			return Comp[component].Secondary;
		}

		void initialise() {
			LEM.calculation_type = cLEM::CalculationType::FORWARDMODEL;
			LEM.rzerotype = cLEM::RZeroMethod::PROPOGATIONMATRIX;
			inverse_fftplan = 0;
		}

		// Setup
		void read_system_descriptor_file(const std::string& systemdescriptorfile) {
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

			WvForm.initialise(b, systemdescriptorfile);

			cBlock rxblock = STM.findblock("Receiver");
			WindScheme = WindowingScheme(rxblock, WvForm);

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

			SaveDiagnosticFiles = STM.getboolvalue("ForwardModelling.SaveDiagnosticFiles");

			if (WvForm.Time.size() <= 2 || WvForm.Time.size() != WvForm.T_Waveform.size()) {
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

			system_setup();
		};

		void system_setup() {
			setup_discrete_frequencies();
			setup_transforms();
			setup_splines();
			setup_scaling();
		}

		void setup_transforms() {
			WvForm.NumFrequencies = WvForm.NumSamples / 2 + 1;
			WvForm.fft_frequency.resize(WvForm.NumFrequencies);

			size_t N = WvForm.NumSamples;
			size_t NC = N;
			size_t NR = 2 * (N / 2 + 1);

			//Forward transform
			WvForm.F_Waveform.resize(NC);
			double* in = (double*)&(WvForm.T_Waveform[0]);
			fftw_complex* out = (fftw_complex*)&(WvForm.F_Waveform[0]);
			fftw_plan fftwplan_forward = fftw_plan_dft_r2c_1d((int)N, in, out, FFTW_ESTIMATE);
			for (size_t k = 0; k < WvForm.NumSamples; k++) {
				WvForm.T_Waveform[k] /= (double)WvForm.NumSamples;
			}
			fftw_execute(fftwplan_forward);
			fftw_destroy_plan(fftwplan_forward);

			for (size_t k = 0; k < WvForm.NumFrequencies; k++) {
				WvForm.fft_frequency[k] = WvForm.calculate_fft_frequency(k);
			}

			bool convert_B_2_dBdT = false;
			bool convert_dBdT_2_B = false;
			if (WvForm.Type == Waveform::Type::TX) {
				if (OutputType == OutputType::DBDT) {
					convert_B_2_dBdT = true;
				}
			}

			if (WvForm.Type == Waveform::Type::RX) {
				if (OutputType == OutputType::BFIELD) {
					convert_dBdT_2_B = true;
				}
			}

			for (size_t k = 0; k < WvForm.NumFrequencies; k += 2) {
				WvForm.F_Waveform[k] = 0.0;
			}

			WvForm.Transfer.resize(WvForm.NumFrequencies);
			WvForm.Transfer = WvForm.F_Waveform;

			for (size_t k = 0; k < WvForm.NumFrequencies; k++) {
				const double& frequency = WvForm.fft_frequency[k];
				if (convert_B_2_dBdT == true) {
					WvForm.Transfer[k] *= cdouble(0.0, -TWOPI<double> *frequency);
				}
				if (convert_dBdT_2_B == true) {
					WvForm.Transfer[k] *= cdouble(0.0, -1.0 / (TWOPI<double> *frequency));
				}

				for (size_t fi = 0; fi < Filters.size(); fi++) {
					const cdouble w = Filters[fi].weight(frequency);
					WvForm.Transfer[k] *= w;
				}
			}

			//Frequencies to be splined
			NumberOfSplinedFrequencies = WvForm.NumFrequencies / 2;
			SplinedFrequencieslog10.resize(NumberOfSplinedFrequencies);
			for (size_t k = 0; k < NumberOfSplinedFrequencies; k++) {
				SplinedFrequencieslog10[k] = log10(fabs(WvForm.fft_frequency[k * 2 + 1]));
			}

			//Setup inverse transform work array	
			WvForm.FFTWork.resize(NR);

#if defined MULTITHREADED
			//FFTW_MEASURE does not seem to be thread safe
			unsigned int FLAGS = FFTW_ESTIMATE;
#else
			unsigned int FLAGS = FFTW_MEASURE;
#endif

			fftw_complex* invin = (fftw_complex*)(WvForm.FFTWork.data());
			double* invout = (double*)(WvForm.FFTWork.data());
			inverse_fftplan = fftw_plan_dft_c2r_1d((int)N, invin, invout, FLAGS);
		}

		void setup_scaling() {
			Tx.PeakdIdT = WvForm.compute_peak_didt();
			double tx_scale = MUZERO<double> *Tx.LoopArea * Tx.NumberOfTurns * Tx.PeakCurrent;
			double xos = STM.getdoublevalue("ForwardModelling.XOutputScaling");
			double yos = STM.getdoublevalue("ForwardModelling.YOutputScaling");
			double zos = STM.getdoublevalue("ForwardModelling.ZOutputScaling");

			Scale.resize(NCOMP);
			Scale[XCOMP] = tx_scale * xos;
			Scale[YCOMP] = tx_scale * yos;
			Scale[ZCOMP] = tx_scale * zos;

			if (NormalisationType == NormalizationType::PPM || NormalisationType == NormalizationType::PPM_PEAKTOPEAK) {
				cBlock b = STM.findblock("ReferenceGeometry");
				if (b.Entries.size() == 0) {
					glog.errormsg(_SRC_, "Must define a ReferenceGeometry for PPM or PPMPEAKTOPEAK normalisation\n");
				}
				NormalizationGeometry = cTDEmGeometry(b);
				setgeometry(NormalizationGeometry);
				setprimaryfields();

				double s = 1.0;
				if (NormalisationType == NormalizationType::PPM) {
					s *= 1.0e6;
				}
				else if (NormalisationType == NormalizationType::PPM_PEAKTOPEAK) {
					s *= 1.0e6;
				}

				RefGeomPrimaryX = PX();
				RefGeomPrimaryY = PY();
				RefGeomPrimaryZ = PZ();

				if (RefGeomPrimaryX == 0.0) Scale[XCOMP] = 0.0;
				else Scale[XCOMP] *= (s / RefGeomPrimaryX);

				if (RefGeomPrimaryY == 0.0) Scale[YCOMP] = 0.0;
				else Scale[YCOMP] *= (s / RefGeomPrimaryY);

				if (RefGeomPrimaryZ == 0.0) Scale[ZCOMP] = 0.0;
				else Scale[ZCOMP] *= (s / RefGeomPrimaryZ);
			}

		}

		void setup_discrete_frequencies() {
			double lf1 = log10(WvForm.BaseFrequency);
			double lf2 = log10(WvForm.SampleFrequency / 2);
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

		void setup_splines() {
			Comp.resize(NCOMP);
			Comp[XCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies);
			Comp[YCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies);
			Comp[ZCOMP].resize(NumberOfDiscreteFrequencies, NumberOfSplinedFrequencies);
			FrequencySpliner.initialise(DiscreteFrequenciesLog10, SplinedFrequencieslog10);
			LEM.init_frequencies(DiscreteFrequencies);
		}

		// Modelling
		void setearthproperties(const cEarth1D& E) {
			LEM.setproperties(E);
		}

		void setconductivitythickness(const size_t nlayers, const double* conductivity, const double* thickness){
			LEM.setconductivitythickness(nlayers, conductivity, thickness);
		}

		void setconductivitythickness(const std::vector<double>& conductivity, const std::vector<double>& thickness){
			LEM.setconductivitythickness(conductivity, thickness);
		}

		void setgeometry(const cTDEmGeometry& G) {
			Geometry = G;
			Vec3 tx_reference_orientation = Vec3::UnitZ();
			const Vec3 sep = Geometry.txrx_separation();
			const double& h = Geometry.tx_height;
			const double& x = sep.x();
			const double& y = sep.y();
			const double& z = h + sep.z();
			const Vec3 tx_orientation = Geometry.tx_orientation(Tx.Reference_Orientation);
			LEM.setgeometry(tx_orientation, h, x, y, z);
		};

		void setup_computations() {
			for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++) {
				LEM.init_frequency(fi);
			}
		}

		void setprimaryfields() {
			LEM.setprimaryfields();
			Vec3 v(LEM.Fields.t.p.x, LEM.Fields.t.p.y, LEM.Fields.t.p.z);

			if (LEM.calculation_type == cLEM::CalculationType::HDERIVATIVE) {
				//This is because when H changes Z also changes
				//and DZ = DH
				//but they should be all zero anyway
				v *= 2.0;
			}

			if (NormalisationType == NormalizationType::PPM_PEAKTOPEAK) {
				v *= 2.0;
			}

			if (OutputType == OutputType::DBDT) {
				//Must convert to dB/dt. This happens implicitly for the secondary via the waveform.
				v *= Tx.PeakdIdT;
			}

			// Rotate field to Rx frame
			const Mat3 RotMatrix = Geometry.inertial_to_rx_frame_rotation_matrix();
			v = RotMatrix * v;

			Comp[XCOMP].Primary = v.x() * Scale[XCOMP];
			Comp[YCOMP].Primary = v.y() * Scale[YCOMP];
			Comp[ZCOMP].Primary = v.z() * Scale[ZCOMP];
		};

		void setsecondaryfields() {
			//Computation for discrete frequencies 	
			const Mat3 RotMatrix = Geometry.inertial_to_rx_frame_rotation_matrix();
			for (size_t fi = 0; fi < NumberOfDiscreteFrequencies; fi++) {
				LEM.dointegrals(fi);
				LEM.setsecondaryfields(fi);
				const cdouble& x = LEM.Fields.t.s.x;
				const cdouble& y = LEM.Fields.t.s.y;
				const cdouble& z = LEM.Fields.t.s.z;
				Vec3 vr = Vec3(x.real(), y.real(), z.real());
				Vec3 vi = Vec3(x.imag(), y.imag(), z.imag());

				// Rotate field to Rx frame
				vr = RotMatrix * vr;
				vi = RotMatrix * vi;

				if (LEM.calculation_type == cLEM::CalculationType::HDERIVATIVE) {
					//This is because when H changes Z also changes
					//and DZ = DH
					vr *= 2.0;
					vi *= 2.0;
				}

				Comp[XCOMP].IR_discrete_real[fi] = vr.x();
				Comp[XCOMP].IR_discrete_imag[fi] = vi.x();
				Comp[YCOMP].IR_discrete_real[fi] = vr.y();
				Comp[YCOMP].IR_discrete_imag[fi] = vi.y();
				Comp[ZCOMP].IR_discrete_real[fi] = vr.z();
				Comp[ZCOMP].IR_discrete_imag[fi] = vi.z();
			};

			//Spline discreet frequencies		
			for (size_t i = 0; i < NCOMP; i++) {
				if (Scale[i] == 0.0) return;
				spline_component(i);
			}

			for (size_t i = 0; i < NCOMP; i++) {
				if (Scale[i] == 0.0) return;
				inverse_fft_window_scale_component(i);
			}

			if (SaveDiagnosticFiles) {
				write_discretefrequencies("diag_discretefrequencies.txt");
				write_splinedfrequencies("diag_splinedfrequencies.txt");
				WvForm.write_frequencydomainwaveform("diag_frequencydomainwaveform.txt");
			}

			if (SaveDiagnosticFiles) {
				WindScheme.write_windows("diag_windows.txt", XS(), YS(), ZS());
			}
		}

		void spline_component(const size_t& component) {
			ComponentWorkStore& C = Comp[component];
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

			ComponentWorkStore& C = Comp[component];
			// Reset to the stored transfer function
			WvForm.FFTWork = WvForm.Transfer;
			// Apply transfer function
			size_t n = 0;
			for (size_t k = 1; k < WvForm.NumFrequencies; k += 2) {
				WvForm.FFTWork[k] *= C.IR_splined[n];
				n++;
			}

			// Inverse FFT
			fftw_execute(inverse_fftplan);

			// Window
			C.Secondary = WindScheme.computewindow((double*)WvForm.FFTWork.data());
			if (SaveDiagnosticFiles) {
				write_timesseries("diag_xtimeseries.txt");
			}

			// Scale
			C.Secondary *= Scale[component];
		}

		void write_discretefrequencies(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumberOfDiscreteFrequencies; i++) {
				ofs << strprint("%15le\t%15le\t%15le\t%15le\t%15le\t%15le\t%15le\n", DiscreteFrequencies[i], 
					Comp[XCOMP].IR_discrete_real[i],
					Comp[XCOMP].IR_discrete_imag[i],
					Comp[YCOMP].IR_discrete_real[i],
					Comp[YCOMP].IR_discrete_imag[i],
					Comp[ZCOMP].IR_discrete_real[i],
					Comp[ZCOMP].IR_discrete_imag[i]);
			}
		}

		void write_splinedfrequencies(const fs::path& path) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t i = 0; i < NumberOfSplinedFrequencies; i++) {
				double f = pow10(SplinedFrequencieslog10[i]);
				ofs << strprint("%15le\t%15le\t%15le\t%15le\t%15le\t%15le\t%15le\n",
					f,
					Comp[XCOMP].IR_splined[i].real(),
					Comp[XCOMP].IR_splined[i].imag(),
					Comp[YCOMP].IR_splined[i].real(),
					Comp[YCOMP].IR_splined[i].imag(),
					Comp[ZCOMP].IR_splined[i].real(),
					Comp[ZCOMP].IR_splined[i].imag());
			}
		}

		void write_timesseries(const std::string& path) const {
			std::ofstream ofs = ofstream_ex(path);
			double* ts = (double*)(WvForm.FFTWork.data());
			for (size_t i = 0; i < WvForm.NumSamples; i++) {
				ofs << strprint("%20.10le\t%20.10le\n", WvForm.Time[i], ts[i]);
			}
		}

		void drx_pitch(double xb, double zb, double p, double& dxbdp, double& dzbdp) {
			//xi = (  xb*cosp  + zb*sinp);Inertial
			//zi = ( -xb*sinp  + zb*cosp);
			//xb = (  xi*cosp  - zi*sinp);As bird sees it
			//zb = (  xi*sinp  + zi*cosp);						

			if (NormalisationType == NormalizationType::PPM || NormalisationType == NormalizationType::PPM_PEAKTOPEAK) {
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

			if (NormalisationType == NormalizationType::PPM || NormalisationType == NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dxbdp /= RefGeomPrimaryX;
				dzbdp /= RefGeomPrimaryZ;
			}
		}

		void drx_pitch(std::vector<double> xb, std::vector<double> zb, double p, std::vector<double>& dxbdp, std::vector<double>& dzbdp) {
			//xi = (  xb*cosp  + zb*sinp);Inertial
			//zi = ( -xb*sinp  + zb*cosp);
			//xb = (  xi*cosp  - zi*sinp);As bird sees it
			//zb = (  xi*sinp  + zi*cosp);						

			if (NormalisationType == NormalizationType::PPM || NormalisationType == NormalizationType::PPM_PEAKTOPEAK) {
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

			if (NormalisationType == NormalizationType::PPM || NormalisationType == NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dxbdp /= RefGeomPrimaryX;
				dzbdp /= RefGeomPrimaryZ;
			}
		}

		void drx_roll(double yb, double zb, double r, double& dybdr, double& dzbdr)		{
			//yi = (  yb*cosr  - zb*sinr);Inertial
			//zi = (  yb*sinr  + zb*cosr);
			//yb = (  yi*cosr  + zi*sinr);As bird sees it
			//zb = ( -yi*sinr  + zi*cosr);						

			if (NormalisationType == NormalizationType::PPM || NormalisationType == NormalizationType::PPM_PEAKTOPEAK) {
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

			if (NormalisationType == NormalizationType::PPM || NormalisationType == NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dybdr /= RefGeomPrimaryY;
				dzbdr /= RefGeomPrimaryZ;
			}
		}

		void  drx_roll(std::vector<double> yb, std::vector<double> zb, double r, std::vector<double>& dybdr, std::vector<double>& dzbdr) {
			//yi = (  yb*cosr  - zb*sinr);Inertial
			//zi = (  yb*sinr  + zb*cosr);
			//yb = (  yi*cosr  + zi*sinr);As bird sees it
			//zb = ( -yi*sinr  + zi*cosr);						

			if (NormalisationType == NormalizationType::PPM || NormalisationType == NormalizationType::PPM_PEAKTOPEAK) {
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

			if (NormalisationType == NormalizationType::PPM || NormalisationType == NormalizationType::PPM_PEAKTOPEAK) {
				//Convert back to PPMS
				dybdr /= RefGeomPrimaryY;
				dzbdr /= RefGeomPrimaryZ;
			}
		}

		void set_response(cTDEmResponse& Response) const {
			Response.PX = PX();
			Response.PY = PY();
			Response.PZ = PZ();
			Response.SX = XS();
			Response.SY = YS();
			Response.SZ = ZS();
		};

		void forwardmodel(const cTDEmGeometry& G, const cEarth1D& E, cTDEmResponse& R)
		{
			setgeometry(G);
			setearthproperties(E);
			setup_computations();
			setprimaryfields();
			setsecondaryfields();
			set_response(R);
		}

		void forwardmodel(const std::vector<double>& conductivity, const std::vector<double>& thickness, const cTDEmGeometry& geometry)
		{
			setconductivitythickness(conductivity, thickness);
			setgeometry(geometry);
			LEM.calculation_type = cLEM::CalculationType::FORWARDMODEL;
			LEM.derivative_layer = undefinedvalue<size_t>();

			setup_computations();
			setprimaryfields();
			setsecondaryfields();
		}

		void forwardmodel(const size_t nlayers, const double* conductivity, const double* thickness, const double* g, double& PX, double& PY, double& PZ, double* SX, double* SY, double* SZ)
		{
			//Order const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw)
			cTDEmGeometry G(g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9]);
			setgeometry(G);
			setconductivitythickness(nlayers, conductivity, thickness);
			setup_computations();
			setprimaryfields();
			setsecondaryfields();
			getfields(PX, PY, PZ, SX, SY, SZ);
		}

		void getfields(double& PX, double& PY, double& PZ, double* SX, double* SY, double* SZ)
		{
			PX = PX;
			PY = PY;
			PZ = PZ;
			const size_t nw = nwindows();
			for (size_t i = 0; i < nw; i++) SX[i] = XS()[i];
			for (size_t i = 0; i < nw; i++) SY[i] = YS()[i];
			for (size_t i = 0; i < nw; i++) SZ[i] = ZS()[i];
		}

	};
};


