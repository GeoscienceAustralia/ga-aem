/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include "general_utils.hpp"
#include "blocklanguage.hpp"
#include "vector_utils.hpp"
#include "aem_types.hpp"
#include "calculation_type.hpp"
#include "layeredearthmodeller.hpp"
#include "lem.hpp"

namespace AEM {

	using namespace CppUtils;

	//using namespace VectorUtils;
	// This allows the VectorUtils operator overloads to be found from this namespace 
	// ... see https://www.reddit.com/r/cpp_questions/comments/17mdzzx/overload_resolution_fails_to_find_operator_from_a/?rdt=36798
	using CppUtils::operator+;
	using CppUtils::operator-;
	using CppUtils::operator*;
	using CppUtils::operator/;
	using CppUtils::operator+=;
	using CppUtils::operator-=;
	using CppUtils::operator<<;

	using namespace LEM2;

	enum class SystemType { TimeDomain, SpectralTimeDomain };

	inline constexpr size_t XCOMP = 0;
	inline constexpr size_t YCOMP = 1;
	inline constexpr size_t ZCOMP = 2;
	inline constexpr size_t NCOMP = 3;

	// Get AEM SystemType from string
	static AEM::SystemType aem_system_type(const std::string& typestr) {
		if (strcasecmp(typestr, "Spectral Time Domain") == 0) return SystemType::SpectralTimeDomain;
		else if (strcasecmp(typestr, "Time Domain") == 0) return SystemType::TimeDomain;
		else {
			glog.errormsg(_SRC_, "The AEM System 'Type' must be 'Spectral Time Domain' or 'Time Domain'.\n");
		}
		//never gets here but suppresses warning
		return aem_system_type(typestr);
	};

	// Get AEM SystemType from stmfile
	static AEM::SystemType aem_system_type(const fs::path& stmfile) {
		cBlock b(stmfile);
		std::string typestr;
		if (b.getvalue("Type", typestr)) return aem_system_type(typestr);
		else glog.errormsg(_SRC_, "The AEM System 'Type' is not specified.\n");
		
		//never gets here but suppresses warning
		return aem_system_type(typestr);
	};

	// Get AEM SystemType from template parameter
	template<typename AEMSystemClass>
	inline AEM::SystemType aem_system_type() {
		return AEMSystemClass::get_type();
	};
	
	template <typename T>
	std::vector<T> real(const std::vector<std::complex<T>>& cv) {
		const size_t n = cv.size();
		std::vector<T> v(n);
		for (size_t i = 0; i < n; i++) {
			v[i] = cv[i].real();
		}
		return v;
	};

	template <typename T>
	std::vector<T> imaginary(const std::vector<std::complex<T>>& cv) {
		const size_t n = cv.size();
		std::vector<T> v(n);
		for (size_t i = 0; i < n; i++) {
			v[i] = cv[i].imag();
		}
		return v;
	};


	// Utilities for handling double and complex<double> the same
	inline double ewise_hypot(const double& x, const double& y) {
		return std::hypot(x, y);
	};

	inline double ewise_hypot(const double& x, const double& y, const double& z) {
		return std::hypot(x, y, z);
	};

	inline cdouble ewise_hypot(const cdouble& x, const cdouble& y) {
		return cdouble(std::hypot(x.real(), y.real()), std::hypot(x.imag(), y.imag()));
	};

	inline cdouble ewise_hypot(const cdouble& x, const cdouble& y, const cdouble& z) {
		return cdouble(std::hypot(x.real(), y.real(), z.real()), std::hypot(x.imag(), y.imag(), z.imag()));
	};

	inline std::vector<double>ewise_hypot(const std::vector<double>& x, const std::vector<double>& y) {
		const size_t n = x.size();
		std::vector<double> h(n);
		for (size_t i = 0; i < n; i++) {
			h[i]  = ewise_hypot(x[i], y[i]);
		}
		return h;
	};

	inline std::vector<cdouble> ewise_hypot(const std::vector<cdouble>& x, const std::vector<cdouble>& y) {
		const size_t n = x.size();
		std::vector<cdouble> h(n);
		for (size_t i = 0; i < n; i++) {
			h[i] = ewise_hypot(x[i], y[i]);
		}
		return h;
	};

	inline double ewise_mul(const double& x, const double& y) {
		return x * y;
	};

	inline cdouble ewise_mul(const cdouble& x, const cdouble& y) {
		return cdouble(x.real() * y.real() , x.imag() * y.imag());
	};

	inline double ewise_div(const double& x, const double& y) {
		return x / y;
	};

	inline cdouble ewise_div(const cdouble& x, const cdouble& y) {
		return cdouble(x.real() / y.real(), x.imag() / y.imag());
	};

	template <typename T>
	void complex_merge(const std::vector<T>& realv, const std::vector<T>& imagv, std::vector<std::complex<T>>& v) {
		assert(realv.size() == imagv.size());
		const size_t n = realv.size();
		v.resize(n);
		for (size_t k = 0; k < n; k++) {
			v[k] = std::complex<T>(realv[k], imagv[k]);
		}
	}

	template <typename T>
	std::vector<std::complex<T>> complex_merge(const std::vector<T>& realv, const std::vector<T>& imagv) {
		std::vector<std::complex<T>> v;
		complexmerge(realv, imagv, v);
		return v;
	};

	template <typename T>
	size_t value_size() { return 1; };
	template<> size_t value_size<double>() { return 1; };
	template<> size_t value_size<cdouble>() { return 2; };

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

	class Transmitter {

	public:
		double BaseFrequency = 0.0;
		size_t nTurns = 1;
		double LoopArea = 1.0;
		double PeakCurrent = 1.0;
		double PeakdIdT = 1.0;
		Vec3d Reference_Orientation = Vec3d::UnitZ();
		Vec3d Orientation = Vec3d::UnitZ();

		Transmitter() {};

		Transmitter(const cBlock& b) {
			if (b.getvalue("BaseFrequency", BaseFrequency) == false) {
				glog.errormsg(_SRC_, "A Transmitter BaseFrequency must be specified.\n");
			};
			if (b.getvalue("PeakCurrent", PeakCurrent) == false) PeakCurrent = 1.0;
			if (b.getvalue("LoopArea", LoopArea) == false) LoopArea = 1.0;
			if (b.getvalue("NumberOfTurns", nTurns) == false) nTurns = 1;
		};
	};

	struct WindowSpecification {
		size_t SampleLow = 0;
		size_t SampleHigh = 0;
		size_t NumberOfSamples = 0;

		double Low = 0.0;
		double High = 0.0;
		double Width = 0.0;

		std::vector<size_t> Sample;
		std::vector<double> Weight;

		double centre() const {
			return (Low + High) / 2.0;
		}
	};

	class WindowingScheme {

	private:


	public:

		enum class WeightingMethod { BoxCar, AreaUnderCurve, LinearTaper };

		size_t nwindows = 0;
		double Shift = 0.0;
		WeightingMethod Method = WeightingMethod::BoxCar;
		std::vector<WindowSpecification> Windows;

		WindowingScheme() {};

		/*
		WindowingScheme(const cBlock& receiverblock, const Waveform& WFM) {
			const cBlock& b = receiverblock;

			if (!b.getvalue("NumberOfWindows", nwindows)) {
				glog.errormsg(_SRC_, "NumberOfWindows is not specified");
			}

			if (!b.getvalue("Shift", Shift)) {
				Shift = 0.0;
			}

			Windows.resize(nwindows);

			//Read window times
			std::vector<std::vector<double>> wt;
			if (!b.getvalue("WindowTimes", wt)) {
				glog.errormsg(_SRC_, "The WindowTimes have not been specified.");
			}
			size_t nw = wt.size();
			if (nw != nwindows) {
				glog.errormsg(_SRC_, "The number of WindowTimes does not match the NumberOfWindows\n");
			}

			for (size_t i = 0; i < nwindows; i++) {
				if (wt[i].size() != 2) {
					glog.errormsg(_SRC_, "The number of WindowTimes must have exactly 2 columns (error in window %lu)\n", i + 1);
				}
				Windows[i].Low = wt[i][0] + Shift;
				Windows[i].High = wt[i][1] + Shift;
			}

			std::string wmethod = b.getstringvalue("WindowWeightingScheme");
			if (strcasecmp(wmethod, "AreaUnderCurve") == 0) {
				initialise_area(WFM.Time);
			}
			else if (strcasecmp(wmethod, "Boxcar") == 0) {
				initialise_boxcar(WFM.Time);
			}
			else if (strcasecmp(wmethod, "LinearTaper") == 0) {
				initialise_lineartaper(WFM.Time);
			}
			else glog.errormsg(_SRC_, "WindowWeightingScheme %s unknown (must be \"AreaUnderCurve\" or  \"Boxcar\" or \"LinearTaper\")\n", wmethod.c_str());
		}
		*/
		WindowingScheme(const cBlock& receiverblock, const std::vector<double> series) {
			const cBlock& b = receiverblock;

			if (!b.getvalue("NumberOfWindows", nwindows)) {
				glog.errormsg(_SRC_, "NumberOfWindows is not specified");
			}

			if (!b.getvalue("Shift", Shift)) {
				Shift = 0.0;
			}

			Windows.resize(nwindows);

			//Read window times
			std::vector<std::vector<double>> w;
			if(b.getvalue("WindowTimes", w)) {
				glog.warningmsg("'WindowTimes' is deprecated. Please use 'Windows' instead.\n");
			}
			else if (!b.getvalue("Windows", w)) {
				glog.errormsg(_SRC_, "The Windows have not been specified.");
			}
			size_t nw = w.size();
			if (nw != nwindows) {
				glog.errormsg(_SRC_, "The number of Windows does not match the NumberOfWindows\n");
			}

			for (size_t i = 0; i < nwindows; i++) {
				if (w[i].size() != 2) {
					glog.errormsg(_SRC_, "The number of WindowFrequencies must have exactly 2 columns (error in window %lu)\n", i + 1);
				}
				Windows[i].Low = w[i][0] + Shift;
				Windows[i].High = w[i][1] + Shift;
			}

			std::string wmethod = b.getstringvalue("WindowWeightingScheme");
			if (strcasecmp(wmethod, "AreaUnderCurve") == 0) {
				initialise_area(series);
			}
			else if (strcasecmp(wmethod, "Boxcar") == 0) {
				initialise_boxcar(series);
			}
			else if (strcasecmp(wmethod, "LinearTaper") == 0) {
				initialise_lineartaper(series);
			}
			else glog.errormsg(_SRC_, "WindowWeightingScheme %s unknown (must be \"AreaUnderCurve\" or  \"Boxcar\" or \"LinearTaper\")\n", wmethod.c_str());
		}

		const size_t& nWindows() const { return nwindows; }

		void initialise_area(const std::vector<double>& series) {
			double dwt = series[1] - series[0];
			double eps = 1.0e-7;

			for (size_t w = 0; w < nwindows; w++) {
				Windows[w].Width = Windows[w].High - Windows[w].Low;
				for (size_t s = 0; s < series.size(); s++) {
					const double& t = series[s];
					if (t + eps >= Windows[w].Low) {
						Windows[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = series.size(); s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					const double& t = series[s];
					if (t - eps <= Windows[w].High) {
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

				const double& tlow = Windows[w].Low;
				const double& thigh = Windows[w].High;
				double wsum = 0;
				for (size_t k = 0; k < Windows[w].NumberOfSamples; k++) {
					size_t s = Windows[w].Sample[k];
					const double& t = series[s];
					const double tp = series[s] - dwt;
					const double tn = series[s] + dwt;
					const double tleft = std::max(tp, tlow);
					const double tright = std::min(tn, thigh);
					Windows[w].Weight[k] = 0.5 * (t - tleft) + 0.5 * (tright - t);
					Windows[w].Weight[k] /= Windows[w].Width;
					wsum += Windows[w].Weight[k];
				}
				for (size_t k = 0; k < Windows[w].NumberOfSamples; k++) {
					Windows[w].Weight[k] /= wsum;
				}
			}
		}

		void initialise_boxcar(const std::vector<double>& series) {
			double eps = 1.0e-7;
			for (size_t w = 0; w < nwindows; w++) {
				Windows[w].Width = Windows[w].High - Windows[w].Low;
				for (size_t s = 0; s < series.size(); s++) {
					double t = series[s];
					if (t + eps >= Windows[w].Low) {
						Windows[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = series.size(); s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					double t = series[s];
					if (t - eps <= Windows[w].High) {
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

		void initialise_lineartaper(const std::vector<double>& series) {
			double eps = 1.0e-7;
			for (size_t w = 0; w < nwindows; w++) {
				Windows[w].Width = Windows[w].High - Windows[w].Low;
				for (size_t s = 0; s < series.size(); s++) {
					double t = series[s];
					if (t + eps >= Windows[w].Low) {
						Windows[w].SampleLow = s;
						break;
					}
				}

				for (size_t s = series.size(); s-- > 0;) {
					//Note the unusual syntax for decrement of unsigned variable
					double t = series[s];
					if (t - eps <= Windows[w].High) {
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

		template<typename T>
		void computewindow(const T* values, std::vector<T>& windowed_values) {
			std::fill(windowed_values.begin(), windowed_values.end(), T(0.0)); // Reset to zero
			for (size_t w = 0; w < nwindows; w++) {
				for (size_t k = 0; k < Windows[w].Sample.size(); k++) {
					windowed_values[w] += values[Windows[w].Sample[k]] * Windows[w].Weight[k];
				}
			}
		}

		template<typename T>
		void computewindow(const T* values, Eigen::Matrix<T,-1,1>& windowed_values) {
			std::fill(windowed_values.begin(), windowed_values.end(), T(0.0)); // Reset to zero
			for (size_t w = 0; w < nwindows; w++) {
				for (size_t k = 0; k < Windows[w].Sample.size(); k++) {
					windowed_values[w] += values[Windows[w].Sample[k]] * Windows[w].Weight[k];
				}
			}
		}

		void printwindows(const double& PX, const double& PY, const double& PZ, const std::vector<double>& SX, const std::vector<double>& SY, const std::vector<double>& SZ) const {
			printf("Primary   %15.8lf%15.8lf%15.8lf\n\n", PX, PY, PZ);
			printf("Window#             X               Y               Z\n");
			for (size_t w = 0; w < nwindows; w++) {
				printf("%2zu        %15.8lf%15.8lf%15.8lf\n", w + 1, SX[w], SY[w], SZ[w]);
			}
		};

		template<typename T, template<typename> typename V>
		void write_windows(const fs::path& path, 
			const V<T>& SX, 
			const V<T>& SY, 
			const V<T>& SZ) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t w = 0; w < nwindows; w++) {
				ofs << strprint("%2zu\t%20e\t%20e\t%15e%15e%15e\n", w + 1, Windows[w].Low, Windows[w].High, SX[w], SY[w], SZ[w]);
			}
		};

		template<typename T, template<typename> typename V>
		void write_windows(const fs::path& path, const V<std::complex<T>>& SX, const V<std::complex<T>>& SY, const V<std::complex<T>>& SZ) const {
			std::ofstream ofs = ofstream_ex(path);
			for (size_t w = 0; w < nwindows; w++) {
				ofs << strprint("%2zu\t%20e\t%20e\t%15e%15e%15e\n", w + 1, Windows[w].Low, Windows[w].High, SX[w], SY[w], SZ[w]);
			}
		};
	};

	class Receiver {

	public:
		double SamplingFrequency = 0.0;
		double Area = 1.0;
		size_t nTurns = 1;
		//Vec3d Reference_Orientation = Vec3d::UnitZ();

		Receiver() {};

		Receiver(const cBlock& b) {
			if (b.getvalue("SamplingFrequency", SamplingFrequency) == false) {
				glog.errormsg(_SRC_,"A Receiver SamplingFrequency must be specified.\n");
			};
			if(b.getvalue("CoilArea", Area) == false) Area = 1.0;
			if (b.getvalue("NumberOfTurns", nTurns) == false) nTurns = 1;
		};
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

};



