#pragma once

#include <cassert>
#include <iostream>
#include "numerical_utils.hpp"
#include "vector_utils.hpp"

namespace AEM {

	template <typename T>
	class TDEmScalarResponse {

	private:
		std::vector<T> v;

	public:

		TDEmScalarResponse() {};

		TDEmScalarResponse(const size_t& nwindows) {
			resize(nwindows);
		}

		inline const size_t size() const { return v.size(); }

		T& operator[](const size_t& i) {
			return v[i];
		}

		T operator[](const size_t& i) const {
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

		std::vector<T> storage() const {
			return v;
		}

	private:

		void  resize(const size_t nwindows) {
			v.resize(nwindows);
		}
	};

	template <typename T>
	class TDEmVectorResponse {

	private:
		size_t nwindows = 0;
		std::array<std::vector<T>, NCOMP> v;
				
	public:

		TDEmVectorResponse(const size_t _nwindows = 0) {
			set_nWindows(_nwindows);
		};

		inline const size_t nWindows() const { return nwindows; }

		void set_nWindows(const size_t _nwindows) {
			nwindows = _nwindows;
			v[0].resize(nwindows);
			v[1].resize(nwindows);
			v[2].resize(nwindows);
		};

		Eigen::Vector<T,3> get_vec3(const size_t window) const {
			return Eigen::Vector<T,3>(v[0][window], v[1][window], v[2][window]);
		};

		void set_vec3(const size_t window, const Eigen::Vector<T,3>& vec) {
			v[0][window] = vec[0];
			v[1][window] = vec[1];
			v[2][window] = vec[2];
		};

		std::vector<T>& operator[](const size_t& component) {
			return v[component];
		}

		const std::vector<T>& operator[](const size_t& component) const {
			return v[component];
		}

		TDEmVectorResponse& operator+=(const TDEmVectorResponse& rhs) {
			v[0] += rhs.v[0];
			v[1] += rhs.v[1];
			v[2] += rhs.v[2];
			return *this;
		}

		TDEmVectorResponse& operator-=(const TDEmVectorResponse& rhs) {
			v[0] -= rhs.v[0];
			v[1] -= rhs.v[1];
			v[2] -= rhs.v[2];
			return *this;
		}

		TDEmVectorResponse& operator*=(const double& s) {
			v[0] *= s;
			v[1] *= s;
			v[2] *= s;
			return *this;
		}

		TDEmVectorResponse& operator/=(const double& s) {
			v[0] /= s;
			v[1] /= s;
			v[2] /= s;
			return *this;
		}

		T& operator()(const size_t& component, const size_t& window) {
			assert(component < NCOMP);
			assert(window < nWindows());
			return v[component][window];
		}

		void scale_components(const Vec3d& scalefactors) {
			const size_t nw = v.size();
			v[0] *= scalefactors[0];
			v[1] *= scalefactors[1];
			v[2] *= scalefactors[2];
		};

		void divide_components(const Vec3d& dividefactors) {
			v[0] /= dividefactors[0];
			v[1] /= dividefactors[1];
			v[2] /= dividefactors[2];
		};

		void plus_components(const Vec3d& values) {
			v[0] += values[0];
			v[1] += values[1];
			v[2] += values[2];
		};

		void minus_components(const Vec3d& values) {
			v[0] -= values[0];
			v[1] -= values[1];
			v[2] -= values[2];
		};

		void set_values(const T& value) {
			set(v[0], value);
			set(v[1], value);
			set(v[2], value);
		};

		void set_values(const Vec3d& value) {
			set(v[0], value[0]);
			set(v[1], value[1]);
			set(v[2], value[2]);
		};

		void set_values(const Vec3cd& value) {
			set(v[0], value[0]);
			set(v[1], value[1]);
			set(v[2], value[2]);
		};

		TDEmScalarResponse<T> component(const size_t& ci) const {
			TDEmScalarResponse<T> c(nwindows);
			for (size_t wi = 0; wi < nwindows; wi++) {
				c[wi] = v[ci][wi];
			}
			return c;
		};

		TDEmScalarResponse<T> xzamp() const {
			TDEmScalarResponse<T> r(nwindows);
			for (size_t i = 0; i < nwindows; i++) {
				r[i] = AEM::ewise_hypot(v[XCOMP][i], v[ZCOMP][i]);
			}
			return r;
		};

		TDEmScalarResponse<T> xyzamp() const {
			TDEmScalarResponse<T> r(nwindows);
			for (size_t i = 0; i < nwindows; i++) {
				r[i] = AEM::ewise_hypot(v[XCOMP][i], v[YCOMP][i], v[ZCOMP][i]);
			}
			return r;
		};

		friend std::ostream& operator<<(std::ostream& os, const TDEmVectorResponse& R) {			
			for (size_t i = 0; i < R.nwindows; i++) { 
			//for (size_t i = 0; i < 2; i++) {
				os  << exd(16, 6) << R[XCOMP][i]
					<< exd(16, 6) << R[YCOMP][i]
					<< exd(16, 6) << R[ZCOMP][i]
					<< std::endl;
			}
			return os;
		}

		const static std::string ss(const double& v) {
			std::ostringstream os;
			os << exd(16, 6) << v;
			return os.str();
		}

		const static std::string ss(const std::complex<double>& v) {
			std::ostringstream os;
			os << exd(16, 6) << v.real() << exd(16, 6) << v.imag();
			return os.str();
		}

		void simple_output(std::ostream& os) const {
			for (size_t i = 0; i < nwindows; i++) {
				os << ss(v[XCOMP][i]);
				os << ss(v[YCOMP][i]);
				os << ss(v[ZCOMP][i]);
				os << std::endl;
			}
		}

	private:

	};

	template <typename T>
	class TDEmResponse {

	public:
		TDEmVectorResponse<T> P;
		TDEmVectorResponse<T> S;

		TDEmResponse() {};

		TDEmResponse(const size_t& _nwindows) {
			set_nWindows(_nwindows);
		}

		void set_nWindows(const size_t& _nwindows) {
			P.set_nWindows(_nwindows);
			S.set_nWindows(_nwindows);
		}

		const size_t& nWindows() const {
			return S.nWindows();
		}

		const T primary(const size_t& component, const size_t& window) const {
			assert(component < NCOMP);
			assert(window < nWindows());
			return P[component][window];
		};

		const std::vector<T> primary(const size_t& component) const {
			assert(component < NCOMP);
			return P[component];
		};


		const T secondary(const size_t& component, const size_t& window) const {
			assert(component < NCOMP);
			assert(window < nWindows());
			return S[component][window];
		};

		const std::vector<T> secondary(const size_t& component) const {
			assert(component < NCOMP);
			return S[component];
		};

		const std::vector<T> total(const size_t& component) const {
			assert(component < NCOMP);
			return primary(component) + secondary(component);
		};

		std::vector<T> psi(const size_t& component, const double& scaled_ga_component) const {
			assert(component < NCOMP);
			std::vector<T> t = primary(component) + secondary(component);
			t /= scaled_ga_component;
			return t;
		};

		std::vector<T> psi_xzamp(const Vec3d& scaled_ga) const {
			std::vector<T> psix = psi(XCOMP, scaled_ga[XCOMP]);
			std::vector<T> psiz = psi(ZCOMP, scaled_ga[ZCOMP]);
			return AEM::ewise_hypot(psix,psiz);
		};

		TDEmVectorResponse<T> totalfield() const {
			TDEmVectorResponse T = S;
			T += P;
			return T;
		};

		TDEmResponse& operator*=(const double& s) {
			P *= s;
			S *= s;
			return *this;
		};

		TDEmResponse& operator/=(const double& s) {
			P /= s;
			S /= s;
			return *this;
		};

		TDEmResponse& operator+=(const TDEmResponse& rhs) {
			P += rhs.P;
			S += rhs.S;
			return *this;
		};

		TDEmResponse& operator-=(const TDEmResponse& rhs) {
			P -= rhs.P;
			S -= rhs.S;
			return *this;
		};

		friend TDEmResponse operator+(const TDEmResponse& a, const TDEmResponse& b) {
			TDEmResponse r = a;
			r += b;
			return r;
		};

		friend TDEmResponse operator-(const TDEmResponse& a, const TDEmResponse& b) {
			TDEmResponse r = a;
			r -= b;
			return r;
		};

		friend TDEmResponse operator*(const TDEmResponse& a, const double& s) {
			TDEmResponse r = a;
			r *= s;
			return r;
		};

		friend TDEmResponse operator*(const double& s, const TDEmResponse& a) {
			TDEmResponse r = a;
			r *= s;
			return r;
		};

		friend TDEmResponse operator/(const TDEmResponse& a, const double& s) {
			TDEmResponse r = a;
			r /= s;
			return r;
		};

		friend TDEmResponse elementwise_div(const TDEmResponse& a, const TDEmResponse& b) {
			TDEmResponse r = a;
			r.P[0] /= b.P[0];
			r.P[1] /= b.P[1];
			r.P[2] /= b.P[2];
			r.S[0] /= b.S[0];
			r.S[1] /= b.S[1];
			r.S[2] /= b.S[2];
			return r;
		};

		friend TDEmResponse percent_difference(const TDEmResponse& a, const TDEmResponse& b, const double tol = std::numeric_limits<double>::epsilon()) {
			TDEmResponse r(a.nWindows());
			for (size_t ci = 0; ci < NCOMP; ci++) {
				for (size_t wi = 0; wi < r.nWindows(); wi++) {
					r.P[ci][wi] = pct_diff_ex(a.P[ci][wi], b.P[ci][wi], tol);
					r.S[ci][wi] = pct_diff_ex(a.S[ci][wi], b.S[ci][wi], tol);
				}
			}
			return r;
		};
		
		static void xxxdisplay_max_abs_percent_difference(const TDEmResponse<double>& PCD) {
			fxd fmt = fxd(12, 6);
			std::cout << "P (%): ";
			std::cout << fmt << std::max(std::abs(min(PCD.P[XCOMP])), std::abs(max(PCD.P[XCOMP]))) << " ";
			std::cout << fmt << std::max(std::abs(min(PCD.P[YCOMP])), std::abs(max(PCD.P[YCOMP]))) << " ";
			std::cout << fmt << std::max(std::abs(min(PCD.P[ZCOMP])), std::abs(max(PCD.P[ZCOMP]))) << std::endl;
			std::cout << "S (%): ";
			std::cout << fmt << std::max(std::abs(min(PCD.S[XCOMP])), std::abs(max(PCD.S[XCOMP]))) << " ";
			std::cout << fmt << std::max(std::abs(min(PCD.S[YCOMP])), std::abs(max(PCD.S[YCOMP]))) << " ";
			std::cout << fmt << std::max(std::abs(min(PCD.S[ZCOMP])), std::abs(max(PCD.S[ZCOMP]))) << std::endl;
			std::cout << std::endl;
		};

		static void display_max_abs_percent_difference(const TDEmResponse<double>& PCD) {
			fxd fmt = fxd(10, 6);
			std::cout << "P (%): ";
			std::cout << fmt << maxabs(PCD.P[XCOMP]) << " ";
			std::cout << fmt << maxabs(PCD.P[YCOMP]) << " ";
			std::cout << fmt << maxabs(PCD.P[ZCOMP]) << std::endl;
			std::cout << "S (%): ";
			std::cout << fmt << maxabs(PCD.S[XCOMP]) << " ";
			std::cout << fmt << maxabs(PCD.S[YCOMP]) << " ";
			std::cout << fmt << maxabs(PCD.S[ZCOMP]) << std::endl;
			std::cout << std::endl;
		};

		static void display_max_abs_percent_difference(const TDEmResponse<cdouble>& PCD) {
			fxd fmt = fxd(12, 8);
			std::cout << "P (%): "
				<< "[" << fmt << maxabs(real(PCD.P[XCOMP])) << "," << fmt << maxabs(imaginary(PCD.P[XCOMP])) << "]  "
				<< "[" << fmt << maxabs(real(PCD.P[YCOMP])) << "," << fmt << maxabs(imaginary(PCD.P[YCOMP])) << "]  "
				<< "[" << fmt << maxabs(real(PCD.P[ZCOMP])) << "," << fmt << maxabs(imaginary(PCD.P[ZCOMP])) << "]" << std::endl;
			std::cout << "S (%): " 
				<< "[" << fmt << maxabs(real(PCD.S[XCOMP])) << "," << fmt << maxabs(imaginary(PCD.S[XCOMP])) << "]  "
				<< "[" << fmt << maxabs(real(PCD.S[YCOMP])) << "," << fmt << maxabs(imaginary(PCD.S[YCOMP])) << "]  "
				<< "[" << fmt << maxabs(real(PCD.S[ZCOMP])) << "," << fmt << maxabs(imaginary(PCD.S[ZCOMP])) << "]" << std::endl;
			std::cout << std::endl;
		};
		

		friend std::ostream& operator<<(std::ostream& os, const TDEmResponse& R) {
			//os << "--Primary--" << std::endl;
			//os << R.P;
			os << "--Secondary--" << std::endl;
			os << R.S;
			return os;
		};
	};
};
