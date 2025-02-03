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

	private:

		void  resize(const size_t nwindows) {
			v.resize(nwindows);
		}
	};

	template <typename T>
	class TDEmVectorResponse {

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

		TDEmScalarResponse<T> xzamp() {
			TDEmScalarResponse<T> r(nwindows);
			for (size_t i = 0; i < nwindows; i++) {
				r[i] = AEM::hypot(v[XCOMP][i], v[ZCOMP][i]);
			}
			return r;
		};

		friend std::ostream& operator<<(std::ostream& os, const TDEmVectorResponse& R) {			
			for (size_t i = 0; i < R.nwindows; i++) { 
				os  << exd(16, 6) << R[XCOMP][i]
					<< exd(16, 6) << R[YCOMP][i]
					<< exd(16, 6) << R[ZCOMP][i]
					<< std::endl;
			}
			return os;
		}

		void simple_output(std::ostream& os) {
			for (size_t i = 0; i < nwindows; i++) {
				os <<  exd(16, 6) << v[XCOMP][i].real()
					<< exd(16, 6) << v[XCOMP][i].imag()
					<< exd(16, 6) << v[YCOMP][i].real()
					<< exd(16, 6) << v[YCOMP][i].imag()
					<< exd(16, 6) << v[ZCOMP][i].real()
					<< exd(16, 6) << v[ZCOMP][i].imag()
					<< std::endl;
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

		friend TDEmResponse percent_difference(const TDEmResponse& a, const TDEmResponse& b) {
			TDEmResponse r = 100.0 * elementwise_div(b - a, a);
			constexpr double eps = std::numeric_limits<double>::epsilon();
			for (size_t ci = 0; ci < NCOMP; ci++) {
				for (size_t wi = 0; wi < r.nWindows(); wi++) {
					// Amend for closeness within numerical precision
					if (nearly_equal_ulps(a.P[ci][wi], b.P[ci][wi])) r.P[ci][wi] = 0.0;
					else if (a.P[ci][wi] == 0.0 && b.P[ci][wi] == 0.0) r.P[ci][wi] = 0.0;
					else if (std::abs(a.P[ci][wi]) <= eps && std::abs(b.P[ci][wi] <= eps)) r.P[ci][wi] = 0.0;

					if (nearly_equal_ulps(a.S[ci][wi], b.S[ci][wi])) r.S[ci][wi] = 0.0;
					else if (a.S[ci][wi] == 0.0 && b.S[ci][wi] == 0.0) r.S[ci][wi] = 0.0;
					else if (std::abs(a.S[ci][wi]) <= eps && std::abs(b.S[ci][wi] <= eps)) r.S[ci][wi] = 0.0;

				}
			}
			return r;
		};

		static void display_max_abs_percent_difference(const TDEmResponse& PCD) {
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

		friend std::ostream& operator<<(std::ostream& os, const TDEmResponse& R) {
			os << "--Primary--" << std::endl;
			os << R.P;
			os << "--Secondary--" << std::endl;
			os << R.S;
			return os;
		};
	};
};
