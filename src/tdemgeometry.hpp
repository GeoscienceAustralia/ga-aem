/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

//#include <vector>
//#include "general_constants.h"
//#include "geometry3d.h"

#include "blocklanguage.hpp"
#include "rollpitchyaw.hpp"
#include "aem_coredefs.hpp"

namespace AEM {

	inline static Mat3d YPR(const double& roll_degrees, const double& pitch_degrees, const double& yaw_degrees) {
		const Mat3d Rot = yawpitchroll_matrix(roll_degrees * D2R<double>, pitch_degrees * D2R<double>, yaw_degrees * D2R<double>);
		return Rot;
	};

	inline static Mat3d invYPR(const double& roll_degrees, const double& pitch_degrees, const double& yaw_degrees) {
		const Mat3d Rot = yawpitchroll_matrix(roll_degrees * D2R<double>, pitch_degrees * D2R<double>, yaw_degrees * D2R<double>);
		Mat3d RotT = Rot.transpose();
		return RotT;
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
		};

		void write(std::string path) const {
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
};