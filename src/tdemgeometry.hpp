/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include <array>
//#include <vector>
//#include "general_constants.h"
//#include "geometry3d.h"

#include "blocklanguage.hpp"
#include "rollpitchyaw.hpp"
#include "aem_coredefs.hpp"
#include "calculation_type.hpp"

namespace AEM {

	static constexpr const char TX_HEIGHT[] = "TX_HEIGHT";
	static constexpr const char TX_ROLL[] = "TX_ROLL";
	static constexpr const char TX_PITCH[] = "TX_PITCH";
	static constexpr const char TX_YAW[] = "TX_YAW";
	static constexpr const char TXRX_DX[] = "TXRX_DX";
	static constexpr const char TXRX_DY[] = "TXRX_DY";
	static constexpr const char TXRX_DZ[] = "TXRX_DZ";
	static constexpr const char RX_ROLL[] = "RX_ROLL";
	static constexpr const char RX_PITCH[] = "RX_PITCH";
	static constexpr const char RX_YAW[] = "RX_YAW";

	class TDEmGeometry {

	public:
		enum class ElementType { tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw, unknown};
		static constexpr size_t NELEM = 10;

	private:
		
		inline static double undefined = undefinedvalue<double>();
		inline static const std::array<std::string, NELEM> defined_names{TX_HEIGHT, TX_ROLL, TX_PITCH, TX_YAW, TXRX_DX, TXRX_DY, TXRX_DZ, RX_ROLL, RX_PITCH, RX_YAW };
		inline static const std::array<std::string, NELEM> defined_units{"m", "degrees", "degrees", "degrees", "m", "m", "m", "degrees", "degrees", "degrees" };
		inline static const std::array<ElementType, NELEM> defined_elementtypes{ ElementType::tx_height, ElementType::tx_roll, ElementType::tx_pitch, ElementType::tx_yaw, ElementType::txrx_dx, ElementType::txrx_dy, ElementType::txrx_dz, ElementType::rx_roll, ElementType::rx_pitch, ElementType::rx_yaw };
		inline static const std::array<CalculationType::Mode, NELEM> defined_derivative_modes{ CalculationType::Mode::DTX_HEIGHT,
			CalculationType::Mode::DTX_ROLL,
			CalculationType::Mode::DTX_PITCH,
			CalculationType::Mode::DTX_YAW,
			CalculationType::Mode::DX,
			CalculationType::Mode::DY,
			CalculationType::Mode::DZ,
			CalculationType::Mode::DRX_ROLL,
			CalculationType::Mode::DRX_PITCH,
			CalculationType::Mode::DRX_YAW };
		inline static const std::array<std::string, NELEM> defined_descriptions{"Tx height above ground level", "Tx roll - left side up + ve", "Tx pitch - nose down + ve", "Tx yaw - turn left + ve", "Tx - Rx horizonatl inline separation", "Tx - Rx horizonatl transverse separation", "Tx - Rx vertical separation", "Rx roll - left side up + ve", "Rx pitch - nose down + ve", "Rx yaw - turn left + ve" };

		std::array<double,NELEM> _elements_;

	public:

		inline static Mat3d YPR(const double& roll_degrees, const double& pitch_degrees, const double& yaw_degrees) {
			const Mat3d M = yawpitchroll_matrix(roll_degrees * D2R<double>, pitch_degrees * D2R<double>, yaw_degrees * D2R<double>);
			return M;
		};

		inline static Mat3d invYPR(const double& roll_degrees, const double& pitch_degrees, const double& yaw_degrees) {
			const Mat3d M = yawpitchroll_matrix(roll_degrees * D2R<double>, pitch_degrees * D2R<double>, yaw_degrees * D2R<double>);
			Mat3d MT = M.transpose();
			return MT;
		};


		double& tx_height() { return _elements_[0]; };
		double& tx_roll() { return _elements_[1]; };
		double& tx_pitch() { return _elements_[2]; };
		double& tx_yaw() { return _elements_[3]; };
		double& txrx_dx() { return _elements_[4]; };
		double& txrx_dy() { return _elements_[5]; };
		double& txrx_dz() { return _elements_[6]; };
		double& rx_roll() { return _elements_[7]; };
		double& rx_pitch() { return _elements_[8]; };
		double& rx_yaw() { return _elements_[9]; };

		const double& tx_height() const { return _elements_[0]; };
		const double& tx_roll() const { return _elements_[1]; };
		const double& tx_pitch() const { return _elements_[2]; };
		const double& tx_yaw() const { return _elements_[3]; };
		const double& txrx_dx() const { return _elements_[4]; };
		const double& txrx_dy() const { return _elements_[5]; };
		const double& txrx_dz() const { return _elements_[6]; };
		const double& rx_roll() const { return _elements_[7]; };
		const double& rx_pitch() const { return _elements_[8]; };
		const double& rx_yaw() const { return _elements_[9]; };


		TDEmGeometry() {};

		TDEmGeometry(const double& _tx_height, const double& _tx_roll, const double& _tx_pitch, const double& _tx_yaw, const double& _txrx_dx, const double& _txrx_dy, const double& _txrx_dz, const double& _rx_roll, const double& _rx_pitch, const double& _rx_yaw) {
			initialise(_tx_height, _tx_roll, _tx_pitch, _tx_yaw, _txrx_dx, _txrx_dy, _txrx_dz, _rx_roll, _rx_pitch, _rx_yaw);
		}

		TDEmGeometry(const double* g) {
			//const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw)
			initialise(g[0], g[1], g[2], g[3], g[4], g[5], g[6], g[7], g[8], g[9]);
		}

		TDEmGeometry(const std::vector<double> gvector) {
			for (size_t i = 0; i < NELEM; i++) {
				_elements_[i] = gvector[i];
			}
		}

		TDEmGeometry(const cBlock& b) {
			set_zero();
			for (size_t gi = 0; gi < NELEM; gi++) {
				b.getvalue(defined_names[gi], _elements_[gi]);
			}
		};

		inline static const size_t nelem() {
			return NELEM;
		}

		void initialise(const double& _tx_height, const double& _tx_roll, const double& _tx_pitch, const double& _tx_yaw, const double& _txrx_dx, const double& _txrx_dy, const double& _txrx_dz, const double& _rx_roll, const double& _rx_pitch, const double& _rx_yaw) {
			tx_height() = _tx_height;
			tx_roll() = _tx_roll; tx_pitch() = _tx_pitch; tx_yaw() = _tx_yaw;
			txrx_dx() = _txrx_dx; txrx_dy() = _txrx_dy; txrx_dz() = _txrx_dz;
			rx_roll() = _rx_roll; rx_pitch() = _rx_pitch; rx_yaw() = _rx_yaw;
		};

		double& operator[](const size_t& index) {
			return _elements_[index];
		};

		const double& operator[](const size_t& index) const {
			return _elements_[index];
		};

		double& operator[](const std::string& gname) {
			const size_t index = element_index(gname);
			return _elements_[index];
		};

		const double& operator[](const std::string& gname) const {
			const size_t index = element_index(gname);
			return _elements_[index];			
		};

		void set_zero() {
			for (size_t i = 0; i < NELEM; i++) {
				_elements_[i] = 0.0;
			}
		}

		void fillundefined(const TDEmGeometry& g) {
			for (size_t i = 0; i < NELEM; i++) {
				if (_elements_[i] == undefined) {
					_elements_[i] = g[i];
				}
			}
		};

		static size_t element_index(const std::string& name) {
			for (size_t i = 0; i < NELEM; i++) {
				if (strcasecmp(name, element_name(i)) == 0) return i;
			}
			glog.errormsg(_SRC_, "Geometry field name %s is bad\n", name.c_str());
			return 0;
		};

		static std::string element_name(const size_t& index) {
			if(index >= TDEmGeometry::NELEM) {
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
			}
			return TDEmGeometry::defined_names[index];
		};

		static std::string units(const size_t& index) {
			if (index >= TDEmGeometry::NELEM) {
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
			}
			return TDEmGeometry::defined_units[index];
		};

		static std::string description(const size_t& index) {
			if (index >= TDEmGeometry::NELEM) {
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
			}
			return TDEmGeometry::defined_descriptions[index];
		};

		static ElementType elementtype(const size_t& index) {
			if (index >= TDEmGeometry::NELEM) {
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
			}
			return TDEmGeometry::defined_elementtypes[index];
		};

		static AEM::CalculationType::Mode derivative_mode(const size_t& index) {
			if (index >= TDEmGeometry::NELEM) {
				glog.errormsg(_SRC_, "Geometry index %zu out of range\n", index);
			}
			return TDEmGeometry::defined_derivative_modes[index];
		};

		void write(std::string path) const {
			std::ofstream ofs(path);
			for (size_t i = 0; i < NELEM; i++) {
				ofs << element_name(i) << "\t" << _elements_[i] << std::endl;
			}
		}

		double txrx_dh() const {
			return std::hypot(txrx_dx(), txrx_dy());
		};

		double txrx_dv() const {
			return std::hypot(txrx_dx(), txrx_dz());
		};

		double txrx_dr() const {
			return std::sqrt(txrx_dx() * txrx_dx() + txrx_dy() * txrx_dy() + txrx_dz() * txrx_dz());
		};

		Vec3d tx_orientation(const Vec3d& tx_reference_orientation) const {
			Mat3d Rot = YPR(tx_roll(), tx_pitch(), tx_yaw());
			return Rot * tx_reference_orientation;
		}

		Vec3d txrx_separation() const {
			Vec3d v = Vec3d(txrx_dx(), txrx_dy(), txrx_dz());
			return v;
		}

		inline Mat3d inertial_to_rx_frame_rotation_matrix() const {
			// Mat3d inertial_to_rx_frame_rotation_matrix() const {
			return invYPR(rx_roll(), rx_pitch(), rx_yaw());
		};

		inline Mat3d tx_roll_derivative_matrix() const {
			// M = R P Y // dM/dr = dR/dr * P * Y
			const Mat3d dR = roll_matrix_derivative_degrees(tx_roll());
			const Mat3d P = pitch_matrix_degrees(tx_pitch());
			const Mat3d Y = yaw_matrix_degrees(tx_yaw());
			const Mat3d m = dR * P * Y;
			return m;
		};

		inline Mat3d tx_pitch_derivative_matrix() const {
			// M = R P Y // dM/dp = R * P/dp * Y
			const Mat3d R = roll_matrix_degrees(tx_roll());
			const Mat3d dP = pitch_matrix_derivative_degrees(tx_pitch());
			const Mat3d Y = yaw_matrix_degrees(tx_yaw());
			const Mat3d m = R * dP * Y;
			return m;
		};

		inline Mat3d tx_yaw_derivative_matrix() const {
			// M = R P Y // dM/dp = R * P * dY/dy
			const Mat3d R = roll_matrix_degrees(tx_roll());
			const Mat3d P = pitch_matrix_degrees(tx_pitch());
			const Mat3d dY = yaw_matrix_derivative_degrees(tx_yaw());
			const Mat3d m = R * P * dY;
			return m;
		};

		inline Mat3d dInt2RxFrame_droll() const {
			const Mat3d dR = roll_matrix_derivative_degrees(rx_roll());
			const Mat3d P = pitch_matrix_degrees(rx_pitch());
			const Mat3d Y = yaw_matrix_degrees(rx_yaw());
			const Mat3d m = (dR * P * Y).transpose();
			return m;
		};

		inline Mat3d dInt2RxFrame_dpitch() const {
			const Mat3d R = roll_matrix_degrees(rx_roll());
			const Mat3d dP = pitch_matrix_derivative_degrees(rx_pitch());
			const Mat3d Y = yaw_matrix_degrees(rx_yaw());
			const Mat3d m = (R * dP * Y).transpose();
			return m;
		};

		inline Mat3d dInt2RxFrame_dyaw() const {
			const Mat3d R = roll_matrix_degrees(rx_roll());
			const Mat3d P = pitch_matrix_degrees(rx_pitch());
			const Mat3d dY = yaw_matrix_derivative_degrees(rx_yaw());
			const Mat3d m = (R * P * dY).transpose();
			return m;
		};
	};
};