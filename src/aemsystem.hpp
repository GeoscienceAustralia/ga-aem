/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once
#include "blocklanguage.hpp"
#include "lem.hpp"
#include "layeredearthmodeller.hpp"
#include "aem_coredefs.hpp"
#include "tdemresponse.hpp"
#include "tdemgeometry.hpp"

namespace AEM {
	template <typename RT>
	class AEMSystem {

		using VectorResponse = TDEmVectorResponse<RT>;
		using Response = TDEmResponse<RT>;

	protected:
		std::string SystemName;
		cBlock STM;
		LEModeller LEM;
		Transmitter Tx;
		Receiver    Rx;
		Mat3d InertialToRxFrame;//Matrix to rotate inertial frame vector to Rx frame
		Vec3d Scale;
		Vec3d RefGeomPrimary; // For PPM systems
		WindowingScheme WindScheme;

	public:
		
		AEMSystem() {};

		AEMSystem(const fs::path& descriptorpath) {
			read_system_descriptor_file(descriptorpath);
		};

		const std::string& name() const { return SystemName; };
		const cBlock& system_descriptor_block() const { return STM; };
		LEModeller& lem() { return LEM; };

		virtual SystemType type() const = 0;
		virtual std::string type_string() const = 0;

		const size_t& nWindows() const {return WindScheme.nWindows();};

		VectorResponse forward_model_primary_field(const TDEmGeometry& G) {
			set_geometry(G);
			set_calculationtype(CMode::FM);
			VectorResponse P(nWindows());
			set_primaryfields(P, Tx.Orientation, InertialToRxFrame);
			return P;
		};

		Response forward_model(const Earth1D& E, const TDEmGeometry& G) {
			set_earth(E);
			set_geometry(G);
			setup_computations();
			set_calculationtype(CMode::FM);

			Response F(nWindows());
			set_primaryfields(F.P, Tx.Orientation, InertialToRxFrame);
			set_secondaryfields(F.S, Tx.Orientation, InertialToRxFrame);
			return F;
		};

		Response derivative(const TDEmGeometry& G, const CalculationType& calc) {
			Response D(nWindows());
			const CMode& cmode = calc.get_mode();
			if (cmode == CMode::DC || cmode == CMode::DT || cmode == CMode::DX || cmode == CMode::DY || cmode == CMode::DZ || cmode == CMode::DH){
				set_calculationtype(calc);
				set_primaryfields(D.P, Tx.Orientation, InertialToRxFrame);
				set_secondaryfields(D.S, Tx.Orientation, InertialToRxFrame);
			}
			else if(cmode == CMode::DTX_HEIGHT) {
				// This is because when H changes Z also changes
				set_calculationtype(CMode::DZ);
				set_primaryfields(D.P, Tx.Orientation, InertialToRxFrame);
				set_secondaryfields(D.S, Tx.Orientation, InertialToRxFrame);

				Response DH(nWindows());
				set_calculationtype(CMode::DH);
				set_primaryfields(DH.P, Tx.Orientation, InertialToRxFrame);
				set_secondaryfields(DH.S, Tx.Orientation, InertialToRxFrame);
				D += DH;
			}
			else if (cmode == CMode::DTX_ROLL) {
				const Mat3d& M = G.tx_roll_derivative_matrix();
				const Vec3d& txvec = M * Tx.Reference_Orientation;
				set_primaryfields(D.P, txvec, InertialToRxFrame);
				set_secondaryfields(D.S, txvec, InertialToRxFrame);
			}
			else if (cmode == CMode::DTX_PITCH) {
				const Mat3d& M = G.tx_pitch_derivative_matrix();
				const Vec3d& txvec = M * Tx.Reference_Orientation;
				set_primaryfields(D.P, txvec, InertialToRxFrame);
				set_secondaryfields(D.S, txvec, InertialToRxFrame);
			}
			else if (cmode == CMode::DTX_YAW) {
				const Mat3d& M = G.tx_yaw_derivative_matrix();
				const Vec3d& txvec = M * Tx.Reference_Orientation;
				set_primaryfields(D.P, txvec, InertialToRxFrame);
				set_secondaryfields(D.S, txvec, InertialToRxFrame);
			}
			else if (cmode == CMode::DRX_ROLL) {
				const Mat3d& rxmat = G.dInt2RxFrame_droll();
				set_primaryfields(D.P, Tx.Orientation, rxmat);
				set_secondaryfields(D.S, Tx.Orientation, rxmat);
			}
			else if (cmode == CMode::DRX_PITCH) {
				const Mat3d& rxmat = G.dInt2RxFrame_dpitch();
				set_primaryfields(D.P, Tx.Orientation, rxmat);
				set_secondaryfields(D.S, Tx.Orientation, rxmat);
			}
			else if (cmode == CMode::DRX_YAW) {
				const Mat3d& rxmat = G.dInt2RxFrame_dyaw();
				set_primaryfields(D.P, Tx.Orientation, rxmat);
				set_secondaryfields(D.S, Tx.Orientation, rxmat);
			}
			else {
				glog.errormsg(_SRC_, "Invalid derivative operation.");
			};
			return D;
		};	

	protected:
		virtual void read_system_descriptor_file(const fs::path& systemdescriptorfile) = 0;
		virtual bool is_ppm_system() const = 0;

		virtual void set_primaryfields(VectorResponse& P, const Vec3d& txvec, const Mat3d& rxmat) = 0;
		virtual void set_secondaryfields(VectorResponse& S, const Vec3d& txvec, const Mat3d& rxmat) = 0;

		void set_earth(const Earth1D& earth) {
			lem().set_earth(earth);
		};

		void set_calculationtype(const CalculationType& _calculationtype) {
			lem().set_calculationtype(_calculationtype);
		};

		void setup_computations() {
			lem().setup_computations();
		};

		void set_geometry(const TDEmGeometry& G) {
			// Set geometry inside the LE Modeller
			const Vec3d sep = G.txrx_separation();
			const double& h = G.tx_height();
			const double& x = sep.x();
			const double& y = sep.y();
			const double  z = h + sep.z();
			
			Tx.Orientation = G.tx_orientation(Tx.Reference_Orientation);
			lem().set_geometry(Tx.Orientation, h, x, y, z);
			// Set the rotation matrix for rotating vector fields to Rx frame of reference
			InertialToRxFrame = G.inertial_to_rx_frame_rotation_matrix();
		};

		/*
		void drx_roll(const TDEmGeometry& G, const VectorResponse& forward_model, VectorResponse& derivatives) const {
			const Mat3d dM = G.rx_roll_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model, derivatives);
		};

		void drx_pitch(const TDEmGeometry& G, const VectorResponse& forward_model, VectorResponse& derivatives) const {
			const Mat3d dM = G.rx_pitch_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model, derivatives);
		};

		void drx_yaw(const TDEmGeometry& G, const VectorResponse& forward_model, VectorResponse& derivatives) const {
			const Mat3d dM = G.rx_yaw_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model, derivatives);
		};

		void drx_roll(const TDEmGeometry& G, const Response& forward_model, Response& derivatives) const {
			const Mat3d dM = G.rx_roll_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model.P, derivatives.P);
			apply_rx_derivative_matrix(dM, forward_model.S, derivatives.S);
		};

		void drx_pitch(const TDEmGeometry& G, const Response& forward_model, Response& derivatives) const {
			const Mat3d dM = G.rx_pitch_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model.P, derivatives.P);
			apply_rx_derivative_matrix(dM, forward_model.S, derivatives.S);
		};

		void drx_yaw(const TDEmGeometry& G, const Response& forward_model, Response& derivatives) const {
			const Mat3d dM = G.rx_yaw_derivative_matrix();
			apply_rx_derivative_matrix(dM, forward_model.P, derivatives.P);
			apply_rx_derivative_matrix(dM, forward_model.S, derivatives.S);
		};*/

		/*
		void apply_rx_derivative_matrix(const Mat3d& dM, const VectorResponse& fields, VectorResponse& derivatives) const {
			const size_t n = fields.nWindows();
			//if (MO.NormalisationType == ModellingOptions::NormalizationType::PPM || MO.NormalisationType == ModellingOptions::NormalizationType::PPM_PEAKTOPEAK) {
			if(is_ppm_system()){
				for (size_t i = 0; i < n; i++) {
					auto ftrue = fields.get_vec3(i);
					//Must work with true field vector directions (not the PPM scaled versinn)
					ftrue[XCOMP] *= RefGeomPrimary[XCOMP];
					ftrue[YCOMP] *= RefGeomPrimary[YCOMP];
					ftrue[ZCOMP] *= RefGeomPrimary[ZCOMP];
					derivatives.set_vec3(i, dM * ftrue);
					//Convert back to PPMS
					derivatives[XCOMP][i] /= RefGeomPrimary[XCOMP];
					derivatives[YCOMP][i] /= RefGeomPrimary[YCOMP];
					derivatives[ZCOMP][i] /= RefGeomPrimary[ZCOMP];
				}
			}
			else {
				for (size_t wi = 0; wi < n; wi++) {
					derivatives.set_vec3(wi, dM * fields.get_vec3(wi));
				}
			}
		};*/

	};
};

