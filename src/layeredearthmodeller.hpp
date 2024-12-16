/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#pragma once

#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

#include "logger.hpp"
#include "general_constants.hpp"
#include "eigen_utils.hpp"
#include "calculation_type.hpp"
#include "earth1d.hpp"

//Formulation mainly from the book 
//Geo-Electromagnetism, Wait James, R. Academic Press 1982.

namespace AEM {
namespace LEM2 {
	using cdouble = std::complex<double>;
	using cvector = std::vector<std::complex<double>>;
	using CalculationType = CT::CalculationType;
	using CMode = CT::CalculationType::Mode;
	using PropogationMatrix = Eigen::Matrix2cd;
	constexpr double DefaultLowerFractionalWidth = 4.44;
	constexpr double DefaultUpperFractionalWidth = 1.84;

	struct AbscissaLayerNode {
		cdouble   U;
		cdouble   Exp2UT;
		PropogationMatrix Matrix;
		PropogationMatrix preMatrix;
		PropogationMatrix postMatrix;
	};

	class AbscissaNode {

	public:
		double Lambda = 0.0;
		double Lambda2 = 0.0;
		double Lambda3 = 0.0;
		double Lambda4 = 0.0;
		double Lambda_r = 0.0;
		double j0Lambda_r = 0.0;
		double j1Lambda_r = 0.0;
		std::vector<AbscissaLayerNode> Layer;
		PropogationMatrix P_Full;
		cdouble P21onP11;

		double LoopFactor(const double& ModellingLoopRadius) const {
			double factor = 1.0;
			if (ModellingLoopRadius > 0.0) {
				const double lambda_a = Lambda * ModellingLoopRadius;
				factor = 2.0 * std::cyl_bessel_j(1, lambda_a) / lambda_a;
			}
			return factor;
		}
	};

	class LESingleFrequencyModeller {

	private:

		double x = 0.0, y = 0.0, z = 0.0, h = 0.0;
		double x2 = 0, y2 = 0, x4 = 0, y4 = 0;
		double zh = 0, zh2 = 0;
		double r = 0, r2 = 0, r3 = 0, r5 = 0, r4 = 0;
		double R = 0, R2 = 0, R5 = 0, R7 = 0;

		bool geometrychanged;
		std::vector<AbscissaNode> Abscissa;
		std::shared_ptr<Earth1D> EarthPtr;
		bool earthchanged;

		Vec3cd ForwardModel;

		//Hankle Stuff
		double meanlog10conductivity;
		double LowerFractionalWidth = LEM2::DefaultLowerFractionalWidth;
		double UpperFractionalWidth = LEM2::DefaultUpperFractionalWidth;

	public:
		double ModellingLoopRadius = 0.0; //dipole by default
		double Frequency;
		double Omega;
		double MuZeroOmega;
		cdouble iMuZeroOmega;
		double ApproximateHalfspace;
		double PeakLambda;
		double LowerBound;
		double UpperBound;
		double AbscissaSpacing;

		LESingleFrequencyModeller() { };
		
		size_t nAbscissa1() const { return Abscissa.size(); }
		
		size_t nLayers() const { return EarthPtr->nlayers(); }
		
		void initialise(const double& frequency, const size_t& numabscissa, const double& modelling_loop_radius, const double _LowerFractionalWidth = LEM2::DefaultLowerFractionalWidth, const double _UpperFractionalWidth = LEM2::DefaultUpperFractionalWidth) {
			set_frequency(frequency);
			Abscissa.resize(numabscissa);
			ModellingLoopRadius = modelling_loop_radius;
			LowerFractionalWidth = _LowerFractionalWidth;
			UpperFractionalWidth = _UpperFractionalWidth;
		};

		void setupcomputations() {
			if (earthchanged || geometrychanged) {
				ForwardModel[0] = std::numeric_limits<double>::max();
				ForwardModel[1] = std::numeric_limits<double>::max();
				ForwardModel[2] = std::numeric_limits<double>::max();
				setfrequencyabscissalayers();
				earthchanged = false;
				geometrychanged = false;
			}
			else{
				glog.warningmsg(_SRC_, "Redoing setupcomputations without earth or geometry change.");
			}
		}

		void set_earth_ptr(const std::shared_ptr<Earth1D>& earthptr, const double _meanlog10conductivity) {
			EarthPtr = earthptr;
			meanlog10conductivity = _meanlog10conductivity;
			earthchanged = true;
			const size_t na = nAbscissa1();
			for (size_t ai = 0; ai < na; ai++) {
				Abscissa[ai].Layer.resize(EarthPtr->nlayers());
			}
		}

		void setxyzh(const double& _x, const double& _y, const double& _z, const double& _h)
		{
			bool update_r = false;
			if (x != _x) {
				x = _x;
				x2 = x * x;
				x4 = x2 * x2;
				update_r = true;
			}

			if (_y != y) {
				y = _y;
				y2 = y * y;
				y4 = y2 * y2;
				update_r = true;
			}

			if (update_r) {
				r2 = x2 + y2;
				r = std::sqrt(r2);
				r3 = r2 * r;
				r4 = r2 * r2;
				r5 = r4 * r;
			}

			bool update_zh = false;
			if (z != _z) {
				z = _z;
				update_zh = true;
			}

			if (h != _h) {
				h = _h;
				update_zh = true;
			}

			if (update_zh) {
				zh = z - h;
				zh2 = zh * zh;
			}

			if (update_r || update_zh) {
				R2 = (x2 + y2 + zh2);
				R = std::sqrt(R2);
				R5 = R2 * R2 * R;
				R7 = R5 * R2;
				geometrychanged = true;
			}
			else {
				geometrychanged = false;
			}
		}

	private:

		void set_frequency(const double& frequency) {
			Frequency = frequency;
			Omega = TWOPI<double> *frequency;
			MuZeroOmega = MUZERO<double> *Omega;
			iMuZeroOmega = cdouble(0.0, MuZeroOmega);
		};

		void setfrequencyabscissalayers() {
			const size_t nl = nLayers();
			const size_t na = nAbscissa1();
			setintegrationnodes();
			for (size_t ai = 0; ai < na; ai++) {
				for (size_t li = 0; li < nl; li++) {
					double gamma2 = EarthPtr->conductivity[li] * MuZeroOmega;
					cdouble u = std::sqrt(cdouble(Abscissa[ai].Lambda2, gamma2));
					Abscissa[ai].Layer[li].U = u;
					if (li < nl - 1) Abscissa[ai].Layer[li].Exp2UT = exp(-2.0 * u * EarthPtr->thickness[li]);
				}
				setlayermatrices(ai);
				setpmatrix(ai);
			}
		};

		double approximatehalfspace() const {
			if (nLayers() == 1) return EarthPtr->conductivity[0];

			//peak lambda
			double   peaklambda = std::sqrt(MuZeroOmega * meanlog10conductivity / 4.0);
			const cdouble  rz = rzero_recursive(peaklambda);
			cdouble  v = (1.0 + rz);
			v = iMuZeroOmega * v * v;
			cdouble a = (-4.0 * peaklambda * peaklambda * rz / v);
			return a.real();
		}

		cdouble rzero_recursive(const double lambda) const {
			//Wait's recursive formulation
			const int nl = (int)nLayers();
			double lambda2 = lambda * lambda;
			double muzeroomega = MuZeroOmega;
			cdouble imuzeroomega(0.0, muzeroomega);

			double gamma2 = muzeroomega * EarthPtr->conductivity[nl - 1];
			cdouble u = std::sqrt(cdouble(lambda2, gamma2));
			cdouble y = u / imuzeroomega;
			for (int li = nl - 2; li >= 0; li--) {
				gamma2 = muzeroomega * EarthPtr->conductivity[li];
				u = std::sqrt(cdouble(lambda2, gamma2));
				const cdouble Nn = u / imuzeroomega;

				const cdouble v1 = u * EarthPtr->thickness[li];

				//Expand - unstable    	
				//tanh(v) = (1.0 - v4)/(1.0 + v4 + 2.0*v2);
				const cdouble v2 = std::exp(-2.0 * v1);
				const cdouble v4 = std::exp(-4.0 * v1);
				const cdouble top = 1.0 - v4;
				const cdouble bot = 1.0 + v4 + 2.0 * v2;
				const cdouble tanhv = top / bot;

				y = Nn * (y + Nn * tanhv) / (Nn + y * tanhv);
			}

			cdouble N0(0.0, -lambda / (muzeroomega));  // minus because dividing by i.muzero.omega
			return (N0 - y) / (N0 + y);

		};

		inline cdouble rzero_propogationmatrix(const size_t ai) {
			//Oldenberg's propogation matrix formulation
			setlayermatrices(ai);
			setpmatrix(ai);
			return Abscissa[ai].P21onP11;
		};

		void setlayermatrices(const size_t ai) {
			cdouble e, eh, e1, e2;
			AbscissaNode& A = Abscissa[ai];

			//M1  
			e = A.Layer[0].U / A.Lambda;
			eh = e / 2.0;
			e1 = 0.5 + eh;
			e2 = 0.5 - eh;

			PropogationMatrix M;
			M(0, 0) = e1;
			M(0, 1) = e2;
			M(1, 0) = e2;
			M(1, 1) = e1;
			A.Layer[0].Matrix = M;

			const size_t nl = nLayers();
			for (size_t li = 1; li < nl; li++) {
				//assumes all pearmabilities are muzero
				e = A.Layer[li].U / A.Layer[li - 1].U;
				eh = e / 2.0;
				e1 = 0.5 + eh;
				e2 = 0.5 - eh;
				M(0, 0) = e1;
				M(0, 1) = e2;
				M(1, 0) = e2 * A.Layer[li - 1].Exp2UT;
				M(1, 1) = e1 * A.Layer[li - 1].Exp2UT;

				A.Layer[li].Matrix = M;
			}

			//Set Prematrices - prematrix for layer 0 does not apply
			if (nl > 1) A.Layer[1].preMatrix = A.Layer[0].Matrix;
			for (size_t li = 2; li < nl; li++) {
				A.Layer[li].preMatrix = A.Layer[li - 1].preMatrix * A.Layer[li - 1].Matrix;
			}

			//Set Postmatrices - postmatrix for layer N-1 does not apply  
			if (nl > 1) {
				A.Layer[nl - 2].postMatrix = A.Layer[nl - 1].Matrix;
			}

			for (int li = (int)nl - 3; li >= 0; li--) {
				A.Layer[li].postMatrix = A.Layer[(size_t)li + 1].Matrix * A.Layer[(size_t)li + 1].postMatrix;
			}

		};

		void setpmatrix(const size_t ai) {
			const size_t nl = nLayers();
			AbscissaNode& A = Abscissa[ai];
			//Set Full matrix
			if (nl == 1) {
				A.P_Full = A.Layer[0].Matrix;
			}
			else {
				A.P_Full = Abscissa[ai].Layer[nl - 1].preMatrix * Abscissa[ai].Layer[nl - 1].Matrix;
			}
			A.P21onP11 = A.P_Full(1, 0) / A.P_Full(0, 0);
		}

		PropogationMatrix dMjdCj(const size_t ai, const size_t li) const {
			PropogationMatrix m;
			if (li == 0) {
				const cdouble a = iMuZeroOmega / (4.0 * Abscissa[ai].Lambda * Abscissa[ai].Layer[li].U);
				m(0, 0) = a;
				m(0, 1) = -a;
				m(1, 0) = -a;
				m(1, 1) = a;
				return m;
			}
			else {
				const cdouble a = iMuZeroOmega / (4.0 * Abscissa[ai].Layer[li - 1].U * Abscissa[ai].Layer[li].U);
				const cdouble ae = a * Abscissa[ai].Layer[li - 1].Exp2UT;
				m(0, 0) = a;
				m(0, 1) = -a;
				m(1, 0) = -ae;
				m(1, 1) = ae;
				return m;
			}
		}

		PropogationMatrix dMjplus1dCj(const size_t ai, const size_t li) const {
			cdouble duds = iMuZeroOmega / (2.0 * Abscissa[ai].Layer[li].U);
			cdouble y = Abscissa[ai].Layer[li + 1].U / Abscissa[ai].Layer[li].U;

			cdouble dydu = -y / Abscissa[ai].Layer[li].U;
			cdouble dyds = dydu * duds;

			cdouble v = Abscissa[ai].Layer[li].Exp2UT;
			cdouble dvdu = -2.0 * EarthPtr->thickness[li] * Abscissa[ai].Layer[li].Exp2UT;

			cdouble dvds = dvdu * duds;

			cdouble ydvds = y * dvds;
			cdouble vdyds = v * dyds;


			PropogationMatrix M;
			M(0, 0) = 0.5 * dyds;
			M(0, 1) = -M(0, 0);
			M(1, 0) = 0.5 * (dvds - (ydvds + vdyds));
			M(1, 1) = 0.5 * (dvds + (ydvds + vdyds));

			return M;

		}

		PropogationMatrix dPdCj(const size_t ai, const size_t li) const {
			const AbscissaNode& A = Abscissa[ai];
			const size_t nl = nLayers();
			//One layer case
			if (nl == 1) return dMjdCj(ai, li);

			//Not last layer
			if (li < nl - 1) {
				PropogationMatrix c = dMjdCj(ai, li) * A.Layer[li + 1].Matrix + A.Layer[li].Matrix * dMjplus1dCj(ai, li);
				//First layer case
				if (li == 0) {
					if (nl == 2) return c;
					return c * A.Layer[li + 1].postMatrix;
				}
				else if (li == nl - 2) {
					return A.Layer[li].preMatrix * c;
				}
				else {
					return A.Layer[li].preMatrix * c * A.Layer[li + 1].postMatrix;
				}
			}
			//Last layer
			else {
				return A.Layer[li].preMatrix * dMjdCj(ai, li);
			}

		}

		cdouble dP21onP11dCj(const size_t ai, const size_t li) const {
			PropogationMatrix m = dPdCj(ai, li);
			return m(1, 0) / Abscissa[ai].P_Full(0, 0) - m(0, 0) * Abscissa[ai].P21onP11 / Abscissa[ai].P_Full(0, 0);
		}

		PropogationMatrix dMjplus1dTj(const size_t ai, const size_t li) const {
			cdouble y = Abscissa[ai].Layer[li + 1].U / Abscissa[ai].Layer[li].U;
			cdouble dvdt = -2.0 * Abscissa[ai].Layer[li].U * Abscissa[ai].Layer[li].Exp2UT;

			PropogationMatrix m;
			m(0, 0) = 0.0;
			m(0, 1) = 0.0;
			m(1, 0) = 0.5 * (1.0 - y) * dvdt;
			m(1, 1) = 0.5 * (1.0 + y) * dvdt;
			return m;

		}

		PropogationMatrix dPdTj(const size_t ai, const size_t li) const {
			const size_t nl = nLayers();
			const AbscissaNode& A = Abscissa[ai];

			//One layer case
			if (nl == 1) return PropogationMatrix::Zero();
			PropogationMatrix b = A.Layer[li].Matrix * dMjplus1dTj(ai, li);

			//First layer case
			if (li == 0) {
				if (nl == 2) return b;
				return b * A.Layer[li + 1].postMatrix;
			}
			else if (li > 0 && li < nl - 1) {
				PropogationMatrix d = A.Layer[li].preMatrix * b;
				if (li == nl - 2) {
					return d;
				}
				return d * A.Layer[li + 1].postMatrix;
			}
			//Last layer case
			else {
				glog.errormsg(_SRC_, "Zero thickness derivative for halfspace layer\n");
				return PropogationMatrix::Zero();
			}
		}

		cdouble dP21onP11dTj(const size_t ai, const size_t li) const {
			PropogationMatrix m = dPdTj(ai, li);
			return m(1, 0) / Abscissa[ai].P_Full(0, 0) - m(0, 0) * Abscissa[ai].P21onP11 / Abscissa[ai].P_Full(0, 0);
		}

		void setintegrationnodes() {
			const size_t na = nAbscissa1();
			double peak_exp2 = 2.0 / (z + h);
			double peak_exp3 = 3.0 / (z + h);

			ApproximateHalfspace = approximatehalfspace();
			PeakLambda = std::sqrt(MuZeroOmega * ApproximateHalfspace / 4.0);

			double lp = std::log(std::min(PeakLambda, peak_exp2));
			double up = std::log(std::max(PeakLambda, peak_exp3));

			LowerBound = lp - LowerFractionalWidth;
			UpperBound = up + UpperFractionalWidth;

			AbscissaSpacing = (UpperBound - LowerBound) / (double)(na - 1);

			double lambda;
			double loglambda = LowerBound;
			for (size_t ai = 0; ai < na; ai++) {
				AbscissaNode& A = Abscissa[ai];
				lambda = std::exp(loglambda);
				A.Lambda = lambda;
				A.Lambda2 = A.Lambda * lambda;
				A.Lambda3 = A.Lambda2 * lambda;
				A.Lambda4 = A.Lambda3 * lambda;
				A.Lambda_r = A.Lambda * r;
				A.j0Lambda_r = std::cyl_bessel_j(0, A.Lambda_r);
				A.j1Lambda_r = std::cyl_bessel_j(1, A.Lambda_r);
				loglambda += AbscissaSpacing;
			}
		}

		Vec3cd compute_integrand(const size_t& ai, const CalculationType& calculationtype) const {
			const AbscissaNode& A = Abscissa[ai];

			const double& lambdar = A.Lambda_r;
			const double& j0 = A.j0Lambda_r;
			const double& j1 = A.j1Lambda_r;
			const double& e = std::exp(-(z + h) * A.Lambda);
			const double& l2e = A.Lambda2 * e;
			const double& l3e = A.Lambda3 * e;
			const double& l4e = A.Lambda4 * e;

			double k0, k1, k2;
			k0 = k1 = k2 = std::numeric_limits<double>::max();
			cdouble earthkernel;

			const size_t& derivativelayer = calculationtype.get_layer();
			switch (calculationtype.get_mode()) {
			case CMode::FM:
				earthkernel = -A.P21onP11;
				k0 = l3e * j0;
				k1 = l3e * j1;
				k2 = l2e * j1;
				break;
			case CMode::DX:
				earthkernel = -A.P21onP11;
				k0 = -l4e * j1 * x / r;
				k1 = l4e * (j0 - j1 / lambdar) * x / r;
				k2 = l3e * (j0 - j1 / lambdar) * x / r;
				break;
			case CMode::DY:
				earthkernel = -A.P21onP11;
				k0 = -l4e * j1 * y / r;
				k1 = l4e * (j0 - j1 / lambdar) * y / r;
				k2 = l3e * (j0 - j1 / lambdar) * y / r;
				break;
			case CMode::DZ:
				earthkernel = -A.P21onP11;
				k0 = -l4e * j0;
				k1 = -l4e * j1;
				k2 = -l3e * j1;
				break;
			case CMode::DH:
				earthkernel = -A.P21onP11;;
				k0 = -l4e * j0;
				k1 = -l4e * j1;
				k2 = -l3e * j1;
				break;
			case CMode::DC:
				earthkernel = -dP21onP11dCj(ai, derivativelayer);
				k0 = l3e * j0;
				k1 = l3e * j1;
				k2 = l2e * j1;
				break;
			case CMode::DT:
				earthkernel = -dP21onP11dTj(ai, derivativelayer);
				k0 = l3e * j0;
				k1 = l3e * j1;
				k2 = l2e * j1;
				break;
			default:
				glog.errormsg(_SRC_, "Error: compute_integrand() unknown calculation type %c\n", calculationtype);
			}

			const double loopfactor = A.LoopFactor(ModellingLoopRadius);

			Vec3cd integrand;
			integrand[0] = earthkernel * (k0 * loopfactor);
			integrand[1] = earthkernel * (k1 * loopfactor);
			integrand[2] = earthkernel * (k2 * loopfactor);
			return integrand;
		};

		Vec3cd integrate_trapezoidal(const CalculationType& calculationtype) const {
			const size_t na = nAbscissa1();
			Vec3cd sum = Vec3cd::Zero();
			// First abscissa
			sum = 0.5 * compute_integrand(0, calculationtype);
			// Cenral abscissas
			for (size_t ai = 1; ai < na - 1; ai++) {
				sum += compute_integrand(ai, calculationtype);
			}
			// Last Abscissa
			sum += 0.5 * compute_integrand(na - 1, calculationtype);

			// Scale for abscissa spacing
			sum *= AbscissaSpacing;
			return sum;
		}

		// Tensors
		Mat3d PTFM() const {
			Mat3d T;
			T(0, 0) = (3.0 * x2 - R2) / R5;
			T(0, 1) = 3.0 * x * y / R5;
			T(0, 2) = 3.0 * x * zh / R5;

			//Note error in Fitterman and Yin paper should no be minus sign at element 2,1
			T(1, 0) = T(0, 1);
			T(1, 1) = (3.0 * y2 - R2) / R5;
			T(1, 2) = 3.0 * y * zh / R5;

			T(2, 0) = T(0, 2);
			T(2, 1) = T(1, 2);
			T(2, 2) = (3.0 * zh2 - R2) / R5;
			return -ONEONFOURPI<double> * T;
		};

		Mat3d dPTdX() const {
			Mat3d T;
			T(0, 0) = -3.0 * x * (2.0 * x2 - 3.0 * y2 - 3.0 * zh2) / R7;
			T(0, 1) = -3.0 * y * (4.0 * x2 - y2 - zh2) / R7;
			T(0, 2) = -3.0 * zh * (4.0 * x2 - y2 - zh2) / R7;
			T(1, 0) = T(0, 1);
			T(1, 1) = 3.0 * x * (x2 - 4.0 * y2 + zh2) / R7;
			T(1, 2) = -15.0 * y * zh / R7 * x;
			T(2, 0) = T(0, 2);
			T(2, 1) = T(1, 2);
			T(2, 2) = 3.0 * x * (x2 + y2 - 4.0 * zh2) / R7;
			return -ONEONFOURPI<double> * T;
		};

		Mat3d dPTdY() const {
			Mat3d T;
			T(0, 0) = -3.0 * y * (4.0 * x2 - y2 - zh2) / R7;
			T(0, 1) = 3.0 * x * (x2 - 4.0 * y2 + zh2) / R7;
			T(0, 2) = -15.0 * x * zh / R7 * y;
			T(1, 0) = T(0, 1);
			T(1, 1) = 3.0 * y * (3.0 * x2 - 2.0 * y2 + 3.0 * zh2) / R7;
			T(1, 2) = 3.0 * zh * (x2 - 4.0 * y2 + zh2) / R7;
			T(2, 0) = T(0, 2);
			T(2, 1) = T(1, 2);
			T(2, 2) = 3.0 * y * (x2 + y2 - 4.0 * zh2) / R7;
			return -ONEONFOURPI<double> *T;
		};

		Mat3d dPTdZ() const {
			Mat3d T;
			T(0, 0) = -3.0 * zh * (4.0 * x2 - y2 - zh2) / R7;
			T(0, 1) = -15.0 * x * y / R7 * zh;
			T(0, 2) = 3.0 * x * (x2 + y2 - 4.0 * zh2) / R7;
			T(1, 0) = T(0, 1);
			T(1, 1) = 3.0 * zh * (x2 - 4.0 * y2 + zh2) / R7;
			T(1, 2) = 3.0 * y * (x2 + y2 - 4.0 * zh2) / R7;
			T(2, 0) = T(0, 2);
			T(2, 1) = T(1, 2);
			T(2, 2) = 3.0 * zh * (3.0 * x2 + 3.0 * y2 - 2.0 * zh2) / R7;
			return -ONEONFOURPI<double> *T;
		};

		Mat3d dPTdH() const {
			Mat3d m;
			//of course this is just minus d/dZ		
			m(0, 0) = 3.0 * zh * (4.0 * x2 - y2 - zh2) / R7;
			m(0, 1) = 15.0 * x * y / R7 * zh;
			m(0, 2) = -3.0 * x * (x2 + y2 - 4.0 * zh2) / R7;
			m(1, 0) = m(0, 1);
			m(1, 1) = -3.0 * zh * (x2 - 4.0 * y2 + zh2) / R7;
			m(1, 2) = -3.0 * y * (x2 + y2 - 4.0 * zh2) / R7;
			m(2, 0) = m(0, 2);
			m(2, 1) = m(1, 2);
			m(2, 2) = -3.0 * zh * (3.0 * x2 + 3.0 * y2 - 2.0 * zh2) / R7;
			return -ONEONFOURPI<double> *m;
		};

		Mat3cd STFM(const Vec3cd& FM) const {
			const cdouble& T0 = FM[0];
			const cdouble& T1 = FM[1];
			const cdouble& T2 = FM[2];
			
			Mat3cd m;
			m(0, 0) = ((x2 / r2 - y2 / r2) * T2 / r - T0 * x2 / r2);
			m(0, 1) = (x * y / r2) * (2.0 * T2 / r - T0);
			m(0, 2) = (-x / r) * T1;

			//Note error in Fitterman and Yin paper should not be minus sign at element 2,1
			m(1, 0) = m(0, 1);
			m(1, 1) = ((y2 / r2 - x2 / r2) * T2 / r - T0 * y2 / r2);
			m(1, 2) = (-y / r) * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> *m;
		};

		Mat3cd dSTdC(const Vec3cd& dC) const {
			const cdouble& T0 = dC[0];
			const cdouble& T1 = dC[1];
			const cdouble& T2 = dC[2];

			Mat3cd m;
			m(0, 0) = ((x2 / r2 - y2 / r2) * T2 / r - T0 * x2 / r2);
			m(0, 1) = (x * y / r2) * (2.0 * T2 / r - T0);
			m(0, 2) = (-x / r) * T1;

			m(1, 0) = m(0, 1);
			m(1, 1) = ((y2 / r2 - x2 / r2) * T2 / r - T0 * y2 / r2);
			m(1, 2) = (-y / r) * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> *m;
		};

		Mat3cd dSTdT(const Vec3cd& dT) const {
			const cdouble& T0 = dT[0];
			const cdouble& T1 = dT[1];
			const cdouble& T2 = dT[2];
			
			Mat3cd m;
			m(0, 0) = ((x2 / r2 - y2 / r2) * T2 / r - T0 * x2 / r2);
			m(0, 1) = (x * y / r2) * (2.0 * T2 / r - T0);
			m(0, 2) = (-x / r) * T1;

			m(1, 0) = m(0, 1);
			m(1, 1) = ((y2 / r2 - x2 / r2) * T2 / r - T0 * y2 / r2);
			m(1, 2) = (-y / r) * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> *m;
		};

		Mat3cd dSTdX(const Vec3cd& dX) const {
			const cdouble& T0 = dX[0];
			const cdouble& T1 = dX[1];
			const cdouble& T2 = dX[2];
			const cdouble& T0FM = ForwardModel[0];
			const cdouble& T1FM = ForwardModel[1];
			const cdouble& T2FM = ForwardModel[2];

			Mat3cd m;
			m(0, 0) = (x4 * T2 - T0 * x4 * r - T2FM * x2 * x - T0 * x2 * r * y2 - 2.0 * T0FM * x * r * y2 + 5.0 * x * y2 * T2FM - y4 * T2) / r5;
			m(0, 1) = 2.0 * y / r3 * T2FM - y / r2 * T0FM - 6.0 * x2 * y / r5 * T2FM + 2.0 * x2 * y / r4 * T0FM + 2.0 * x * y / r3 * T2 - x * y / r2 * T0;
			m(0, 2) = -(T1FM * y2 + x2 * x * T1 + x * T1 * y2) / r3;

			m(1, 0) = m(0, 1);
			m(1, 1) = -(T2 * x4 - T2FM * x2 * x + T0 * y2 * r * x2 - 2.0 * T0FM * y2 * x * r + 5.0 * x * y2 * T2FM + T0 * y4 * r - y4 * T2) / r5;
			m(1, 2) = y / r3 * T1FM * x - y / r * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> * m;
		};

		Mat3cd dSTdY(const Vec3cd& dY) const {
			const cdouble& T0 = dY[0];
			const cdouble& T1 = dY[1];
			const cdouble& T2 = dY[2];
			const cdouble& T0FM = ForwardModel[0];
			const cdouble& T1FM = ForwardModel[1];
			const cdouble& T2FM = ForwardModel[2];

			Mat3cd m;
			m(0, 0) = (x4 * T2 - T0 * x4 * r + 2.0 * T0FM * x2 * y * r - T0 * x2 * r * y2 - 5.0 * x2 * y * T2FM + y2 * y * T2FM - y4 * T2) / r5;
			m(0, 1) = 2.0 * x / r3 * T2FM - x / r2 * T0FM - 6.0 * x * y2 / r5 * T2FM + 2.0 * x * y2 / r4 * T0FM + 2.0 * x * y / r3 * T2 - x * y / r2 * T0;
			m(0, 2) = x / r3 * T1FM * y - x / r * T1;

			m(1, 0) = m(0, 1);
			m(1, 1) = -(T2 * x4 + T0 * y2 * r * x2 + 2.0 * T0 * y * r * x2 - 5.0 * x2 * y * T2FM + T0 * y4 * r + y2 * y * T2 - y4 * T2) / r5;
			m(1, 2) = -(T1 * x2 + y * T1 * x2 + y2 * y * T1) / r3;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> * m;
		};

		Mat3cd dSTdZ(const Vec3cd& dZ) const {
			const cdouble& T0 = dZ[0];
			const cdouble& T1 = dZ[1];
			const cdouble& T2 = dZ[2];

			Mat3cd m;
			m(0, 0) = ((x2 / r2 - y2 / r2) * T2 / r - T0 * x2 / r2);
			m(0, 1) = x * y / r2 * (2.0 * T2 / r - T0);
			m(0, 2) = -x / r * T1;

			m(1, 0) = m(0, 1);
			m(1, 1) = ((y2 / r2 - x2 / r2) * T2 / r - T0 * y2 / r2);
			m(1, 2) = -y / r * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> *m;
		};

		Mat3cd dSTdH(const Vec3cd& dH) const {
			const cdouble& T0 = dH[0];
			const cdouble& T1 = dH[1];
			const cdouble& T2 = dH[2];

			Mat3cd m;
			m(0, 0) = ((x2 / r2 - y2 / r2) * T2 / r - T0 * x2 / r2);
			m(0, 1) = x * y / r2 * (2.0 * T2 / r - T0);
			m(0, 2) = -x / r * T1;

			m(1, 0) = m(0, 1);
			m(1, 1) = ((y2 / r2 - x2 / r2) * T2 / r - T0 * y2 / r2);
			m(1, 2) = -y / r * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> * m;
		};

		Mat3d PrimaryTensor(const CalculationType& calculationtype) const {
			switch (calculationtype.get_mode()) {
			case CMode::FM: return PTFM();
			case CMode::DX: return dPTdX();
			case CMode::DY: return dPTdY();
			case CMode::DZ: return dPTdZ();
			case CMode::DH: return dPTdH();
			case CMode::DC: return Mat3d::Zero();
			case CMode::DT: return Mat3d::Zero();
			default:
				glog.errormsg(_SRC_, "Unknown calculation type %c\n", calculationtype);
			}
		};

		Mat3cd SecondaryTensor(const CalculationType& calculationtype) {
			const Vec3cd integral = integrate_trapezoidal(calculationtype);

			switch (calculationtype.get_mode()) {
			case CMode::FM: {
							// Save the forward result for the dX and dY derivative calculations
							ForwardModel = integral;
							return STFM(integral);
							};
			case CMode::DC: return dSTdC(integral);
			case CMode::DT: return dSTdT(integral);
			case CMode::DX: return dSTdX(integral);
			case CMode::DY: return dSTdY(integral);
			case CMode::DZ: return dSTdZ(integral);
			case CMode::DH: return dSTdH(integral);
			default:
				glog.errormsg(_SRC_, "Unknown calculation type %s\n", calculationtype.string().c_str());
			}
		};

	public:

		Vec3d primary_inertial_frame(const CalculationType& calculationtype, const Vec3d& txdir) const {
			return PrimaryTensor(calculationtype) * txdir;
		}

		Vec3cd secondary_inertial_frame(const CalculationType& calculationtype, const Vec3d& txdir) {
			//std::cout << txdir << std::endl;
			//std::cout << SecondaryTensor(calculationtype) << std::endl;
			return SecondaryTensor(calculationtype) * txdir;
		}

		double primary(const Vec3d& txdir, const Vec3d& rxdir)
		{
			Vec3d v = primary_inertial_frame(CMode::FM, txdir);
			return v.dot(rxdir);
		}

		cdouble secondary(const Vec3d& txdir, const Vec3d& rxdir) {
			Vec3cd v = secondary_inertial_frame(CMode::FM, txdir);
			return v.dot(rxdir);
		}

		double  dp(const CalculationType& calculationtype, const Vec3d& txdir, const Vec3d& rxdir)
		{
			Vec3d v = primary_inertial_frame(calculationtype, txdir);
			return v.dot(rxdir);
		}

		cdouble ds(const CalculationType& calculationtype, const Vec3d& txdir, const Vec3d& rxdir)
		{
			Vec3cd v = secondary_inertial_frame(calculationtype, txdir);
			return v.dot(rxdir);
		}

		cdouble dsdx(const Vec3d& txdir, const Vec3d& rxdir)
		{
			Vec3cd sf = secondary_inertial_frame(CMode::DX, txdir);
			return sf.dot(rxdir);
		}

		cdouble dsdy(const Vec3d& txdir, const Vec3d& rxdir)
		{
			Vec3cd sf = secondary_inertial_frame(CMode::DY, txdir);
			return sf.dot(rxdir);
		}

		cdouble dsdz(const Vec3d& txdir, const Vec3d& rxdir)
		{
			Vec3cd sf = secondary_inertial_frame(CMode::DZ, txdir);
			return sf.dot(rxdir);
		}

		cdouble dsdh(const Vec3d& txdir, const Vec3d& rxdir)
		{
			Vec3cd sf = secondary_inertial_frame(CMode::DH, txdir);
			return sf.dot(rxdir);
		}

		cdouble dsdc(const size_t dlayer, const Vec3d& txdir, const Vec3d& rxdir)
		{
			const CalculationType calct(CMode::DC, dlayer);
			Vec3cd sf = secondary_inertial_frame(calct, txdir);
			return sf.dot(rxdir);
		}

		cdouble dsdt(const size_t dlayer, const Vec3d& txdir, const Vec3d& rxdir)
		{
			const CalculationType calct(CMode::DT, dlayer);
			Vec3cd sf = secondary_inertial_frame(calct, txdir);
			return sf.dot(rxdir);
		}

		cdouble ppm(const Vec3d& txdir, const Vec3d& rxdir)
		{
			double  pf = primary(txdir, rxdir);
			cdouble sf = secondary(txdir, rxdir);

			sf = std::complex<double>(sf.real(), -sf.imag());
			return 1.0e6 * (sf / pf);
		}

		cdouble dppm(const CalculationType& calculationtype, const size_t& derivativelayer, const Vec3d& txdir, const Vec3d& rxdir)
		{
			double  pf = primary(txdir, rxdir);
			cdouble dsf = ds(calculationtype, txdir, rxdir);
			dsf = std::complex<double>(dsf.real(), -dsf.imag());
			return 1.0e6 * (dsf / pf);
		}

	};

	class LEModeller {

	private:
		double ModellingLoopRadius=0;
		Vec3d Source_Orientation;
		CalculationType calculationtype;
		std::shared_ptr<Earth1D> EarthPtr;

	public:
		std::vector<LESingleFrequencyModeller> FM;

		LEModeller() {};

		size_t nFrequencies() const { return FM.size(); };
		
		size_t nLayers() const { 
			// Todo
			return FM[0].nLayers(); 
		};

		void initialise(const std::vector<double>& discrete_frequencies, const size_t& numabscissa, const double& modelling_loop_radius) {
			const size_t nf = discrete_frequencies.size();
			FM.resize(nf);
			for (size_t fi = 0; fi < nf; fi++) {
				FM[fi].initialise(discrete_frequencies[fi], numabscissa, modelling_loop_radius);
			}
		}

		void set_earth(const Earth1D& earth) {
			EarthPtr = std::make_shared<Earth1D>(earth);
			const double meanlog10conductivity = EarthPtr->mean_weighted_conductivity_log10_calculation();
			const size_t nf = nFrequencies();
			for (size_t i = 0; i < nf; i++) {
				FM[i].set_earth_ptr(EarthPtr, meanlog10conductivity);
			}
		};
		
		void set_geometry(const Vec3d& source_orientation, double h, double x, double y, double z) {
			Source_Orientation = source_orientation;
			const size_t nf = nFrequencies();
			for (size_t i = 0; i < nf; i++) {
				FM[i].setxyzh(x, y, z, h);
			}
		};

		void setup_computations() {
			const size_t nf = nFrequencies();
			for (size_t i = 0; i < nf; i++) {
				FM[i].setupcomputations();
			}
		};

		void set_calculationtype(const CalculationType& _calculationtype) {
			calculationtype = _calculationtype;
		};

		const CalculationType::Mode& cmode() const {
			return calculationtype.get_mode();
		};

		Vec3d primaryfield_inertial() {
			Vec3d pf = FM[0].primary_inertial_frame(calculationtype, Source_Orientation);
			//std::cout << pf << std::endl;
			return pf;
		};

		Vec3cd secondaryfield_inertial(const size_t fi) {
			//std::cout << Source_Orientation;
			Vec3cd sf = FM[fi].secondary_inertial_frame(calculationtype, Source_Orientation);
			//std::cout << sf << std::endl;
			return sf;
		};

		void set_iptype(const IPType& _iptype) {
			//Todo
			//Earth.set_iptype(_iptype);
		};

	private:

	};
};
};