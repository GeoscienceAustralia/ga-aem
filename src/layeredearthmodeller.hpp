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
	using sPropogationMatrix = Eigen::Matrix2cd;

	struct sHankelTransform {
		Vec3cd FM;
		Vec3cd dC;
		Vec3cd dT;
		Vec3cd dX;
		Vec3cd dY;
		Vec3cd dZ;
		Vec3cd dH;
	};

	struct sAbscissaLayerNode {
		cdouble   U;
		cdouble   Exp2UT;
		sPropogationMatrix LayerMatrix;
		sPropogationMatrix LayerPreMatrix;
		sPropogationMatrix LayerPostMatrix;
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
		std::vector<sAbscissaLayerNode> Layer;
		sPropogationMatrix P_Full;
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

	struct sLayerNode {
		double Conductivity;
		double Thickness;
	};

	class LESingleFrequencyModeller {

	private:

		double X = 0.0, Y = 0.0, Z = 0.0, H = 0.0;
		double X2 = 0, Y2 = 0, X4 = 0, Y4 = 0, ZH = 0, ZH2 = 0;
		double r = 0, r2 = 0, r3 = 0, r5 = 0, r4 = 0, R = 0, R2 = 0, R5 = 0, R7 = 0;

		bool geometrychanged;
		bool earthchanged;

		std::vector<AbscissaNode> Abscissa;
		std::vector<double> Conductivity;
		std::vector<double> Thickness;


		//Hankle Stuff
		sHankelTransform HT;
		//sHankelTransform T0;
		//sHankelTransform T1;
		//sHankelTransform T2;

		//Vec3cd trapezoid_result;
		//Vec3cd integrand_result;
		//cdouble trapezoid_result[3];
		//cdouble integrand_result[3];

		double meanconductivity;
		double meanlog10conductivity;

		double LowerFractionalWidth;
		double UpperFractionalWidth;

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

		LESingleFrequencyModeller() { initialise(); };

		size_t na() const { return Abscissa.size(); }
		size_t nl() const { return Conductivity.size(); }
		const std::vector<double>& getConductivity() {
			return Conductivity;
		};

		const std::vector<double>& getThickness() {
			return Thickness;
		};

		void initialise() {
			LowerFractionalWidth = 4.44;
			UpperFractionalWidth = 1.84;
			//Ground EM 31 and 38
			//NumLayers()     = 0;  
			//NumAbscissa   = 1581;
			//LowerFractionalWidth  = 4.44;
			//UpperFractionalWidth  = 10.84;
			setnumabscissa(17);
		};

		void setupcomputations() {
			if (earthchanged || geometrychanged) {
				setfrequencyabscissalayers();
				earthchanged = false;
				geometrychanged = false;
			}
			else{
				int dummy = 0;
			}
		}

		void setnumabscissa(const size_t nabscissa)
		{
			Abscissa.resize(nabscissa);
		}

		void setconductivitythickness(const std::vector<double>& conductivity, const std::vector<double>& thickness)
		{
			const size_t nlayers = conductivity.size();
			Conductivity = conductivity;
			Thickness = thickness;
			for (size_t ai = 0; ai < na(); ai++) {
				Abscissa[ai].Layer.resize(nlayers);
			}
			setmeanconductivity();
			setmeanlog10conductivity();
			earthchanged = true;
		};

		void printearth()
		{
			for (size_t i = 0; i < nl() - 1; i++) {
				printf("Layer %02zu:\t%10lf mS/m\t%10lf m\n", i + 1, 1000.0 * Conductivity[i], Thickness[i]);
			}
			printf("Layer %02zu:\t%10lf mS/m\n", nl(), 1000.0 * Conductivity[nl() - 1]);
		}

		void setmeanlog10conductivity()
		{
			//Returns the thickness weighted mean (in linear space) conductivity (but calculated in 10g10 space)
			double sumc = 0.0;
			double sumt = 0.0;
			for (size_t i = 0; i < nl(); i++) {
				if (i < (nl() - 1)) {
					sumc += log10(Conductivity[i]) * Thickness[i];
					sumt += Thickness[i];
				}
				else {
					sumc += log10(Conductivity[i]) * sumt; //make basement layer as thick as sum of all overlying
					sumt += sumt;
				}
			}
			meanlog10conductivity = pow(10.0, sumc / sumt);
		}

		void setmeanconductivity()
		{
			//Returns the thickness weighted mean (in linear space) conductivity (but calculated in linear space)
			double sumc = 0.0;
			double sumt = 0.0;
			for (size_t i = 0; i < nl(); i++) {
				if (i < (nl() - 1)) {
					sumc += Conductivity[i] * Thickness[i];
					sumt += Thickness[i];
				}
				else {
					sumc += Conductivity[i] * sumt; //make basement layer as thick as sum of all overlying
					sumt += sumt;
				}
			}
			meanconductivity = sumc / sumt;
		}

		void setxyzh(const double& x, const double& y, const double& z, const double& h)
		{
			bool update_r = false;
			if (X != x) {
				X = x;
				X2 = X * X;
				X4 = X2 * X2;
				update_r = true;
			}

			if (y != Y) {
				Y = y;
				Y2 = Y * Y;
				Y4 = Y2 * Y2;
				update_r = true;
			}

			if (update_r) {
				r2 = X2 + Y2;
				r = std::sqrt(r2);
				r3 = r2 * r;
				r4 = r2 * r2;
				r5 = r4 * r;
			}

			bool update_zh = false;
			if (Z != z) {
				Z = z;
				update_zh = true;
			}

			if (H != h) {
				H = h;
				update_zh = true;
			}

			if (update_zh) {
				ZH = Z - H;
				ZH2 = ZH * ZH;
			}

			if (update_r || update_zh) {
				R2 = (X2 + Y2 + ZH2);
				R = std::sqrt(R2);
				R5 = R2 * R2 * R;
				R7 = R5 * R2;
				geometrychanged = true;
			}
			else {
				geometrychanged = false;
			}
		}

		void setfrequency(const double& frequency)
		{
			Frequency = frequency;
			Omega = TWOPI<double> * frequency;
			MuZeroOmega = MUZERO<double> * Omega;
			iMuZeroOmega = cdouble(0.0, MuZeroOmega);
		};

	private:
		void setfrequencyabscissalayers()
		{
			setintegrationnodes();
			for (size_t ai = 0; ai < na(); ai++) {
				for (size_t li = 0; li < nl(); li++) {
					double gamma2 = Conductivity[li] * MuZeroOmega;
					cdouble u = std::sqrt(cdouble(Abscissa[ai].Lambda2, gamma2));
					Abscissa[ai].Layer[li].U = u;
					if (li < nl() - 1)Abscissa[ai].Layer[li].Exp2UT = exp(-2.0 * u * Thickness[li]);
				}
				setlayermatrices(ai);
				setpmatrix(ai);
			}
		};

		double approximatehalfspace()
		{
			if (nl() == 1)return Conductivity[0];

			//peak lambda
			double   peaklambda = std::sqrt(MuZeroOmega * meanlog10conductivity / 4.0);
			cdouble  rz = rzero_recursive(peaklambda);
			cdouble  v = (1.0 + rz);

			v = iMuZeroOmega * v * v;
			cdouble a = (-4.0 * peaklambda * peaklambda * rz / v);
			return a.real();
		}

		cdouble rzero_recursive(const double lambda) const {
			//Wait's recursive formulation

			double lambda2 = lambda * lambda;
			double muzeroomega = MuZeroOmega;
			cdouble imuzeroomega(0.0, muzeroomega);

			double gamma2 = muzeroomega * Conductivity[nl() - 1];
			cdouble u = std::sqrt(cdouble(lambda2, gamma2));
			cdouble y = u / imuzeroomega;
			cdouble Nn, tanhv, v1, v2, v4, top, bot;

			for (int i = (int)nl() - 2; i >= 0; i--) {
				gamma2 = muzeroomega * Conductivity[i];
				u = std::sqrt(cdouble(lambda2, gamma2));
				Nn = u / imuzeroomega;

				v1 = u * Thickness[i];

				//Expand - unstable    	
				//tanh(v) = (1.0 - v4)/(1.0 + v4 + 2.0*v2);

				v2 = std::exp(-2.0 * v1);
				v4 = std::exp(-4.0 * v1);
				top = 1.0 - v4;
				bot = 1.0 + v4 + 2.0 * v2;
				tanhv = top / bot;

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

		void setlayermatrices(const size_t ai)
		{
			cdouble e, eh, e1, e2;

			AbscissaNode& A = Abscissa[ai];

			//M1  
			e = A.Layer[0].U / A.Lambda;
			eh = e / 2.0;
			e1 = 0.5 + eh;
			e2 = 0.5 - eh;

			sPropogationMatrix M;
			M(0, 0) = e1;
			M(0, 1) = e2;
			M(1, 0) = e2;
			M(1, 1) = e1;
			A.Layer[0].LayerMatrix = M;

			for (size_t li = 1; li < nl(); li++) {
				//assumes all pearmabilities are muzero
				e = A.Layer[li].U / A.Layer[li - 1].U;
				eh = e / 2.0;
				e1 = 0.5 + eh;
				e2 = 0.5 - eh;
				M(0, 0) = e1;
				M(0, 1) = e2;
				M(1, 0) = e2 * A.Layer[li - 1].Exp2UT;
				M(1, 1) = e1 * A.Layer[li - 1].Exp2UT;

				A.Layer[li].LayerMatrix = M;
			}

			//Set Prematrices - prematrix for layer 0 does not apply
			if (nl() > 1) A.Layer[1].LayerPreMatrix = A.Layer[0].LayerMatrix;
			for (size_t li = 2; li < nl(); li++) {
				A.Layer[li].LayerPreMatrix = A.Layer[li - 1].LayerPreMatrix * A.Layer[li - 1].LayerMatrix;
			}

			//Set Postmatrices - postmatrix for layer N-1 does not apply  
			if (nl() > 1) {
				A.Layer[nl() - 2].LayerPostMatrix = A.Layer[nl() - 1].LayerMatrix;
			}

			for (int li = (int)nl() - 3; li >= 0; li--) {
				A.Layer[li].LayerPostMatrix = A.Layer[(size_t)li + 1].LayerMatrix * A.Layer[(size_t)li + 1].LayerPostMatrix;
			}

		};

		void setpmatrix(const size_t ai)
		{
			AbscissaNode& A = Abscissa[ai];
			//Set Full matrix
			if (nl() == 1) {
				A.P_Full = A.Layer[0].LayerMatrix;
			}
			else {
				A.P_Full = Abscissa[ai].Layer[nl() - 1].LayerPreMatrix * Abscissa[ai].Layer[nl() - 1].LayerMatrix;
			}
			A.P21onP11 = A.P_Full(1, 0) / A.P_Full(0, 0);
		}

		sPropogationMatrix dMjdCj(const size_t ai, const size_t li)
		{
			sPropogationMatrix m;

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

		sPropogationMatrix dMjplus1dCj(const size_t ai, const size_t li)
		{
			cdouble duds = iMuZeroOmega / (2.0 * Abscissa[ai].Layer[li].U);
			cdouble y = Abscissa[ai].Layer[li + 1].U / Abscissa[ai].Layer[li].U;

			cdouble dydu = -y / Abscissa[ai].Layer[li].U;
			cdouble dyds = dydu * duds;

			cdouble v = Abscissa[ai].Layer[li].Exp2UT;
			cdouble dvdu = -2.0 * Thickness[li] * Abscissa[ai].Layer[li].Exp2UT;

			cdouble dvds = dvdu * duds;

			cdouble ydvds = y * dvds;
			cdouble vdyds = v * dyds;


			sPropogationMatrix M;
			M(0, 0) = 0.5 * dyds;
			M(0, 1) = -M(0, 0);
			M(1, 0) = 0.5 * (dvds - (ydvds + vdyds));
			M(1, 1) = 0.5 * (dvds + (ydvds + vdyds));

			return M;

		}

		sPropogationMatrix dPdCj(const size_t ai, const size_t li)
		{
			AbscissaNode& A = Abscissa[ai];
			//One layer case
			if (nl() == 1) return dMjdCj(ai, li);

			//Not last layer
			if (li < nl() - 1) {
				sPropogationMatrix c = dMjdCj(ai, li) * A.Layer[li + 1].LayerMatrix + A.Layer[li].LayerMatrix * dMjplus1dCj(ai, li);
				//First layer case
				if (li == 0) {
					if (nl() == 2) return c;
					return c * A.Layer[li + 1].LayerPostMatrix;
				}
				else if (li == nl() - 2) {
					return A.Layer[li].LayerPreMatrix * c;
				}
				else {
					return A.Layer[li].LayerPreMatrix * c * A.Layer[li + 1].LayerPostMatrix;
				}
			}
			//Last layer
			else {
				return A.Layer[li].LayerPreMatrix * dMjdCj(ai, li);
			}

		}

		cdouble dP21onP11dCj(const size_t ai, const size_t li)
		{
			sPropogationMatrix m = dPdCj(ai, li);
			return m(1, 0) / Abscissa[ai].P_Full(0, 0) - m(0, 0) * Abscissa[ai].P21onP11 / Abscissa[ai].P_Full(0, 0);
		}

		sPropogationMatrix dMjplus1dTj(const size_t ai, const size_t li)
		{
			cdouble y = Abscissa[ai].Layer[li + 1].U / Abscissa[ai].Layer[li].U;
			cdouble dvdt = -2.0 * Abscissa[ai].Layer[li].U * Abscissa[ai].Layer[li].Exp2UT;

			sPropogationMatrix m;
			m(0, 0) = 0.0;
			m(0, 1) = 0.0;
			m(1, 0) = 0.5 * (1.0 - y) * dvdt;
			m(1, 1) = 0.5 * (1.0 + y) * dvdt;
			return m;

		}

		sPropogationMatrix dPdTj(const size_t ai, const size_t li) {
			AbscissaNode& A = Abscissa[ai];

			//One layer case
			if (nl() == 1) return sPropogationMatrix::Zero();
			sPropogationMatrix b = A.Layer[li].LayerMatrix * dMjplus1dTj(ai, li);

			//First layer case
			if (li == 0) {
				if (nl() == 2) return b;
				return b * A.Layer[li + 1].LayerPostMatrix;
			}
			else if (li > 0 && li < nl() - 1) {
				sPropogationMatrix d = A.Layer[li].LayerPreMatrix * b;
				if (li == nl() - 2) {
					return d;
				}
				return d * A.Layer[li + 1].LayerPostMatrix;
			}
			//Last layer case
			else {
				glog.errormsg(_SRC_, "Zero thickness derivative for halfspace layer\n");
				return sPropogationMatrix::Zero();
			}
		}

		cdouble dP21onP11dTj(const size_t ai, const size_t li)
		{
			sPropogationMatrix m = dPdTj(ai, li);
			return m(1, 0) / Abscissa[ai].P_Full(0, 0) - m(0, 0) * Abscissa[ai].P21onP11 / Abscissa[ai].P_Full(0, 0);
		}

		void setintegrationnodes() {
			double peak_exp2 = 2.0 / (Z + H);
			double peak_exp3 = 3.0 / (Z + H);

			ApproximateHalfspace = approximatehalfspace();
			PeakLambda = std::sqrt(MuZeroOmega * ApproximateHalfspace / 4.0);

			double lp = std::log(std::min(PeakLambda, peak_exp2));
			double up = std::log(std::max(PeakLambda, peak_exp3));

			LowerBound = lp - LowerFractionalWidth;
			UpperBound = up + UpperFractionalWidth;

			AbscissaSpacing = (UpperBound - LowerBound) / (double)(na() - 1);

			double lambda;
			double loglambda = LowerBound;
			for (size_t ai = 0; ai < na(); ai++) {
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

		void dointegrals(const CalculationType& calculationtype) {
			Vec3cd integral = integrate_trapezoidal(calculationtype);

			switch (calculationtype.get_mode()){
			case CMode::FM:
				HT.FM = integral;
				break;
			case CMode::DC:
				HT.dC = integral;
				break;
			case CMode::DT:
				HT.dT = integral;
				break;
			case CMode::DZ:
				HT.dZ = integral;
				break;
			case CMode::DH:
				HT.dH = integral;
				break;
			case CMode::DX:
				HT.dX = integral;
				break;
			case CMode::DY:
				HT.dY = integral;
				break;
			default:
				glog.errormsg(_SRC_,"Unknown calculation type %s\n", calculationtype.string().c_str());
			}
		}

		Vec3cd integrate_trapezoidal(const CalculationType& calculationtype) {
			Vec3cd sum = Vec3cd::Zero();
			// First abscissa
			sum = 0.5 * compute_integrand(0, calculationtype);
			// Cenral abscissas
			for (size_t ai = 1; ai < na() - 1; ai++) {
				sum += compute_integrand(ai, calculationtype);
			}
			// Last Abscissa
			sum += 0.5 * compute_integrand(na() - 1, calculationtype);

			// Scale for abscissa spacing
			sum *= AbscissaSpacing;
			return sum;
		}

		Vec3cd compute_integrand(const size_t& ai, const CalculationType& calculationtype) {
			AbscissaNode& A = Abscissa[ai];

			const double& lambdar = A.Lambda_r;
			const double& j0 = A.j0Lambda_r;
			const double& j1 = A.j1Lambda_r;
			const double& e = exp(-(Z + H) * A.Lambda);
			const double& l2e = A.Lambda2 * e;
			const double& l3e = A.Lambda3 * e;
			const double& l4e = A.Lambda4 * e;

			double k0, k1, k2;
			k0 = k1 = k2 = std::numeric_limits<double>::max();
			cdouble earthkernel;

			const size_t& derivativelayer = calculationtype.get_layer();
			switch (calculationtype.get_mode()){
			case CMode::FM:
				earthkernel = -A.P21onP11;
				k0 = l3e * j0;
				k1 = l3e * j1;
				k2 = l2e * j1;
				break;
			case CMode::DX:
				earthkernel = -A.P21onP11;
				k0 = -l4e * j1 * X / r;
				k1 = l4e * (j0 - j1 / lambdar) * X / r;
				k2 = l3e * (j0 - j1 / lambdar) * X / r;
				break;
			case CMode::DY:
				earthkernel = -A.P21onP11;
				k0 = -l4e * j1 * Y / r;
				k1 = l4e * (j0 - j1 / lambdar) * Y / r;
				k2 = l3e * (j0 - j1 / lambdar) * Y / r;
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
				glog.errormsg(_SRC_,"Error: compute_integrand() unknown calculation type %c\n", calculationtype);
			}

			const double loopfactor = A.LoopFactor(ModellingLoopRadius);

			Vec3cd integrand;
			integrand[0] = earthkernel * (k0 * loopfactor);
			integrand[1] = earthkernel * (k1 * loopfactor);
			integrand[2] = earthkernel * (k2 * loopfactor);
			return integrand;
		}

		// Tensors
		Mat3d PTFM() const {
			Mat3d T;
			T(0, 0) = (3.0 * X2 - R2) / R5;
			T(0, 1) = 3.0 * X * Y / R5;
			T(0, 2) = 3.0 * X * ZH / R5;

			//Note error in Fitterman and Yin paper should no be minus sign at element 2,1
			T(1, 0) = T(0, 1);
			T(1, 1) = (3.0 * Y2 - R2) / R5;
			T(1, 2) = 3.0 * Y * ZH / R5;

			T(2, 0) = T(0, 2);
			T(2, 1) = T(1, 2);
			T(2, 2) = (3.0 * ZH2 - R2) / R5;
			return -ONEONFOURPI<double> * T;
		};

		Mat3d dPTdX() const {
			Mat3d T;
			T(0, 0) = -3.0 * X * (2.0 * X2 - 3.0 * Y2 - 3.0 * ZH2) / R7;
			T(0, 1) = -3.0 * Y * (4.0 * X2 - Y2 - ZH2) / R7;
			T(0, 2) = -3.0 * ZH * (4.0 * X2 - Y2 - ZH2) / R7;
			T(1, 0) = T(0, 1);
			T(1, 1) = 3.0 * X * (X2 - 4.0 * Y2 + ZH2) / R7;
			T(1, 2) = -15.0 * Y * ZH / R7 * X;
			T(2, 0) = T(0, 2);
			T(2, 1) = T(1, 2);
			T(2, 2) = 3.0 * X * (X2 + Y2 - 4.0 * ZH2) / R7;
			return -ONEONFOURPI<double> * T;
		};

		Mat3d dPTdY() const {
			Mat3d T;
			T(0, 0) = -3.0 * Y * (4.0 * X2 - Y2 - ZH2) / R7;
			T(0, 1) = 3.0 * X * (X2 - 4.0 * Y2 + ZH2) / R7;
			T(0, 2) = -15.0 * X * ZH / R7 * Y;
			T(1, 0) = T(0, 1);
			T(1, 1) = 3.0 * Y * (3.0 * X2 - 2.0 * Y2 + 3.0 * ZH2) / R7;
			T(1, 2) = 3.0 * ZH * (X2 - 4.0 * Y2 + ZH2) / R7;
			T(2, 0) = T(0, 2);
			T(2, 1) = T(1, 2);
			T(2, 2) = 3.0 * Y * (X2 + Y2 - 4.0 * ZH2) / R7;
			return -ONEONFOURPI<double> *T;
		};

		Mat3d dPTdZ() const {
			Mat3d T;
			T(0, 0) = -3.0 * ZH * (4.0 * X2 - Y2 - ZH2) / R7;
			T(0, 1) = -15.0 * X * Y / R7 * ZH;
			T(0, 2) = 3.0 * X * (X2 + Y2 - 4.0 * ZH2) / R7;
			T(1, 0) = T(0, 1);
			T(1, 1) = 3.0 * ZH * (X2 - 4.0 * Y2 + ZH2) / R7;
			T(1, 2) = 3.0 * Y * (X2 + Y2 - 4.0 * ZH2) / R7;
			T(2, 0) = T(0, 2);
			T(2, 1) = T(1, 2);
			T(2, 2) = 3.0 * ZH * (3.0 * X2 + 3.0 * Y2 - 2.0 * ZH2) / R7;
			return -ONEONFOURPI<double> *T;
		};

		Mat3d dPTdH() const {
			Mat3d m;
			//of course this is just minus d/dZ		
			m(0, 0) = 3.0 * ZH * (4.0 * X2 - Y2 - ZH2) / R7;
			m(0, 1) = 15.0 * X * Y / R7 * ZH;
			m(0, 2) = -3.0 * X * (X2 + Y2 - 4.0 * ZH2) / R7;
			m(1, 0) = m(0, 1);
			m(1, 1) = -3.0 * ZH * (X2 - 4.0 * Y2 + ZH2) / R7;
			m(1, 2) = -3.0 * Y * (X2 + Y2 - 4.0 * ZH2) / R7;
			m(2, 0) = m(0, 2);
			m(2, 1) = m(1, 2);
			m(2, 2) = -3.0 * ZH * (3.0 * X2 + 3.0 * Y2 - 2.0 * ZH2) / R7;
			return -ONEONFOURPI<double> *m;
		};

		Mat3cd STFM() const {
			const cdouble& T0 = HT.FM[0];
			const cdouble& T1 = HT.FM[1];
			const cdouble& T2 = HT.FM[2];
			
			Mat3cd m;
			m(0, 0) = ((X2 / r2 - Y2 / r2) * T2 / r - T0 * X2 / r2);
			m(0, 1) = (X * Y / r2) * (2.0 * T2 / r - T0);
			m(0, 2) = (-X / r) * T1;

			//Note error in Fitterman and Yin paper should not be minus sign at element 2,1
			m(1, 0) = m(0, 1);
			m(1, 1) = ((Y2 / r2 - X2 / r2) * T2 / r - T0 * Y2 / r2);
			m(1, 2) = (-Y / r) * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> *m;
		};

		Mat3cd dSTdC() const {
			const cdouble& T0 = HT.dC[0];
			const cdouble& T1 = HT.dC[1];
			const cdouble& T2 = HT.dC[2];

			Mat3cd m;
			m(0, 0) = ((X2 / r2 - Y2 / r2) * T2 / r - T0 * X2 / r2);
			m(0, 1) = (X * Y / r2) * (2.0 * T2 / r - T0);
			m(0, 2) = (-X / r) * T1;

			m(1, 0) = m(0, 1);
			m(1, 1) = ((Y2 / r2 - X2 / r2) * T2 / r - T0 * Y2 / r2);
			m(1, 2) = (-Y / r) * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> *m;
		};

		Mat3cd dSTdT() const {
			const cdouble& T0 = HT.dT[0];
			const cdouble& T1 = HT.dT[1];
			const cdouble& T2 = HT.dT[2];
			
			Mat3cd m;
			m(0, 0) = ((X2 / r2 - Y2 / r2) * T2 / r - T0 * X2 / r2);
			m(0, 1) = (X * Y / r2) * (2.0 * T2 / r - T0);
			m(0, 2) = (-X / r) * T1;

			m(1, 0) = m(0, 1);
			m(1, 1) = ((Y2 / r2 - X2 / r2) * T2 / r - T0 * Y2 / r2);
			m(1, 2) = (-Y / r) * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> *m;
		};

		Mat3cd dSTdX() const {
			const cdouble& T0 = HT.dX[0];
			const cdouble& T1 = HT.dX[1];
			const cdouble& T2 = HT.dX[2];
			const cdouble& T0FM = HT.FM[0];
			const cdouble& T1FM = HT.FM[1];
			const cdouble& T2FM = HT.FM[2];

			Mat3cd m;
			m(0, 0) = (X4 * T2 - T0 * X4 * r - T2FM * X2 * X - T0 * X2 * r * Y2 - 2.0 * T0FM * X * r * Y2 + 5.0 * X * Y2 * T2FM - Y4 * T2) / r5;
			m(0, 1) = 2.0 * Y / r3 * T2FM - Y / r2 * T0FM - 6.0 * X2 * Y / r5 * T2FM + 2.0 * X2 * Y / r4 * T0FM + 2.0 * X * Y / r3 * T2 - X * Y / r2 * T0;
			m(0, 2) = -(T1FM * Y2 + X2 * X * T1 + X * T1 * Y2) / r3;

			m(1, 0) = m(0, 1);
			m(1, 1) = -(T2 * X4 - T2FM * X2 * X + T0 * Y2 * r * X2 - 2.0 * T0FM * Y2 * X * r + 5.0 * X * Y2 * T2FM + T0 * Y4 * r - Y4 * T2) / r5;
			m(1, 2) = Y / r3 * T1FM * X - Y / r * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> * m;
		};

		Mat3cd dSTdY() const {
			const cdouble& T0 = HT.dY[0];
			const cdouble& T1 = HT.dY[1];
			const cdouble& T2 = HT.dY[2];
			const cdouble& T0FM = HT.FM[0];
			const cdouble& T1FM = HT.FM[1];
			const cdouble& T2FM = HT.FM[2];

			Mat3cd m;
			m(0, 0) = (X4 * T2 - T0 * X4 * r + 2.0 * T0FM * X2 * Y * r - T0 * X2 * r * Y2 - 5.0 * X2 * Y * T2FM + Y2 * Y * T2FM - Y4 * T2) / r5;
			m(0, 1) = 2.0 * X / r3 * T2FM - X / r2 * T0FM - 6.0 * X * Y2 / r5 * T2FM + 2.0 * X * Y2 / r4 * T0FM + 2.0 * X * Y / r3 * T2 - X * Y / r2 * T0;
			m(0, 2) = X / r3 * T1FM * Y - X / r * T1;

			m(1, 0) = m(0, 1);
			m(1, 1) = -(T2 * X4 + T0 * Y2 * r * X2 + 2.0 * T0 * Y * r * X2 - 5.0 * X2 * Y * T2FM + T0 * Y4 * r + Y2 * Y * T2 - Y4 * T2) / r5;
			m(1, 2) = -(T1 * X2 + Y * T1 * X2 + Y2 * Y * T1) / r3;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> * m;
		};

		Mat3cd dSTdZ() const {
			const cdouble& T0 = HT.dZ[0];
			const cdouble& T1 = HT.dZ[1];
			const cdouble& T2 = HT.dZ[2];

			Mat3cd m;
			m(0, 0) = ((X2 / r2 - Y2 / r2) * T2 / r - T0 * X2 / r2);
			m(0, 1) = X * Y / r2 * (2.0 * T2 / r - T0);
			m(0, 2) = -X / r * T1;

			m(1, 0) = m(0, 1);
			m(1, 1) = ((Y2 / r2 - X2 / r2) * T2 / r - T0 * Y2 / r2);
			m(1, 2) = -Y / r * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> *m;
		};

		Mat3cd dSTdH() const {
			const cdouble& T0 = HT.dH[0];
			const cdouble& T1 = HT.dH[1];
			const cdouble& T2 = HT.dH[2];

			Mat3cd m;
			m(0, 0) = ((X2 / r2 - Y2 / r2) * T2 / r - T0 * X2 / r2);
			m(0, 1) = X * Y / r2 * (2.0 * T2 / r - T0);
			m(0, 2) = -X / r * T1;

			m(1, 0) = m(0, 1);
			m(1, 1) = ((Y2 / r2 - X2 / r2) * T2 / r - T0 * Y2 / r2);
			m(1, 2) = -Y / r * T1;

			m(2, 0) = -m(0, 2);
			m(2, 1) = -m(1, 2);
			m(2, 2) = -T0;
			return -ONEONFOURPI<double> * m;
		};

		Mat3d PrimaryTensor(const CalculationType& calculationtype) {
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
			dointegrals(calculationtype);

			switch (calculationtype.get_mode()) {
			case CMode::FM: return STFM();
			case CMode::DX: return dSTdX();
			case CMode::DY: return dSTdY();
			case CMode::DZ: return dSTdZ();
			case CMode::DH: return dSTdH();
			case CMode::DC: return dSTdC();
			case CMode::DT: return dSTdT();
			default:
				glog.errormsg(_SRC_, "Unknown calculation type %c\n", calculationtype);
			}
		};

	public:

		Vec3d primary_inertial_frame(const CalculationType& calculationtype, const Vec3d& txdir) {
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

	public:
		std::vector<LESingleFrequencyModeller> FM;

		LEModeller() {};

		size_t nFrequencies() const { return FM.size(); };
		size_t nLayers() const { 
			// Todo
			return FM[0].nl(); 
		};

		void set_numabscissa(const size_t& na) {
			//Todo
		};

		void set_modellingloopradius(const double& radius) {
			ModellingLoopRadius = radius;
		};

		void initialise_frequencies(
			const std::vector<double> discrete_frequencies
		){
			const size_t nf = discrete_frequencies.size();
			FM.resize(nf);
			for (size_t i = 0; i < nf; i++) {
				FM[i].initialise();
				FM[i].setfrequency(discrete_frequencies[i]);
				FM[i].ModellingLoopRadius = ModellingLoopRadius;
			}
		};

		void set_earth(const AEM::Earth1D& E) {
			const size_t nf = nFrequencies();
			for (size_t i = 0; i < nf; i++) {
				FM[i].setconductivitythickness(E.conductivity, E.thickness);
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