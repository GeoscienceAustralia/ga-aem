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
#include <Eigen/Dense>

#include "general_constants.hpp"
#include "general_utils.hpp"
#include "earth1d.hpp"
#include "calculation_type.hpp"
#include "eigen_utils.hpp"

//Formulation mainly from the book 
//Geo-Electromagnetism, Wait James, R. Academic Press 1982
namespace AEM {
namespace LEM1 {
	using cdouble = std::complex<double>;
	using cvector = std::vector<std::complex<double>>;
	using CalculationType = CT::CalculationType;
	using CMode = CT::CalculationType::Mode;
	enum class RZeroMethod { PROPOGATIONMATRIX, RECURSIVE };
	constexpr double DefaultLowerFractionalWidth = 4.44;
	constexpr double DefaultUpperFractionalWidth = 1.84;

	inline cdouble ip_colecole_conductivity(const double& conductivity, const double& chargeability, const double& timeconstant, const double& frequencydependence, const double& omega) {
		//c = c0 - c0*(N / (1 + (1 - N)*(j*omega*T) ^ K));
		if (chargeability == 0.0) {
			return cdouble(conductivity, 0.0);
		}
		else {
			cdouble c = conductivity - conductivity * (chargeability / (1.0 + (1.0 - chargeability) * (std::pow(cdouble(0.0, omega * timeconstant), frequencydependence))));
			return  c;
		}
	};

	inline cdouble ip_pelton_conductivity(const double& conductivity, const double& chargeability, const double& timeconstant, const double& frequencydependence, const double& omega) {
		//p = p0[1 - m*(1 - (1 - 1/(1 + (j*omega*T) ^ K));
		if (chargeability == 0.0) {
			return cdouble(conductivity, 0.0);
		}
		else {
			double  rho0 = 1.0 / conductivity;
			cdouble rho = rho0 * (1.0 - chargeability * (1.0 - (1.0 / (1.0 + std::pow(cdouble(0.0, omega * timeconstant), frequencydependence)))));
			return  1.0 / rho;
		}
	};

	inline cdouble ip_complex_conductivity(const IPType& iptype, const double& conductivity, const double& chargeability, const double& timeconstant, const double& frequencydependence, const double& omega) {
		cdouble complex_conductivity;
		if (iptype == IPType::NONE) complex_conductivity = conductivity;
		else if (iptype == IPType::COLECOLE) complex_conductivity = ip_colecole_conductivity(conductivity, chargeability, timeconstant, frequencydependence, omega);
		else complex_conductivity = ip_pelton_conductivity(conductivity, chargeability, timeconstant, frequencydependence, omega);
		return complex_conductivity;
	};

	struct HankelTransform {
		cdouble FM = 0.0;
		cdouble dC = 0.0;
		cdouble dT = 0.0;
		cdouble dZ = 0.0;
		cdouble dH = 0.0;
		cdouble dR = 0.0;
	};

	struct HankelTransforms {
		HankelTransform I0;
		HankelTransform I1;
		HankelTransform I2;
	};

	struct RealField {
		double x = 0.0;
		double y = 0.0;
		double z = 0.0;
	};

	struct ComplexField {
		cdouble x = 0.0;
		cdouble y = 0.0;
		cdouble z = 0.0;
	};

	struct TotalField {
		RealField    p = {};//primary filed component
		ComplexField s = {};//secondary field component  
	};

	struct ResponseField {
		TotalField   v = {};//due to vertical dipole
		TotalField   h = {};//due to y directed horizontal dipole
		TotalField   t = {};//due to total dipole
	};

	class PropogationMatrix {

	public:
		cdouble e11 = 0.0;
		cdouble e12 = 0.0;
		cdouble e21 = 0.0;
		cdouble e22 = 0.0;

		PropogationMatrix& operator+=(const PropogationMatrix& rhs)
		{
			e11 += rhs.e11;
			e12 += rhs.e12;
			e21 += rhs.e21;
			e22 += rhs.e22;
			return *this;
		}

		PropogationMatrix operator+(const PropogationMatrix& rhs)
		{
			PropogationMatrix m = *this;
			m.e11 += rhs.e11;
			m.e12 += rhs.e12;
			m.e21 += rhs.e21;
			m.e22 += rhs.e22;
			return m;
		}

		PropogationMatrix operator*(const PropogationMatrix& b)
		{
			PropogationMatrix c;
			c.e11 = e11 * b.e11 + e12 * b.e21;
			c.e12 = e11 * b.e12 + e12 * b.e22;
			c.e21 = e21 * b.e11 + e22 * b.e21;
			c.e22 = e21 * b.e12 + e22 * b.e22;
			return c;
		}

	};

	struct AbscissaLayerNode {
		cdouble   U = 0.0;
		cdouble   Exp2UT = 0.0;
		PropogationMatrix LayerMatrix;
		PropogationMatrix LayerPreMatrix;
		PropogationMatrix LayerPostMatrix;
	};

	struct AbscissaNode {
		double Lambda = 0.0;
		double Lambda2 = 0.0;
		double Lambda3 = 0.0;
		double Lambda4 = 0.0;
		double LambdaR = 0.0;
		double j0LambdaR = 0.0;
		double j1LambdaR = 0.0;
		std::vector<AbscissaLayerNode> Layer;
		PropogationMatrix P_Full;
		cdouble P21onP11 = 0.0;
	};

	struct FrequencyNode {
		double Frequency = 0.0;
		double Omega = 0.0;
		double MuZeroOmega = 0.0;
		cdouble iMuZeroOmega = 0.0;
		double ApproximateHalfspace = 0.0;
		double PeakLambda = 0.0;
		double LowerBound = 0.0;
		double UpperBound = 0.0;
		double AbscissaSpacing = 0.0;
		std::vector<AbscissaNode> Abscissa;
	};

	class LEModeller {

	private:
		//Geometry Stuff
		Vec3d Source_Orientation;
		double Xunrotated, Yunrotated, X, Y, Z, H;
		double R, R2, R3, R4, R5;
		double XonR, YonR;
		double BigR, BigR2, BigR3, BigR5, BigR7;
		double xyrotation, cosxyrotation, sinxyrotation;

		//Hankle Stuff
		inline static constexpr size_t NumIntegrands = 3;
		Vec3cd trapezoid_result;
		Vec3cd integrand_result;
		double mean_conductivity;
		double mean_log10conductivity;

		size_t NumAbscissa = 17;
		double LowerFractionalWidth = DefaultLowerFractionalWidth;
		double UpperFractionalWidth = DefaultUpperFractionalWidth;
		size_t number_integrand_calls;
		std::vector<HankelTransforms> Hankel;

		double ModellingLoopRadius = 0.0; //dipole by default	
		RZeroMethod rzerotype = RZeroMethod::PROPOGATIONMATRIX;
		CalculationType calculationtype;

		Earth1D Earth;
		std::vector<FrequencyNode> Frequencies;
		ResponseField Fields;

	public:
		size_t nFrequencies() const { return Frequencies.size(); };
		size_t nLayers() const { return Earth.nlayers(); };

		LEModeller() {};

		void initialise(const std::vector<double>& discrete_frequencies, const size_t& numabscissa, const double& modelling_loop_radius) {
			//const size_t nf = discrete_frequencies.size();
			//FM.resize(nf);
			//for (size_t fi = 0; fi < nf; fi++) {
			//	FM[fi].initialise(discrete_frequencies[fi], numabscissa, modelling_loop_radius);
			//}
			initialise_frequencies(discrete_frequencies);
			set_modellingloopradius(modelling_loop_radius);
			set_numabscissa(numabscissa);


		}


		void set_modellingloopradius(const double& radius) {
			ModellingLoopRadius = radius;
		};

		void set_numabscissa(const size_t& numabscissa) {
			NumAbscissa = numabscissa;
		};

		void set_earth(const Earth1D& E) {
			Earth = E;
			//set_earth_properties(E.conductivity, E.thickness, E.chargeability, E.timeconstant, E.frequencydependence);
			mean_conductivity = Earth.mean_weighted_conductivity();
			mean_log10conductivity = Earth.mean_weighted_conductivity_log10_calculation();
		};

		void set_iptype(const IPType& _iptype) {
			Earth.set_iptype(_iptype);
		};

		void set_calculationtype(const CalculationType& _calculationtype) {
			calculationtype = _calculationtype;
		};

		const CalculationType::Mode& cmode() const {
			return calculationtype.get_mode();
		};

		const size_t& derivative_layer() const {
			return calculationtype.get_layer();
		};
		
		void initialise_frequencies(const std::vector<double>& frequencies) {
			const size_t nf = frequencies.size();
			Frequencies.resize(nf);
			Hankel.resize(nf);
			for (size_t fi = 0; fi < nf; fi++) {
				double omega = TWOPI<double> * frequencies[fi];
				double muzeroomega = MUZERO<double> *omega;
				Frequencies[fi].Frequency = frequencies[fi];
				Frequencies[fi].Omega = omega;
				Frequencies[fi].MuZeroOmega = muzeroomega;
				Frequencies[fi].iMuZeroOmega = cdouble(0.0, muzeroomega);
			}
		};

		void set_geometry(const Vec3d& source_orientation, double h, double x, double y, double z) {
			Xunrotated = x;
			Yunrotated = y;
			Source_Orientation = source_orientation;
			setxyrotation();
			xyrotate(Xunrotated, Yunrotated, &X, &Y);
			R = sqrt(X * X + Y * Y);
			R2 = R * R;
			R3 = R2 * R;
			R4 = R3 * R;
			R5 = R4 * R;
			Z = z;
			H = h;

			BigR = sqrt((Z - H) * (Z - H) + R * R);
			BigR2 = BigR * BigR;
			BigR3 = BigR * BigR2;
			BigR5 = BigR3 * BigR2;
			BigR7 = BigR5 * BigR2;

			if (R == 0.0) {
				XonR = 0.0;
				YonR = 0.0;
			}
			else {
				XonR = X / R;
				YonR = Y / R;
			}
		}

		void setup_computations() {
			for (size_t fi = 0; fi < nFrequencies(); fi++) {
				initialise_frequency_for_computation(fi);
			}
		}

		Vec3d primaryfield_inertial() {
			set_primaryfields();
			return Vec3d(Fields.t.p.x, Fields.t.p.y, Fields.t.p.z);
		};

		Vec3cd secondaryfield_inertial(const size_t fi) {
			set_secondaryfields(fi);
			return Vec3cd(Fields.t.s.x, Fields.t.s.y, Fields.t.s.z);
		};

	private:
		void set_primaryfields() {
			sethorizontaldipoleprimaryfields();
			setverticaldipoleprimaryfields();
			Fields.t.p.x = Fields.v.p.x + Fields.h.p.x;
			Fields.t.p.y = Fields.v.p.y + Fields.h.p.y;
			Fields.t.p.z = Fields.v.p.z + Fields.h.p.z;
		};

		void set_secondaryfields(const size_t& fi) {
			dointegrals(fi);
			sethorizontaldipolesecondaryfields(fi);
			setverticaldipolesecondaryfields(fi);
			Fields.t.s.x = Fields.v.s.x + Fields.h.s.x;
			Fields.t.s.y = Fields.v.s.y + Fields.h.s.y;
			Fields.t.s.z = Fields.v.s.z + Fields.h.s.z;
		};
		
		void setxyrotation() {
			//xyrotation is the anticlockwise angle (in degrees) 
			//that the horizontal coordinate system has to be rotated
			//so the horizontal dipole is all Y directed. 
			if (Source_Orientation.x() == 0.0 && Source_Orientation.y() == 0.0) {
				xyrotation = 0.0;
			}
			else {
				xyrotation = atan2(Source_Orientation.y(), Source_Orientation.x());
				xyrotation = xyrotation * R2D<double> -90.0;
			}
			cosxyrotation = cos(xyrotation * D2R<double>);
			sinxyrotation = sin(xyrotation * D2R<double>);
		}

		void xyrotate(const double& xin, const double& yin, double* xout, double* yout) const {
			*xout = xin * cosxyrotation + yin * sinxyrotation;
			*yout = -xin * sinxyrotation + yin * cosxyrotation;
		}

		void unxyrotate(const double& xin, const double& yin, double* xout, double* yout) const {
			*xout = xin * cosxyrotation - yin * sinxyrotation;
			*yout = xin * sinxyrotation + yin * cosxyrotation;
		}

		void unxyrotateandscale(RealField& f, double scalefactor) {
			double x, y;
			unxyrotate(f.x, f.y, &x, &y);
			f.x = x * scalefactor;
			f.y = y * scalefactor;
			f.z = f.z * scalefactor;
		}

		void unxyrotateandscale(ComplexField& f, double scalefactor) {
			RealField r;
			r.x = f.x.real();
			r.y = f.y.real();
			r.z = f.z.real();
			unxyrotateandscale(r, scalefactor);

			RealField i;
			i.x = f.x.imag();
			i.y = f.y.imag();
			i.z = f.z.imag();
			unxyrotateandscale(i, scalefactor);

			f.x = cdouble(r.x, i.x);
			f.y = cdouble(r.y, i.y);
			f.z = cdouble(r.z, i.z);

		}

		void setR(double r){
			R = r;
			R2 = R * R;
			R3 = R2 * R;
			R4 = R3 * R;
			R5 = R4 * R;
		}

		void initialise_frequency_for_computation(const size_t& fi)	{
			initialise_integration_nodes(fi);
			for (size_t ai = 0; ai < NumAbscissa; ai++) {
				initialise_abscissa(fi, ai);
			}
		};

		void initialise_abscissa(const size_t& fi, const size_t& ai) {
			const size_t nl = nLayers();
			Frequencies[fi].Abscissa[ai].Layer.resize(nl);
			for (size_t li = 0; li < nl; li++) {
				cdouble c = Earth.conductivity[li];
				if (Earth.chargeability.size() > 0) {
					c = ip_complex_conductivity(Earth.get_iptype(), Earth.conductivity[li], Earth.chargeability[li], Earth.timeconstant[li], Earth.frequencydependence[li], Frequencies[fi].Omega);
				};
				const cdouble gamma2 = c * cdouble(0.0, Frequencies[fi].MuZeroOmega);
				const cdouble u = std::sqrt(Frequencies[fi].Abscissa[ai].Lambda2 + gamma2);
				Frequencies[fi].Abscissa[ai].Layer[li].U = u;

				if (li < nl-1) Frequencies[fi].Abscissa[ai].Layer[li].Exp2UT = exp(-2.0 * u * Earth.thickness[li]);
			}
			init_layer_matrices(fi, ai);
			init_pmatrix(fi, ai);
		};

		double approximatehalfspace(const size_t& fi) {
			//approximate halfspace for the frequency at index fi
			if (nLayers() == 1) return Earth.conductivity[0];

			//peak lambda
			double  peaklambda = sqrt(Frequencies[fi].MuZeroOmega * mean_log10conductivity / 4.0);
			cdouble rz = rzero_recursive(fi, peaklambda);

			cdouble  v = (1.0 + rz);

			v = Frequencies[fi].iMuZeroOmega * v * v;
			return (-4.0 * peaklambda * peaklambda * rz / v).real();
		}

		inline cdouble rzero_recursive(const size_t& fi, const double& lambda) const {
			//Wait's recursive formulation
			const size_t nl = Earth.nlayers();
			const double lambda2 = lambda * lambda;
			const double muzeroomega = Frequencies[fi].MuZeroOmega;
			const cdouble imuzeroomega(0.0, muzeroomega);

			double gamma2 = muzeroomega * Earth.conductivity[nl - 1];
			cdouble u = std::sqrt(cdouble(lambda2, gamma2));
			cdouble y = u / imuzeroomega;

			int i = nl - 2;
			while (i >= 0) {
				gamma2 = muzeroomega * Earth.conductivity[i];
				u = std::sqrt(cdouble(lambda2, gamma2));
				const cdouble Nn = u / imuzeroomega;
				const cdouble v = u * Earth.thickness[i];

				//Expand - unstable    	
				//tanh(v) = (1.0 - v4)/(1.0 + v4 + 2.0*v2);
				const cdouble v2 = std::exp(-2.0 * v);
				const cdouble v4 = v2 * v2;
				const cdouble tanhv = (1.0 - v4) / (1.0 + v4 + 2.0 * v2);
				y = Nn * (y + Nn * tanhv) / (Nn + y * tanhv);
				i--;
			}
			const cdouble N0(0.0, -lambda / (muzeroomega));  // minus because dividing by i.muzero.omega
			return (N0 - y) / (N0 + y);
		};

		inline cdouble rzero_propogationmatrix(const size_t& fi, const size_t& ai)	{
			//Oldenberg's propogation matrix formulation
			init_layer_matrices(fi, ai);
			init_pmatrix(fi, ai);
			return Frequencies[fi].Abscissa[ai].P21onP11;
		};

		inline void init_layer_matrices(const size_t& fi, const size_t& ai)	{
			cdouble e, eh, e1, e2;
			AbscissaNode& A = Frequencies[fi].Abscissa[ai];

			//M1  
			e = A.Layer[0].U / A.Lambda;
			eh = e / 2.0;
			e1 = 0.5 + eh;
			e2 = 0.5 - eh;

			A.Layer[0].LayerMatrix.e11 = e1;
			A.Layer[0].LayerMatrix.e12 = e2;
			A.Layer[0].LayerMatrix.e21 = e2;
			A.Layer[0].LayerMatrix.e22 = e1;

			for (size_t li = 1; li < nLayers(); li++) {
				//assumes all pearmabilities are muzero
				e = A.Layer[li].U / A.Layer[li - 1].U;
				eh = e / 2.0;
				e1 = 0.5 + eh;
				e2 = 0.5 - eh;
				A.Layer[li].LayerMatrix.e11 = e1;
				A.Layer[li].LayerMatrix.e12 = e2;
				A.Layer[li].LayerMatrix.e21 = e2 * A.Layer[li - 1].Exp2UT;
				A.Layer[li].LayerMatrix.e22 = e1 * A.Layer[li - 1].Exp2UT;
			}

			//Set Prematrices - prematrix for first layer does not apply
			if (nLayers() > 1) {
				A.Layer[1].LayerPreMatrix = A.Layer[0].LayerMatrix;
				for (size_t li = 2; li < nLayers(); li++) {
					A.Layer[li].LayerPreMatrix = A.Layer[li - 1].LayerPreMatrix * A.Layer[li - 1].LayerMatrix;
				}
			}

			//Set Postmatrices - postmatrix for last layer does not apply  
			if (nLayers() > 1) {
				A.Layer[nLayers() - 2].LayerPostMatrix = A.Layer[nLayers() - 1].LayerMatrix;
				if (nLayers() > 2) {
					for (size_t li = nLayers() - 2; li-- > 0;) {
						A.Layer[li].LayerPostMatrix = A.Layer[li + 1].LayerMatrix * A.Layer[li + 1].LayerPostMatrix;
					}
				}
			}
		};

		inline void init_pmatrix(const size_t& fi, const size_t& ai)
		{
			AbscissaNode& A = Frequencies[fi].Abscissa[ai];
			//Set Full matrix
			if (nLayers() == 1) {
				A.P_Full = A.Layer[0].LayerMatrix;
			}
			else {
				A.P_Full = Frequencies[fi].Abscissa[ai].Layer[nLayers() - 1].LayerPreMatrix * Frequencies[fi].Abscissa[ai].Layer[nLayers() - 1].LayerMatrix;
			}
			A.P21onP11 = A.P_Full.e21 / A.P_Full.e11;
		}

		inline PropogationMatrix dMjdCj(const size_t& fi, const size_t& ai, const size_t& li)
		{
			PropogationMatrix m;

			if (li == 0) {
				cdouble a = Frequencies[fi].iMuZeroOmega / (4.0 * Frequencies[fi].Abscissa[ai].Lambda * Frequencies[fi].Abscissa[ai].Layer[li].U);
				m.e11 = a;
				m.e12 = -a;
				m.e21 = -a;
				m.e22 = a;
				return m;
			}
			else {
				cdouble a = Frequencies[fi].iMuZeroOmega / (4.0 * Frequencies[fi].Abscissa[ai].Layer[li - 1].U * Frequencies[fi].Abscissa[ai].Layer[li].U);
				cdouble ae = a * Frequencies[fi].Abscissa[ai].Layer[li - 1].Exp2UT;
				m.e11 = a;
				m.e12 = -a;
				m.e21 = -ae;
				m.e22 = ae;
				return m;
			}

		}

		inline PropogationMatrix dMjplus1dCj(const size_t& fi, const size_t& ai, const size_t& li) {
			cdouble duds = Frequencies[fi].iMuZeroOmega / (2.0 * Frequencies[fi].Abscissa[ai].Layer[li].U);
			cdouble y = Frequencies[fi].Abscissa[ai].Layer[li + 1].U / Frequencies[fi].Abscissa[ai].Layer[li].U;

			cdouble dydu = -y / Frequencies[fi].Abscissa[ai].Layer[li].U;
			cdouble dyds = dydu * duds;

			cdouble v = Frequencies[fi].Abscissa[ai].Layer[li].Exp2UT;
			cdouble dvdu = -2.0 * Earth.thickness[li] * Frequencies[fi].Abscissa[ai].Layer[li].Exp2UT;

			cdouble dvds = dvdu * duds;

			cdouble ydvds = y * dvds;
			cdouble vdyds = v * dyds;


			PropogationMatrix M;
			M.e11 = 0.5 * dyds;
			M.e12 = -M.e11;
			M.e21 = 0.5 * (dvds - (ydvds + vdyds));
			M.e22 = 0.5 * (dvds + (ydvds + vdyds));

			return M;

		}

		inline PropogationMatrix dPdCj(const size_t& fi, const size_t& ai, const size_t& li)
		{
			AbscissaNode& A = Frequencies[fi].Abscissa[ai];

			PropogationMatrix tmp;
			//One layer case
			if (nLayers() == 1)return tmp = dMjdCj(fi, ai, li);

			//Not last layer
			if (li < nLayers() - 1) {
				PropogationMatrix M =
					(dMjdCj(fi, ai, li) * A.Layer[li + 1].LayerMatrix)
					+ (A.Layer[li].LayerMatrix * dMjplus1dCj(fi, ai, li));

				//First layer case
				if (li == 0) {
					if (nLayers() == 2) return M;
					return M * A.Layer[li + 1].LayerPostMatrix;
				}
				else if (li == nLayers() - 2) {
					return A.Layer[li].LayerPreMatrix * M;
				}
				else {
					return A.Layer[li].LayerPreMatrix * M * A.Layer[li + 1].LayerPostMatrix;
				}
			}
			//Last layer
			else {
				return A.Layer[li].LayerPreMatrix * dMjdCj(fi, ai, li);
			}

		}

		inline cdouble dP21onP11dCj(const size_t& fi, const size_t& ai, const size_t& li)
		{
			PropogationMatrix m = dPdCj(fi, ai, li);
			return m.e21 / Frequencies[fi].Abscissa[ai].P_Full.e11 - m.e11 * Frequencies[fi].Abscissa[ai].P21onP11 / Frequencies[fi].Abscissa[ai].P_Full.e11;
		}

		inline PropogationMatrix dMjplus1dTj(const size_t& fi, const size_t& ai, const size_t& li)
		{
			cdouble y = Frequencies[fi].Abscissa[ai].Layer[li + 1].U / Frequencies[fi].Abscissa[ai].Layer[li].U;
			cdouble dvdt = -2.0 * Frequencies[fi].Abscissa[ai].Layer[li].U * Frequencies[fi].Abscissa[ai].Layer[li].Exp2UT;

			PropogationMatrix m;
			m.e11 = 0.0;
			m.e12 = 0.0;
			m.e21 = 0.5 * (1.0 - y) * dvdt;
			m.e22 = 0.5 * (1.0 + y) * dvdt;
			return m;

		}

		inline PropogationMatrix dPdTj(const size_t& fi, const size_t& ai, const size_t& li)
		{
			AbscissaNode& A = Frequencies[fi].Abscissa[ai];

			PropogationMatrix tmp;
			tmp.e11 = tmp.e12 = tmp.e21 = tmp.e22 = 0;

			//One layer case
			if (nLayers() == 1) return tmp;

			PropogationMatrix M = A.Layer[li].LayerMatrix * dMjplus1dTj(fi, ai, li);

			//First layer case
			if (li == 0) {
				if (nLayers() == 2)return M;
				return M * A.Layer[li + 1].LayerPostMatrix;
			}
			else if (li > 0 && li < nLayers() - 1) {
				if (li == nLayers() - 2) return A.Layer[li].LayerPreMatrix * M;
				return A.Layer[li].LayerPreMatrix * M * A.Layer[li + 1].LayerPostMatrix;
			}
			//Last layer case
			else {
				glog.warningmsg(_SRC_, "Zero thickness derivative for halfspace layer\n");
				tmp.e11 = tmp.e12 = tmp.e21 = tmp.e22 = 0;
				return tmp;
			}
		}

		inline cdouble dP21onP11dTj(const size_t& fi, const size_t& ai, const size_t& li)
		{
			PropogationMatrix m = dPdTj(fi, ai, li);
			return m.e21 / Frequencies[fi].Abscissa[ai].P_Full.e11 - m.e11 * Frequencies[fi].Abscissa[ai].P21onP11 / Frequencies[fi].Abscissa[ai].P_Full.e11;
		}

		void initialise_integration_nodes(const size_t& fi)
		{
			FrequencyNode& F = Frequencies[fi];;
			double peak_exp2 = 2.0 / (Z + H);
			double peak_exp3 = 3.0 / (Z + H);

			F.ApproximateHalfspace = approximatehalfspace(fi);
			F.PeakLambda = sqrt(F.MuZeroOmega * F.ApproximateHalfspace / 4.0);

			double lp = log(std::min(F.PeakLambda, peak_exp2));
			double up = log(std::max(F.PeakLambda, peak_exp3));

			F.LowerBound = lp - LowerFractionalWidth;
			F.UpperBound = up + UpperFractionalWidth;

			F.AbscissaSpacing = (F.UpperBound - F.LowerBound) / (double)(NumAbscissa - 1);
			F.Abscissa.resize(NumAbscissa);

			double loglambda = F.LowerBound;
			for (size_t ai = 0; ai < NumAbscissa; ai++) {
				AbscissaNode& A = F.Abscissa[ai];
				double lambda = exp(loglambda);
				A.Lambda = lambda;
				A.Lambda2 = A.Lambda * lambda;
				A.Lambda3 = A.Lambda2 * lambda;
				A.Lambda4 = A.Lambda3 * lambda;
				A.LambdaR = A.Lambda * R;
				A.j0LambdaR = std::cyl_bessel_j(0, A.LambdaR);
				A.j1LambdaR = std::cyl_bessel_j(1, A.LambdaR);
				loglambda += F.AbscissaSpacing;
			}
		}

		void dointegrals(const size_t& fi) {
			dointegrals_trapezoid(fi);
		}

		inline void dointegrals_trapezoid(const size_t& fi) {
			HankelTransforms& H = Hankel[fi];

			number_integrand_calls = 0;

			trapezoid(fi);//the results go into the variable trapezoid_result

			if (cmode() == CMode::FM) {
				H.I0.FM = trapezoid_result[0];
				H.I1.FM = trapezoid_result[1];
				H.I2.FM = trapezoid_result[2];
			}
			else if (cmode() == CMode::DC) {
				H.I0.dC = trapezoid_result[0];
				H.I1.dC = trapezoid_result[1];
				H.I2.dC = trapezoid_result[2];
			}
			else if (cmode() == CMode::DT) {
				H.I0.dT = trapezoid_result[0];
				H.I1.dT = trapezoid_result[1];
				H.I2.dT = trapezoid_result[2];
			}
			else if (cmode() == CMode::DZ) {
				H.I0.dZ = trapezoid_result[0];
				H.I1.dZ = trapezoid_result[1];
				H.I2.dZ = trapezoid_result[2];
			}
			else if (cmode() == CMode::DH) {
				H.I0.dH = trapezoid_result[0];
				H.I1.dH = trapezoid_result[1];
				H.I2.dH = trapezoid_result[2];
			}
			else if (cmode() == CMode::DR || cmode() == CMode::DX || cmode() == CMode::DY) {
				H.I0.dR = trapezoid_result[0];
				H.I1.dR = trapezoid_result[1];
				H.I2.dR = trapezoid_result[2];
			}
			else {
				glog.errormsg(_SRC_, "LE::dointegrals_trapezoid Calculation type %s not yet implemented\n", calculationtype.string().c_str());
			}
		}
		
		inline void trapezoid(const size_t& fi) {
			std::vector<cdouble> integrand1(3);
			std::vector<cdouble> integrand2(3);
			std::vector<cdouble> integrand3(3);

			trapezoid_result[0] = cdouble(0.0, 0.0);
			trapezoid_result[1] = cdouble(0.0, 0.0);
			trapezoid_result[2] = cdouble(0.0, 0.0);

			//First and last abscissa
			integrand(fi, 0);
			for (size_t ii = 0; ii < NumIntegrands; ii++) {
				trapezoid_result[ii] += integrand_result[ii];
			}

			integrand(fi, NumAbscissa - 1);
			for (size_t ii = 0; ii < NumIntegrands; ii++) {
				trapezoid_result[ii] += integrand_result[ii];
			}

			for (size_t ii = 0; ii < NumIntegrands; ii++) {
				trapezoid_result[ii] *= 0.5;
			}

			//Cenral Abscissas
			for (size_t ai = 1; ai < NumAbscissa - 1; ai++) {
				integrand(fi, ai);
				for (size_t ii = 0; ii < NumIntegrands; ii++) {
					trapezoid_result[ii] += integrand_result[ii];
				}
			}

			for (size_t ii = 0; ii < NumIntegrands; ii++) {
				trapezoid_result[ii] *= Frequencies[fi].AbscissaSpacing;
			}
		}
		
		inline void integrand(const size_t& fi, const size_t& ai) {
			AbscissaNode& A = Frequencies[fi].Abscissa[ai];
			number_integrand_calls++;

			double loopfactor = 1.0;
			if (ModellingLoopRadius > 0.0) {
				double lambda_a = A.Lambda * ModellingLoopRadius;
				loopfactor = 2.0 * std::cyl_bessel_j(1, lambda_a) / lambda_a;
			}

			double& lambdar = A.LambdaR;
			double& j0 = A.j0LambdaR;
			double& j1 = A.j1LambdaR;
			const double e = exp(-(Z + H) * A.Lambda);
			const double l2e = A.Lambda2 * e;
			const double l3e = A.Lambda3 * e;
			const double l4e = A.Lambda4 * e;

			cdouble k;
			switch (cmode()) {
			case CMode::FM:
				k = loopfactor * A.P21onP11;
				integrand_result[0] = k * l3e * j0;
				integrand_result[1] = k * l3e * j1;
				integrand_result[2] = k * l2e * j1;
				break;
			case CMode::DC:
				k = loopfactor * dP21onP11dCj(fi, ai, derivative_layer());
				integrand_result[0] = k * l3e * j0;
				integrand_result[1] = k * l3e * j1;
				integrand_result[2] = k * l2e * j1;
				break;
			case CMode::DT:
				k = loopfactor * dP21onP11dTj(fi, ai, derivative_layer());
				integrand_result[0] = k * l3e * j0;
				integrand_result[1] = k * l3e * j1;
				integrand_result[2] = k * l2e * j1;
				break;
			case CMode::DZ:
				k = loopfactor * A.P21onP11;
				integrand_result[0] = k * -l4e * j0;
				integrand_result[1] = k * -l4e * j1;
				integrand_result[2] = k * -l3e * j1;
				break;
			case CMode::DH:
				k = loopfactor * A.P21onP11;;
				integrand_result[0] = k * -l4e * j0;
				integrand_result[1] = k * -l4e * j1;
				integrand_result[2] = k * -l3e * j1;
				break;
			case CMode::DR:
			case CMode::DX:
			case CMode::DY:
				k = loopfactor * A.P21onP11;
				if (R != 0.0) {
					integrand_result[0] = k * (-l4e * j1);
					integrand_result[1] = k * (l4e * (j0 - j1 / lambdar));
					integrand_result[2] = k * (l3e * (j0 - j1 / lambdar));
				}
				else {
					integrand_result[0] = 0.0;
					integrand_result[1] = 0.0;
					integrand_result[2] = 0.0;
				}
				break;
			default:
				glog.errormsg(_SRC_, "LE::integrands Calculation type %s not yet implemented", calculationtype.string().c_str());
				break;
			}
		}

		void setverticaldipolefields(const size_t& fi) {
			setverticaldipoleprimaryfields();
			setverticaldipolesecondaryfields(fi);
		}

		void setverticaldipoleprimaryfields() {
			Fields.v.p.x = 0.0; Fields.v.p.y = 0.0; Fields.v.p.z = 0.0;

			if (BigR == 0)return;
			if (Source_Orientation.z() == 0.0)return;//ie no vertical dipole contribution

			if (cmode() == CMode::FM) {
				Fields.v.p.x = THREEONFOURPI<double>*X * (Z - H) / BigR5;
				Fields.v.p.y = THREEONFOURPI<double>*Y * (Z - H) / BigR5;
				Fields.v.p.z = THREEONFOURPI<double>*(Z - H) * (Z - H) / BigR5 - ONEONFOURPI<double> / BigR3;
			}
			else if (cmode() == CMode::DC || cmode() == CMode::DT) {
				Fields.v.p.x = 0.0;
				Fields.v.p.y = 0.0;
				Fields.v.p.z = 0.0;
			}
			else if (cmode() == CMode::DH) {
				Fields.v.p.x = 0.0;
				Fields.v.p.y = 0.0;
				Fields.v.p.z = 0.0;
			}
			else if (cmode() == CMode::DZ) {
				Fields.v.p.x = THREEONFOURPI<double>*X * (1.0 / BigR5 - 5.0 * (Z - H) * (Z - H) / BigR7);
				Fields.v.p.y = THREEONFOURPI<double>*Y * (1.0 / BigR5 - 5.0 * (Z - H) * (Z - H) / BigR7);
				Fields.v.p.z = THREEONFOURPI<double>*(3.0 * (Z - H) / BigR5 - 5.0 * (Z - H) * (Z - H) * (Z - H) / BigR7);
			}
			else if (cmode() == CMode::DX || cmode() == CMode::DY || cmode() == CMode::DR) {
				double dxdX = THREEONFOURPI<double>*(Z - H) * (1.0 / BigR5 - 5.0 * X * X / BigR7);
				double dydX = THREEONFOURPI<double>*Y * (Z - H) * -5.0 * X / BigR7;
				double dzdX = THREEONFOURPI<double>*(Z - H) * (Z - H) * -5.0 * X / BigR7 - ONEONFOURPI<double>*-3.0 * X / BigR5;

				double dxdY = THREEONFOURPI<double>*X * (Z - H) * -5.0 * Y / BigR7;
				double dydY = THREEONFOURPI<double>*(Z - H) * (1.0 / BigR5 - 5.0 * Y * Y / BigR7);
				double dzdY = THREEONFOURPI<double>*(Z - H) * (Z - H) * -5.0 * Y / BigR7 - ONEONFOURPI<double>*-3.0 * Y / BigR5;

				if (cmode() == CMode::DX) {
					double dXdXo = cosxyrotation;
					double dYdXo = -sinxyrotation;
					Fields.v.p.x = dxdX * dXdXo + dxdY * dYdXo;
					Fields.v.p.y = dydX * dXdXo + dydY * dYdXo;
					Fields.v.p.z = dzdX * dXdXo + dzdY * dYdXo;
				}
				else if (cmode() == CMode::DY) {
					double dXdYo = sinxyrotation;
					double dYdYo = cosxyrotation;
					Fields.v.p.x = dxdX * dXdYo + dxdY * dYdYo;
					Fields.v.p.y = dydX * dXdYo + dydY * dYdYo;
					Fields.v.p.z = dzdX * dXdYo + dzdY * dYdYo;
				}
				else if (cmode() == CMode::DR) {
					Fields.v.p.x = dxdX * XonR + dxdY * YonR;
					Fields.v.p.y = dydX * XonR + dydY * YonR;
					Fields.v.p.z = dzdX * XonR + dzdY * YonR;
				}
			}
			else {
				glog.errormsg(_SRC_, "LE::setverticaldipoleprimaryfields Calculation type %s not yet implemented", calculationtype.string().c_str());
			}
			unxyrotateandscale(Fields.v.p, Source_Orientation.z());
		}

		void setverticaldipolesecondaryfields(const size_t& fi)
		{
			Fields.v.s.x = cdouble(0.0, 0.0); Fields.v.s.y = cdouble(0.0, 0.0); Fields.v.s.z = cdouble(0.0, 0.0);

			if (Source_Orientation.z() == 0.0)return;//ie no vertical dipole contribution

			if (cmode() == CMode::FM) {
				Fields.v.s.x = -ONEONFOURPI<double> * XonR * Hankel[fi].I1.FM;
				Fields.v.s.y = -ONEONFOURPI<double> * YonR * Hankel[fi].I1.FM;
				Fields.v.s.z = -ONEONFOURPI<double> * Hankel[fi].I0.FM;
			}
			else if (cmode() == CMode::DC) {
				Fields.v.s.x = -ONEONFOURPI<double> *XonR * Hankel[fi].I1.dC;
				Fields.v.s.y = -ONEONFOURPI<double> *YonR * Hankel[fi].I1.dC;
				Fields.v.s.z = -ONEONFOURPI<double> *Hankel[fi].I0.dC;
			}
			else if (cmode() == CMode::DT) {
				Fields.v.s.x = -ONEONFOURPI<double> *XonR * Hankel[fi].I1.dT;
				Fields.v.s.y = -ONEONFOURPI<double> *YonR * Hankel[fi].I1.dT;
				Fields.v.s.z = -ONEONFOURPI<double> *Hankel[fi].I0.dT;
			}
			else if (cmode() == CMode::DH) {
				//these are negative of d/dz derivatives
				Fields.v.s.x = -ONEONFOURPI<double> *XonR * Hankel[fi].I1.dH;
				Fields.v.s.y = -ONEONFOURPI<double> *YonR * Hankel[fi].I1.dH;
				Fields.v.s.z = -ONEONFOURPI<double> *Hankel[fi].I0.dH;
			}
			else if (cmode() == CMode::DZ) {
				Fields.v.s.x = -ONEONFOURPI<double> *XonR * Hankel[fi].I1.dZ;
				Fields.v.s.y = -ONEONFOURPI<double> *YonR * Hankel[fi].I1.dZ;
				Fields.v.s.z = -ONEONFOURPI<double> *Hankel[fi].I0.dZ;
			}
			else if (cmode() == CMode::DX || cmode() == CMode::DY || cmode() == CMode::DR) {

				cdouble dxdX = 0.0; cdouble dydX = 0.0; cdouble dzdX = 0.0;
				cdouble dxdY = 0.0; cdouble dydY = 0.0; cdouble dzdY = 0.0;

				if (R != 0.0) {
					dxdX = -ONEONFOURPI<double> *(Hankel[fi].I1.FM * (1.0 / R - X * X / R3) + XonR * Hankel[fi].I1.dR * XonR);
					dydX = -ONEONFOURPI<double> *Y * (Hankel[fi].I1.FM * (-X / R3) + (1.0 / R) * Hankel[fi].I1.dR * XonR);
					dzdX = -ONEONFOURPI<double> *Hankel[fi].I0.dR * XonR;

					dxdY = -ONEONFOURPI<double> *X * (Hankel[fi].I1.FM * (-Y / R3) + (1.0 / R) * Hankel[fi].I1.dR * YonR);
					dydY = -ONEONFOURPI<double> *(Hankel[fi].I1.FM * (1.0 / R - Y * Y / R3) + YonR * Hankel[fi].I1.dR * YonR);
					dzdY = -ONEONFOURPI<double> *Hankel[fi].I0.dR * YonR;
				}

				if (cmode() == CMode::DX) {
					double dXdXo = cosxyrotation;
					double dYdXo = -sinxyrotation;
					Fields.v.s.x = dxdX * dXdXo + dxdY * dYdXo;
					Fields.v.s.y = dydX * dXdXo + dydY * dYdXo;
					Fields.v.s.z = dzdX * dXdXo + dzdY * dYdXo;
				}
				else if (cmode() == CMode::DY) {
					double dXdYo = sinxyrotation;
					double dYdYo = cosxyrotation;
					Fields.v.s.x = dxdX * dXdYo + dxdY * dYdYo;
					Fields.v.s.y = dydX * dXdYo + dydY * dYdYo;
					Fields.v.s.z = dzdX * dXdYo + dzdY * dYdYo;
				}
				else if (cmode() == CMode::DR) {
					Fields.v.s.x = dxdX * XonR + dxdY * YonR;
					Fields.v.s.y = dydX * XonR + dydY * YonR;
					Fields.v.s.z = dzdX * XonR + dzdY * YonR;
				}
			}
			else {
				glog.errormsg(_SRC_, "LE::setverticaldipolesecondaryfields Calculation type %s not yet implemented", calculationtype.string().c_str());
			}

			unxyrotateandscale(Fields.v.s, Source_Orientation.z());

		}
		
		void sethorizontaldipolefields(const size_t& fi) {
			sethorizontaldipoleprimaryfields();
			sethorizontaldipolesecondaryfields(fi);
		}
		
		void sethorizontaldipoleprimaryfields() {

			Fields.h.p.x = 0.0; Fields.h.p.y = 0.0; Fields.h.p.z = 0.0;

			if (BigR == 0)return;
			if (Source_Orientation.x() == 0.0 && Source_Orientation.y() == 0.0)return;//ie. not horizontal dipole contribution

			if (cmode() == CMode::FM) {
				Fields.h.p.x = THREEONFOURPI<double>*X * Y / BigR5;
				Fields.h.p.y = THREEONFOURPI<double>*Y * Y / BigR5 - ONEONFOURPI<double> / BigR3;
				Fields.h.p.z = THREEONFOURPI<double>*Y * (Z - H) / BigR5;
			}
			else if (cmode() == CMode::DC || cmode() == CMode::DT) {
				Fields.h.p.x = 0.0;
				Fields.h.p.y = 0.0;
				Fields.h.p.z = 0.0;
			}
			else if (cmode() == CMode::DH) {
				Fields.h.p.x = 0.0;
				Fields.h.p.y = 0.0;
				Fields.h.p.z = 0.0;
			}
			else if (cmode() == CMode::DZ) {
				Fields.h.p.x = THREEONFOURPI<double>*X * Y * (-5.0 * (Z - H) / BigR7);
				Fields.h.p.y = THREEONFOURPI<double>*Y * Y * (-5.0 * (Z - H) / BigR7) - ONEONFOURPI<double>*(-3.0 * (Z - H) / BigR5);
				Fields.h.p.z = THREEONFOURPI<double>*Y * (1.0 / BigR5 - 5.0 * (Z - H) * (Z - H) / BigR7);
			}
			else if (cmode() == CMode::DX || cmode() == CMode::DY || cmode() == CMode::DR) {
				double dxdX = THREEONFOURPI<double>*Y * (1.0 / BigR5 - 5.0 * X * X / BigR7);
				double dydX = THREEONFOURPI<double>*Y * Y * -5.0 * X / BigR7 + ONEONFOURPI<double>*3.0 * X / BigR5;
				double dzdX = THREEONFOURPI<double>*Y * (Z - H) * -5.0 * X / BigR7;

				double dxdY = THREEONFOURPI<double>*X * (1.0 / BigR5 - 5.0 * Y * Y / BigR7);
				double dydY = THREEONFOURPI<double>*(3.0 * Y / BigR5 - 5.0 * Y * Y * Y / BigR7);
				double dzdY = THREEONFOURPI<double>*(Z - H) * (1.0 / BigR5 - 5.0 * Y * Y / BigR7);

				if (cmode() == CMode::DX) {
					double dXdXo = cosxyrotation;
					double dYdXo = -sinxyrotation;
					Fields.h.p.x = dxdX * dXdXo + dxdY * dYdXo;
					Fields.h.p.y = dydX * dXdXo + dydY * dYdXo;
					Fields.h.p.z = dzdX * dXdXo + dzdY * dYdXo;
				}
				else if (cmode() == CMode::DY) {
					double dXdYo = sinxyrotation;
					double dYdYo = cosxyrotation;
					Fields.h.p.x = dxdX * dXdYo + dxdY * dYdYo;
					Fields.h.p.y = dydX * dXdYo + dydY * dYdYo;
					Fields.h.p.z = dzdX * dXdYo + dzdY * dYdYo;
				}
				else if (cmode() == CMode::DR) {
					Fields.h.p.x = dxdX * XonR + dxdY * YonR;
					Fields.h.p.y = dydX * XonR + dydY * YonR;
					Fields.h.p.z = dzdX * XonR + dzdY * YonR;
				}
			}
			else glog.errormsg(_SRC_, "LE::sethorizontaldipoleprimaryfields Calculation type %s not yet implemented", calculationtype.string().c_str());


			double scalefactor = std::hypot(Source_Orientation.x(), Source_Orientation.y());
			unxyrotateandscale(Fields.h.p, scalefactor);

		}
		
		void sethorizontaldipolesecondaryfields(const size_t& fi) {

			Fields.h.s.x = cdouble(0.0, 0.0); Fields.h.s.y = cdouble(0.0, 0.0); Fields.h.s.z = cdouble(0.0, 0.0);

			if (R == 0)return;
			if (Source_Orientation.x() == 0.0 && Source_Orientation.y() == 0.0)return;//ie. not horizontal dipole contribution

			if (cmode() == CMode::FM) {
				Fields.h.s.x = ONEONFOURPI<double> *(X * Y) / (R2) * (2.0 * Hankel[fi].I2.FM / R - Hankel[fi].I0.FM);
				Fields.h.s.y = ONEONFOURPI<double> *((Y * Y - X * X) * Hankel[fi].I2.FM / R3 - Y * Y * Hankel[fi].I0.FM / R2);
				Fields.h.s.z = ONEONFOURPI<double> *Y / R * Hankel[fi].I1.FM;
			}
			else if (cmode() == CMode::DC) {
				Fields.h.s.x = ONEONFOURPI<double> *(X * Y) / (R2) * (2.0 * Hankel[fi].I2.dC / R - Hankel[fi].I0.dC);
				Fields.h.s.y = ONEONFOURPI<double> *((Y * Y - X * X) * Hankel[fi].I2.dC / R3 - Y * Y * Hankel[fi].I0.dC / R2);
				Fields.h.s.z = ONEONFOURPI<double> *Y / R * Hankel[fi].I1.dC;
			}
			else if (cmode() == CMode::DT) {
				Fields.h.s.x = ONEONFOURPI<double> *(X * Y) / (R2) * (2.0 * Hankel[fi].I2.dT / R - Hankel[fi].I0.dT);
				Fields.h.s.y = ONEONFOURPI<double> *((Y * Y - X * X) * Hankel[fi].I2.dT / R3 - Y * Y * Hankel[fi].I0.dT / R2);
				Fields.h.s.z = ONEONFOURPI<double> *Y / R * Hankel[fi].I1.dT;
			}
			else if (cmode() == CMode::DH) {
				Fields.h.s.x = ONEONFOURPI<double> *(X * Y) / (R2) * (2.0 * Hankel[fi].I2.dH / R - Hankel[fi].I0.dH);
				Fields.h.s.y = ONEONFOURPI<double> *((Y * Y - X * X) * Hankel[fi].I2.dH / R3 - Y * Y * Hankel[fi].I0.dH / R2);
				Fields.h.s.z = ONEONFOURPI<double> *Y / R * Hankel[fi].I1.dH;
			}
			else if (cmode() == CMode::DZ) {
				Fields.h.s.x = ONEONFOURPI<double> *(X * Y) / (R2) * (2.0 * Hankel[fi].I2.dZ / R - Hankel[fi].I0.dZ);
				Fields.h.s.y = ONEONFOURPI<double> *((Y * Y - X * X) * Hankel[fi].I2.dZ / R3 - Y * Y * Hankel[fi].I0.dZ / R2);
				Fields.h.s.z = ONEONFOURPI<double> *Y / R * Hankel[fi].I1.dZ;
			}
			else if (cmode() == CMode::DX || cmode() == CMode::DY || cmode() == CMode::DR) {
				cdouble a, c, d, e, f, h;
				cdouble dadx, dcdx, dddx, dedx, dfdx, dhdx;
				cdouble dady, dcdy, dddy, dedy, dfdy, dhdy;

				cdouble b, dbdx, dbdy, dbdr;
				cdouble g, dgdx, dgdy;

				a = X * Y / R2;
				dadx = (R2 * Y - X * Y * 2.0 * X) / R4;
				dady = (R2 * X - X * Y * 2.0 * Y) / R4;

				b = 2.0 * Hankel[fi].I2.FM / R - Hankel[fi].I0.FM;
				dbdr = 2.0 * (R * Hankel[fi].I2.dR - Hankel[fi].I2.FM) / R2 - Hankel[fi].I0.dR;
				dbdx = dbdr * XonR;
				dbdy = dbdr * YonR;

				c = X * X / R3;
				dcdx = (R2 * 2.0 * X - X * X * 3.0 * X) / R5;
				dcdy = X * X * -3.0 * Y / R5;

				d = Y * Y / R3;
				dddx = Y * Y * -3.0 * X / R5;
				dddy = (R2 * 2.0 * Y - Y * Y * 3.0 * Y) / R5;

				e = d - c;
				dedx = dddx - dcdx;
				dedy = dddy - dcdy;

				f = Y * Y / R2;
				dfdx = Y * Y * -2.0 * X / R4;
				dfdy = (R2 * 2.0 * Y - Y * Y * 2.0 * Y) / R4;

				g = e * Hankel[fi].I2.FM - f * Hankel[fi].I0.FM;
				dgdx = e * Hankel[fi].I2.dR * XonR + Hankel[fi].I2.FM * dedx - (f * Hankel[fi].I0.dR * XonR + Hankel[fi].I0.FM * dfdx);
				dgdy = e * Hankel[fi].I2.dR * YonR + Hankel[fi].I2.FM * dedy - (f * Hankel[fi].I0.dR * YonR + Hankel[fi].I0.FM * dfdy);

				cdouble dxdX = ONEONFOURPI<double>*(a * dbdx + b * dadx);
				cdouble dxdY = ONEONFOURPI<double>*(a * dbdy + b * dady);

				cdouble dydX = ONEONFOURPI<double>*(dgdx);
				cdouble dydY = ONEONFOURPI<double>*(dgdy);

				h = Y / R;
				dhdx = -X * Y / R3;
				dhdy = 1.0 / R - Y * Y / R3;
				cdouble dzdX = ONEONFOURPI<double>*(Hankel[fi].I1.FM * dhdx + h * Hankel[fi].I1.dR * XonR);
				cdouble dzdY = ONEONFOURPI<double>*(Hankel[fi].I1.FM * dhdy + h * Hankel[fi].I1.dR * YonR);

				if (cmode() == CMode::DX) {
					double dXdXo = cosxyrotation;
					double dYdXo = -sinxyrotation;
					Fields.h.s.x = dxdX * dXdXo + dxdY * dYdXo;
					Fields.h.s.y = dydX * dXdXo + dydY * dYdXo;
					Fields.h.s.z = dzdX * dXdXo + dzdY * dYdXo;
				}
				else if (cmode() == CMode::DY) {
					double dXdYo = sinxyrotation;
					double dYdYo = cosxyrotation;
					Fields.h.s.x = dxdX * dXdYo + dxdY * dYdYo;
					Fields.h.s.y = dydX * dXdYo + dydY * dYdYo;
					Fields.h.s.z = dzdX * dXdYo + dzdY * dYdYo;
				}
				else if (cmode() == CMode::DR) {
					Fields.h.s.x = dxdX * XonR + dxdY * YonR;
					Fields.h.s.y = dydX * XonR + dydY * YonR;
					Fields.h.s.z = dzdX * XonR + dzdY * YonR;
				}
			}
			else glog.errormsg(_SRC_, "LE::sethorizontaldipolesecondaryfields Calculation type %lu not yet implemented", calculationtype.string().c_str());

			double scalefactor = std::hypot(Source_Orientation.x(), Source_Orientation.y());
			unxyrotateandscale(Fields.h.s, scalefactor);
		}

		//Horizontal coplanar
		cdouble ppmHCP(const size_t& fi) {
			return 1.0e6 * (R3 * Hankel[fi].I0.FM);
		}
		cdouble dppmHCPdC(const size_t& fi) {
			return 1.0e6 * (R3 * Hankel[fi].I0.dC);
		}
		cdouble dppmHCPdT(const size_t& fi) {
			return 1.0e6 * (R3 * Hankel[fi].I0.dT);
		}
		cdouble dppmHCPdZ(const size_t& fi) {
			return 1.0e6 * (R3 * Hankel[fi].I0.dZ);
		}
		cdouble dppmHCPdH(const size_t& fi) {
			return 1.0e6 * (R3 * Hankel[fi].I0.dH);
		}
		cdouble dppmHCPdR(const size_t& fi) {
			return 1.0e6 * (R3 * Hankel[fi].I0.dR + 3.0 * R2 * Hankel[fi].I0.FM);
		}

		//Perpendicular
		cdouble ppmPER(const size_t& fi) {
			return -1.0e6 * (1.0 - 0.5 * R3 * Hankel[fi].I1.FM);
		}
		cdouble dppmPERdC(const size_t& fi) {
			return 0.5e6 * (R3 * Hankel[fi].I1.dC);
		}
		cdouble dppmPERdT(const size_t& fi) {
			return 0.5e6 * (R3 * Hankel[fi].I1.dT);
		}
		cdouble dppmPERdZ(const size_t& fi) {
			return 0.5e6 * (R3 * Hankel[fi].I1.dZ);
		}
		cdouble dppmPERdH(const size_t& fi) {
			return 0.5e6 * (R3 * Hankel[fi].I1.dH);
		}
		cdouble dppmPERdR(const size_t& fi) {
			return 0.5e6 * (R3 * Hankel[fi].I1.dR + 3.0 * R2 * Hankel[fi].I1.FM);
		}

		//Vertical coaxial
		cdouble ppmVCX(const size_t& fi) {
			return -0.5e6 * (R2 * Hankel[fi].I2.FM - R3 * Hankel[fi].I0.FM);
		}
		cdouble dppmVCXdC(const size_t& fi) {
			return -0.5e6 * (R2 * Hankel[fi].I2.dC - R3 * Hankel[fi].I0.dC);
		}
		cdouble dppmVCXdT(const size_t& fi) {
			return -0.5e6 * (R2 * Hankel[fi].I2.dT - R3 * Hankel[fi].I0.dT);
		}
		cdouble dppmVCXdZ(const size_t& fi) {
			return -0.5e6 * (R2 * Hankel[fi].I2.dZ - R3 * Hankel[fi].I0.dZ);
		}
		cdouble dppmVCXdH(const size_t& fi) {
			return -0.5e6 * (R2 * Hankel[fi].I2.dH - R3 * Hankel[fi].I0.dH);
		}
		cdouble dppmVCXdR(const size_t& fi) {
			return -0.5e6 * (R2 * Hankel[fi].I2.dR + 2.0 * R * Hankel[fi].I2.FM - R3 * Hankel[fi].I0.dR - 3.0 * R2 * Hankel[fi].I0.FM);
		}

		//Vertical coplanar
		cdouble ppmVCP(const size_t& fi) {
			return 1.0e6 * (R2 * Hankel[fi].I2.FM);
		}
		cdouble dppmVCPdC(const size_t& fi) {
			return 1.0e6 * (R2 * Hankel[fi].I2.dC);
		}
		cdouble dppmVCPdT(const size_t& fi) {
			return 1.0e6 * (R2 * Hankel[fi].I2.dT);
		}
		cdouble dppmVCPdZ(const size_t& fi) {
			return 1.0e6 * (R2 * Hankel[fi].I2.dZ);
		}
		cdouble dppmVCPdH(const size_t& fi) {
			return 1.0e6 * (R2 * Hankel[fi].I2.dH);
		}
		cdouble dppmVCPdR(const size_t& fi) {
			return 1.0e6 * (R2 * Hankel[fi].I2.dR + 2.0 * R * Hankel[fi].I2.FM);
		}
	};
};
};
