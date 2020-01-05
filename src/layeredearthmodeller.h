/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef layeredearthmodeller_H
#define layeredearthmodeller_H

#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>

#include "logger.h"
#include "general_constants.h"

//Formulation mainly from the book 
//Geo-Electromagnetism, Wait James, R. Academic Press 1982.

#include <Eigen/Dense>
typedef Eigen::Vector3d  Vec3d;
typedef Eigen::Vector3cd Vec3cd;
typedef Eigen::Matrix2cd  sPropogationMatrix;
typedef Eigen::Matrix3d  Mat3d;
typedef Eigen::Matrix3cd Mat3cd;
typedef std::complex<double> cdouble;
typedef std::vector<std::complex<double>> cvector;

struct sHankelTransform {
	cdouble FM;
	cdouble dC;
	cdouble dT;
	cdouble dX;
	cdouble dY;
	cdouble dZ;
	cdouble dH;
};

struct sAbscissaLayerNode {
	cdouble   U;
	cdouble   Exp2UT;
	sPropogationMatrix LayerMatrix;
	sPropogationMatrix LayerPreMatrix;
	sPropogationMatrix LayerPostMatrix;
};

struct sAbscissaNode {
	double Lambda   = 0.0;
	double Lambda2  = 0.0;
	double Lambda3  = 0.0;
	double Lambda4  = 0.0;
	double Lambda_r = 0.0;
	double j0Lambda_r = 0.0;
	double j1Lambda_r = 0.0;
	std::vector<sAbscissaLayerNode> Layer;
	sPropogationMatrix P_Full;
	cdouble P21onP11;
};

struct sLayerNode {
	double Conductivity;
	double Thickness;
};

class cLayeredEarthModeller {

public:

	enum class CalculationType { FM, DX, DY, DZ, DH, DC, DT, DB };
	
private:

	double X = 0.0, Y = 0.0, Z = 0.0, H = 0.0;
	double X2 = 0, Y2 = 0, X4 = 0, Y4 = 0, ZH = 0, ZH2 = 0;
	double r = 0, r2 = 0, r3 = 0, r5 = 0, r4 = 0, R = 0, R2 = 0, R5 = 0, R7 = 0;

	bool geometrychanged;
	bool earthchanged;
	
	std::vector<sAbscissaNode> Abscissa;	
	std::vector<double> Conductivity;
	std::vector<double> Thickness;
	
	
	//Hankle Stuff    
	sHankelTransform I0;
	sHankelTransform I1;
	sHankelTransform I2;

	cdouble trapezoid_result[3];
	cdouble integrand_result[3];

	double meanconductivity;
	double meanlog10conductivity;

	double LowerFractionalWidth;
	double UpperFractionalWidth;

public:

	double Frequency;
	double Omega;
	double MuZeroOmega;
	cdouble iMuZeroOmega;
	double ApproximateHalfspace;
	double PeakLambda;
	double LowerBound;
	double UpperBound;
	double AbscissaSpacing;

	cLayeredEarthModeller() { initialise(); };

	~cLayeredEarthModeller() {  };

	
	size_t na() const { return Abscissa.size(); }
	size_t nl() const { return Conductivity.size(); }
	const std::vector<double>& getConductivity() {
		return Conductivity;
	};

	const std::vector<double>& getThickness() {
		return Thickness;
	};

	void initialise()
	{				
		LowerFractionalWidth = 4.44;
		UpperFractionalWidth = 1.84;
		//Ground EM 31 and 38
		//NumLayers()     = 0;  
		//NumAbscissa   = 1581;
		//LowerFractionalWidth  = 4.44;
		//UpperFractionalWidth  = 10.84;
		setnumabscissa(17);
	};

	void setupcomputations(){
		if(earthchanged || geometrychanged){
			setfrequencyabscissalayers();
			earthchanged = false;
			geometrychanged = false;
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
		for (size_t i = 0; i < nl() ; i++) {
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

		if (update_r || update_zh){
			R2 = (X2 + Y2 + ZH2);
			R  = std::sqrt(R2);
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
		Omega = TWOPI * frequency;
		MuZeroOmega = MUZERO * Omega;
		iMuZeroOmega = cdouble(0.0, MuZeroOmega);
	};

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

	cdouble rzero_recursive(const double lambda)
	{
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

	inline cdouble rzero_propogationmatrix(const size_t ai)
	{
		//Oldenberg's propogation matrix formulation
		setlayermatrices(ai);
		setpmatrix(ai);
		return Abscissa[ai].P21onP11;
	};

	void setlayermatrices(const size_t ai)
	{
		cdouble e, eh, e1, e2;

		sAbscissaNode& A = Abscissa[ai];

		//M1  
		e = A.Layer[0].U / A.Lambda;
		eh = e / 2.0;
		e1 = 0.5 + eh;
		e2 = 0.5 - eh;

		sPropogationMatrix M;
		M(0,0) = e1;
		M(0,1) = e2;
		M(1,0) = e2;
		M(1,1) = e1;
		A.Layer[0].LayerMatrix = M;

		for (size_t li = 1; li < nl(); li++) {
			//assumes all pearmabilities are muzero
			e = A.Layer[li].U / A.Layer[li - 1].U;
			eh = e / 2.0;
			e1 = 0.5 + eh;
			e2 = 0.5 - eh;
			M(0,0) = e1;
			M(0,1) = e2;
			M(1,0) = e2 * A.Layer[li - 1].Exp2UT;
			M(1,1) = e1 * A.Layer[li - 1].Exp2UT;

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

		for (int li = (int) nl() - 3; li >= 0; li--) {			
			A.Layer[li].LayerPostMatrix = A.Layer[(size_t)li + 1].LayerMatrix * A.Layer[(size_t)li + 1].LayerPostMatrix;
		}

	};

	void setpmatrix(const size_t ai)
	{
		sAbscissaNode& A = Abscissa[ai];
		//Set Full matrix
		if (nl() == 1) {
			A.P_Full = A.Layer[0].LayerMatrix;
		}
		else {
			A.P_Full = Abscissa[ai].Layer[nl() - 1].LayerPreMatrix * Abscissa[ai].Layer[nl() - 1].LayerMatrix;
		}
		A.P21onP11 = A.P_Full(1,0) / A.P_Full(0,0);
	}

	sPropogationMatrix dMjdCj(const size_t ai, const size_t li)
	{
		sPropogationMatrix m;

		if (li == 0) {
			const cdouble a = iMuZeroOmega / (4.0 * Abscissa[ai].Lambda * Abscissa[ai].Layer[li].U);
			m(0,0) = a;
			m(0,1) = -a;
			m(1,0) = -a;
			m(1,1) = a;
			return m;
		}
		else {
			const cdouble a  = iMuZeroOmega / (4.0 * Abscissa[ai].Layer[li - 1].U * Abscissa[ai].Layer[li].U);
			const cdouble ae = a * Abscissa[ai].Layer[li - 1].Exp2UT;
			m(0,0) = a;
			m(0,1) = -a;
			m(1,0) = -ae;
			m(1,1) = ae;
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
		M(0,0) = 0.5 * dyds;
		M(0,1) = -M(0,0);
		M(1,0) = 0.5 * (dvds - (ydvds + vdyds));
		M(1,1) = 0.5 * (dvds + (ydvds + vdyds));

		return M;

	}

	sPropogationMatrix dPdCj(const size_t ai, const size_t li)
	{
		sAbscissaNode& A = Abscissa[ai];
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
				return A.Layer[li].LayerPreMatrix * c *	A.Layer[li + 1].LayerPostMatrix;
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
		return m(1,0) / Abscissa[ai].P_Full(0,0) - m(0,0) * Abscissa[ai].P21onP11 / Abscissa[ai].P_Full(0,0);
	}

	sPropogationMatrix dMjplus1dTj(const size_t ai, const size_t li)
	{
		cdouble y = Abscissa[ai].Layer[li + 1].U / Abscissa[ai].Layer[li].U;
		cdouble dvdt = -2.0 * Abscissa[ai].Layer[li].U * Abscissa[ai].Layer[li].Exp2UT;

		sPropogationMatrix m;
		m(0,0) = 0.0;
		m(0,1) = 0.0;
		m(1,0) = 0.5 * (1.0 - y) * dvdt;
		m(1,1) = 0.5 * (1.0 + y) * dvdt;
		return m;

	}

	sPropogationMatrix dPdTj(const size_t ai, const size_t li)
	{
		sAbscissaNode& A = Abscissa[ai];

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
			glog.errormsg(_SRC_,"Zero thickness derivative for halfspace layer\n");
			return sPropogationMatrix::Zero();			
		}
	}

	cdouble dP21onP11dTj(const size_t ai, const size_t li)
	{
		sPropogationMatrix m = dPdTj(ai, li);
		return m(1,0) / Abscissa[ai].P_Full(0,0) - m(0,0) * Abscissa[ai].P21onP11 / Abscissa[ai].P_Full(0,0);
	}

	void setintegrationnodes()
	{
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
			sAbscissaNode& A = Abscissa[ai];
			lambda = std::exp(loglambda);
			A.Lambda = lambda;
			A.Lambda2 = A.Lambda * lambda;
			A.Lambda3 = A.Lambda2 * lambda;
			A.Lambda4 = A.Lambda3 * lambda;
			A.Lambda_r = A.Lambda * r;			
			A.j0Lambda_r = j0(A.Lambda_r);
			A.j1Lambda_r = j1(A.Lambda_r);
			loglambda += AbscissaSpacing;
		}
	}

	void dointegrals(const CalculationType& calculationtype, const size_t& derivativelayer)
	{
		trapezoid(calculationtype, derivativelayer);//the results go into the variable trapezoid_result			

		switch (calculationtype)
		{
		case CalculationType::FM:
			I0.FM = trapezoid_result[0];
			I1.FM = trapezoid_result[1];
			I2.FM = trapezoid_result[2];
			break;
		case CalculationType::DC:
			I0.dC = trapezoid_result[0];
			I1.dC = trapezoid_result[1];
			I2.dC = trapezoid_result[2];
			break;
		case CalculationType::DT:
			I0.dT = trapezoid_result[0];
			I1.dT = trapezoid_result[1];
			I2.dT = trapezoid_result[2];
			break;
		case CalculationType::DZ:
			I0.dZ = trapezoid_result[0];
			I1.dZ = trapezoid_result[1];
			I2.dZ = trapezoid_result[2];
			break;
		case CalculationType::DH:
			I0.dH = trapezoid_result[0];
			I1.dH = trapezoid_result[1];
			I2.dH = trapezoid_result[2];
			break;
		case CalculationType::DX:
			I0.dX = trapezoid_result[0];
			I1.dX = trapezoid_result[1];
			I2.dX = trapezoid_result[2];
			break;
		case CalculationType::DY:
			I0.dY = trapezoid_result[0];
			I1.dY = trapezoid_result[1];
			I2.dY = trapezoid_result[2];
			break;
		default:
			glog.logmsg("Error: dointegrals() unknown calculation type %c\n", calculationtype);
			exit(1);
		}
	}

	void trapezoid(const CalculationType& calculationtype, const size_t& derivativelayer)
	{
		trapezoid_result[0] = cdouble(0.0, 0.0);
		trapezoid_result[1] = cdouble(0.0, 0.0);
		trapezoid_result[2] = cdouble(0.0, 0.0);

		//First and last abscissa
		compute_integrand(0,calculationtype,derivativelayer);
		trapezoid_result[0] += integrand_result[0];
		trapezoid_result[1] += integrand_result[1];
		trapezoid_result[2] += integrand_result[2];

		compute_integrand(na() - 1, calculationtype, derivativelayer);
		trapezoid_result[0] += integrand_result[0];
		trapezoid_result[1] += integrand_result[1];
		trapezoid_result[2] += integrand_result[2];

		trapezoid_result[0] *= 0.5;
		trapezoid_result[1] *= 0.5;
		trapezoid_result[2] *= 0.5;


		//Cenral Abscissas
		for (size_t ai = 1; ai < na() - 1; ai++) {
			compute_integrand(ai, calculationtype, derivativelayer);
			trapezoid_result[0] += integrand_result[0];
			trapezoid_result[1] += integrand_result[1];
			trapezoid_result[2] += integrand_result[2];
		}

		trapezoid_result[0] *= AbscissaSpacing;
		trapezoid_result[1] *= AbscissaSpacing;
		trapezoid_result[2] *= AbscissaSpacing;
	}

	void compute_integrand(const size_t& ai, const CalculationType& calculationtype, const size_t& derivativelayer)
	{
		sAbscissaNode& A = Abscissa[ai];

		double lambdar = A.Lambda_r;
		double j0 = A.j0Lambda_r;
		double j1 = A.j1Lambda_r;
		double e = exp(-(Z + H) * A.Lambda);
		double l2e = A.Lambda2 * e;
		double l3e = A.Lambda3 * e;
		double l4e = A.Lambda4 * e;


		double k0, k1, k2;

		cdouble earthkernel;

		switch (calculationtype)
		{
		case CalculationType::FM:
			earthkernel = -A.P21onP11;
			k0 = l3e * j0;
			k1 = l3e * j1;
			k2 = l2e * j1;
			break;
		case CalculationType::DX:
			earthkernel = -A.P21onP11;
			k0 = -l4e * j1 * X / r;
			k1 = l4e * (j0 - j1 / lambdar) * X / r;
			k2 = l3e * (j0 - j1 / lambdar) * X / r;
			break;
		case CalculationType::DY:
			earthkernel = -A.P21onP11;
			k0 = -l4e * j1 * Y / r;
			k1 = l4e * (j0 - j1 / lambdar) * Y / r;
			k2 = l3e * (j0 - j1 / lambdar) * Y / r;
			break;
		case CalculationType::DZ:
			earthkernel = -A.P21onP11;
			k0 = -l4e * j0;
			k1 = -l4e * j1;
			k2 = -l3e * j1;
			break;
		case CalculationType::DH:
			earthkernel = -A.P21onP11;;
			k0 = -l4e * j0;
			k1 = -l4e * j1;
			k2 = -l3e * j1;
			break;
		case CalculationType::DC:
			earthkernel = -dP21onP11dCj(ai, derivativelayer);
			k0 = l3e * j0;
			k1 = l3e * j1;
			k2 = l2e * j1;
			break;
		case CalculationType::DT:
			earthkernel = -dP21onP11dTj(ai, derivativelayer);
			k0 = l3e * j0;
			k1 = l3e * j1;
			k2 = l2e * j1;
			break;
		default:
			glog.logmsg("Error: compute_integrand() unknown calculation type %c\n", calculationtype);
			exit(1);
		}

		integrand_result[0] = earthkernel * k0;
		integrand_result[1] = earthkernel * k1;
		integrand_result[2] = earthkernel * k2;
	}

	Mat3d PTFM()
	{
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
		return T;
	}

	Mat3d dPTdX()
	{
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
		return T;
	}

	Mat3d dPTdY()
	{
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
		return T;
	}

	Mat3d dPTdZ()
	{
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
		return T;
	}

	Mat3d dPTdH()
	{
		Mat3d T;
		//of course this is just minus d/dZ		
		T(0, 0) = 3.0 * ZH * (4.0 * X2 - Y2 - ZH2) / R7;
		T(0, 1) = 15.0 * X * Y / R7 * ZH;
		T(0, 2) = -3.0 * X * (X2 + Y2 - 4.0 * ZH2) / R7;
		T(1, 0) = T(0, 1);
		T(1, 1) = -3.0 * ZH * (X2 - 4.0 * Y2 + ZH2) / R7;
		T(1, 2) = -3.0 * Y * (X2 + Y2 - 4.0 * ZH2) / R7;
		T(2, 0) = T(0, 2);
		T(2, 1) = T(1, 2);
		T(2, 2) = -3.0 * ZH * (3.0 * X2 + 3.0 * Y2 - 2.0 * ZH2) / R7;
		return T;
	}

	Mat3cd STFM()
	{
		Mat3cd T;
		T(0, 0) = ((X2 / r2 - Y2 / r2) * I2.FM / r - I0.FM * X2 / r2);
		T(0, 1) = (X * Y / r2) * (2.0 * I2.FM / r - I0.FM);
		T(0, 2) = (-X / r) * I1.FM;

		//Note error in Fitterman and Yin paper should not be minus sign at element 2,1
		T(1, 0) = T(0, 1);
		T(1, 1) = ((Y2 / r2 - X2 / r2) * I2.FM / r - I0.FM * Y2 / r2);
		T(1, 2) = (-Y / r) * I1.FM;

		T(2, 0) = -T(0, 2);
		T(2, 1) = -T(1, 2);
		T(2, 2) = -I0.FM;
		return T;
	}

	Mat3cd dSTdC()
	{
		Mat3cd T;
		T(0, 0) = ((X2 / r2 - Y2 / r2) * I2.dC / r - I0.dC * X2 / r2);
		T(0, 1) = (X * Y / r2) * (2.0 * I2.dC / r - I0.dC);
		T(0, 2) = (-X / r) * I1.dC;

		T(1, 0) = T(0, 1);
		T(1, 1) = ((Y2 / r2 - X2 / r2) * I2.dC / r - I0.dC * Y2 / r2);
		T(1, 2) = (-Y / r) * I1.dC;

		T(2, 0) = -T(0, 2);
		T(2, 1) = -T(1, 2);
		T(2, 2) = -I0.dC;
		return T;
	}

	Mat3cd dSTdT()
	{
		Mat3cd T;
		T(0, 0) = ((X2 / r2 - Y2 / r2) * I2.dT / r - I0.dT * X2 / r2);
		T(0, 1) = (X * Y / r2) * (2.0 * I2.dT / r - I0.dT);
		T(0, 2) = (-X / r) * I1.dT;

		T(1, 0) = T(0, 1);
		T(1, 1) = ((Y2 / r2 - X2 / r2) * I2.dT / r - I0.dT * Y2 / r2);
		T(1, 2) = (-Y / r) * I1.dT;

		T(2, 0) = -T(0, 2);
		T(2, 1) = -T(1, 2);
		T(2, 2) = -I0.dT;
		return T;
	}

	Mat3cd dSTdX()
	{
		Mat3cd T;
		T(0, 0) = (X4 * I2.dX - I0.dX * X4 * r - I2.FM * X2 * X - I0.dX * X2 * r * Y2 - 2.0 * I0.FM * X * r * Y2 + 5.0 * X * Y2 * I2.FM - Y4 * I2.dX) / r5;
		T(0, 1) = 2.0 * Y / r3 * I2.FM - Y / r2 * I0.FM - 6.0 * X2 * Y / r5 * I2.FM + 2.0 * X2 * Y / r4 * I0.FM + 2.0 * X * Y / r3 * I2.dX - X * Y / r2 * I0.dX;
		T(0, 2) = -(I1.FM * Y2 + X2 * X * I1.dX + X * I1.dX * Y2) / r3;

		T(1, 0) = T(0, 1);
		T(1, 1) = -(I2.dX * X4 - I2.FM * X2 * X + I0.dX * Y2 * r * X2 - 2.0 * I0.FM * Y2 * X * r + 5.0 * X * Y2 * I2.FM + I0.dX * Y4 * r - Y4 * I2.dX) / r5;
		T(1, 2) = Y / r3 * I1.FM * X - Y / r * I1.dX;

		T(2, 0) = -T(0, 2);
		T(2, 1) = -T(1, 2);
		T(2, 2) = -I0.dX;
		return T;
	}

	Mat3cd dSTdY()
	{
		Mat3cd T;
		T(0, 0) = (X4 * I2.dY - I0.dY * X4 * r + 2.0 * I0.FM * X2 * Y * r - I0.dY * X2 * r * Y2 - 5.0 * X2 * Y * I2.FM + Y2 * Y * I2.FM - Y4 * I2.dY) / r5;
		T(0, 1) = 2.0 * X / r3 * I2.FM - X / r2 * I0.FM - 6.0 * X * Y2 / r5 * I2.FM + 2.0 * X * Y2 / r4 * I0.FM + 2.0 * X * Y / r3 * I2.dY - X * Y / r2 * I0.dY;
		T(0, 2) = X / r3 * I1.FM * Y - X / r * I1.dY;

		T(1, 0) = T(0, 1);
		T(1, 1) = -(I2.dY * X4 + I0.dY * Y2 * r * X2 + 2.0 * I0.FM * Y * r * X2 - 5.0 * X2 * Y * I2.FM + I0.dY * Y4 * r + Y2 * Y * I2.FM - Y4 * I2.dY) / r5;
		T(1, 2) = -(I1.FM * X2 + Y * I1.dY * X2 + Y2 * Y * I1.dY) / r3;

		T(2, 0) = -T(0, 2);
		T(2, 1) = -T(1, 2);
		T(2, 2) = -I0.dY;
		return T;
	}

	Mat3cd dSTdZ()
	{
		Mat3cd T;
		T(0, 0) = ((X2 / r2 - Y2 / r2) * I2.dZ / r - I0.dZ * X2 / r2);
		T(0, 1) = X * Y / r2 * (2.0 * I2.dZ / r - I0.dZ);
		T(0, 2) = -X / r * I1.dZ;

		T(1, 0) = T(0, 1);
		T(1, 1) = ((Y2 / r2 - X2 / r2) * I2.dZ / r - I0.dZ * Y2 / r2);
		T(1, 2) = -Y / r * I1.dZ;

		T(2, 0) = -T(0, 2);
		T(2, 1) = -T(1, 2);
		T(2, 2) = -I0.dZ;
		return T;
	}

	Mat3cd dSTdH()
	{
		Mat3cd T;
		T(0, 0) = ((X2 / r2 - Y2 / r2) * I2.dH / r - I0.dH * X2 / r2);
		T(0, 1) = X * Y / r2 * (2.0 * I2.dH / r - I0.dH);
		T(0, 2) = -X / r * I1.dH;

		T(1, 0) = T(0, 1);
		T(1, 1) = ((Y2 / r2 - X2 / r2) * I2.dH / r - I0.dH * Y2 / r2);
		T(1, 2) = -Y / r * I1.dH;

		T(2, 0) = -T(0, 2);
		T(2, 1) = -T(1, 2);
		T(2, 2) = -I0.dH;
		return T;
	}

	Mat3d PrimaryTensor(const CalculationType& calculationtype)
	{		
		switch (calculationtype)
		{
		case CalculationType::FM: return PTFM();
		case CalculationType::DX: return dPTdX();
		case CalculationType::DY: return dPTdY();
		case CalculationType::DZ: return dPTdZ();
		case CalculationType::DH: return dPTdH();
		case CalculationType::DC: return Mat3d::Zero();
		case CalculationType::DT: return Mat3d::Zero();
		default:
			glog.errormsg(_SRC_,"Unknown calculation type %c\n", calculationtype);
		}
	}

	Mat3cd SecondaryTensor(const CalculationType& calculationtype, const size_t& derivativelayer)
	{		
		dointegrals(calculationtype,derivativelayer);

		switch (calculationtype)
		{
		case CalculationType::FM: return STFM();
		case CalculationType::DX: return dSTdX();
		case CalculationType::DY: return dSTdY();
		case CalculationType::DZ: return dSTdZ();
		case CalculationType::DH: return dSTdH();
		case CalculationType::DC: return dSTdC();
		case CalculationType::DT: return dSTdT();
		default:
			glog.errormsg(_SRC_,"Unknown calculation type %c\n", calculationtype);
		}
	}

	Vec3d primary_inertial_frame(const CalculationType& calculationtype, const Vec3d& txdir)
	{		
		return PrimaryTensor(calculationtype) * txdir;		
	}

	Vec3cd secondary_inertial_frame(const CalculationType& calculationtype, const size_t& derivativelayer, const Vec3d& txdir)
	{		
		return SecondaryTensor(calculationtype, derivativelayer) * txdir;		
	}

	double primary(const Vec3d& txdir, const Vec3d& rxdir)
	{
		Vec3d v = primary_inertial_frame(CalculationType::FM, txdir);
		return v.dot(rxdir);
	}

	cdouble secondary(const Vec3d& txdir, const Vec3d& rxdir)
	{
		Vec3cd v = secondary_inertial_frame(CalculationType::FM, 0, txdir);
		return v.dot(rxdir);
	}

	double  dp(const CalculationType& calculationtype, const Vec3d& txdir, const Vec3d& rxdir)
	{
		Vec3d v = primary_inertial_frame(calculationtype, txdir);
		return v.dot(rxdir);
	}

	cdouble ds(const CalculationType& calculationtype, const size_t& derivativelayer, const Vec3d& txdir, const Vec3d& rxdir)
	{
		Vec3cd v = secondary_inertial_frame(calculationtype, derivativelayer, txdir);
		return v.dot(rxdir);
	}

	cdouble dsdx(const Vec3d& txdir, const Vec3d& rxdir)
	{
		Vec3cd sf = secondary_inertial_frame(CalculationType::DX, 0, txdir);
		return sf.dot(rxdir);
	}

	cdouble dsdy(const Vec3d& txdir, const Vec3d& rxdir)
	{
		Vec3cd sf = secondary_inertial_frame(CalculationType::DY, 0, txdir);
		return sf.dot(rxdir);
	}

	cdouble dsdz(const Vec3d& txdir, const Vec3d& rxdir)
	{
		Vec3cd sf = secondary_inertial_frame(CalculationType::DZ, 0, txdir);
		return sf.dot(rxdir);
	}

	cdouble dsdh(const Vec3d& txdir, const Vec3d& rxdir)
	{
		Vec3cd sf = secondary_inertial_frame(CalculationType::DH, 0, txdir);
		return sf.dot(rxdir);
	}

	cdouble dsdc(const size_t dlayer, const Vec3d& txdir, const Vec3d& rxdir)
	{
		Vec3cd sf = secondary_inertial_frame(CalculationType::DC, dlayer, txdir);
		return sf.dot(rxdir);
	}

	cdouble dsdt(const size_t dlayer, const Vec3d& txdir, const Vec3d& rxdir)
	{
		Vec3cd sf = secondary_inertial_frame(CalculationType::DT, dlayer, txdir);
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
		double  pf  = primary(txdir, rxdir);
		cdouble dsf = ds(calculationtype, derivativelayer, txdir, rxdir);
		dsf = std::complex<double>(dsf.real(), -dsf.imag());
		return 1.0e6 * (dsf / pf);
	}

};
#endif

