/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef _lem_H_
#define _lem_H_

#include <cmath>
#include <complex>
#include <vector>
#include "geometry3d.h"
#include "general_constants.h"
#include "general_utils.h"
#include "earth1d.h"

//Formulation mainly from the book 
//Geo-Electromagnetism, Wait James, R. Academic Press 1982

typedef std::complex<double> cdouble;

struct HankelTransform{
	cdouble FM = 0.0;
	cdouble dC = 0.0;
	cdouble dT = 0.0;
	cdouble dZ = 0.0;
	cdouble dH = 0.0;
	cdouble dR = 0.0;
};

struct HankelTransforms{
	HankelTransform I0;
	HankelTransform I1;
	HankelTransform I2;
};

struct RealField{
	double x=0.0;
	double y=0.0;
	double z=0.0;
};

struct ComplexField{
	cdouble x=0.0;
	cdouble y=0.0;
	cdouble z=0.0;
};

struct TotalField{
	RealField    p = {};//primary filed component
	ComplexField s = {};//secondary field component  
};

struct ResponseField{
	TotalField   v = {};//due to vertical dipole
	TotalField   h = {};//due to y directed horizontal dipole
	TotalField   t = {};//due to total dipole
};

class PropogationMatrix{

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

struct AbscissaLayerNode{
	cdouble   U = 0.0;
	cdouble   Exp2UT = 0.0;
	PropogationMatrix LayerMatrix;
	PropogationMatrix LayerPreMatrix;
	PropogationMatrix LayerPostMatrix;
};

struct AbscissaNode{
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

struct FrequencyNode{
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

struct LayerNode{
	double Thickness = 0.0;  // m
	double Conductivity = 0.0; // S/m
	double Chargeability = 0.0; // 
	double TimeConstant = 0.0; // s
	double FrequencyDependence = 0.0;	//unit less
};


class cLEM {

public:

	enum class RZeroMethod {
		PROPOGATIONMATRIX,
		RECURSIVE
	};

	enum class CalculationType {
		FORWARDMODEL,
		CONDUCTIVITYDERIVATIVE,
		THICKNESSDERIVATIVE,
		HDERIVATIVE,
		RDERIVATIVE,
		XDERIVATIVE,
		YDERIVATIVE,
		ZDERIVATIVE,
		NONE
	};

	enum class IPType {
		NONE,
		COLECOLE,
		PELTON
	};

private:
	size_t NumFrequencies;
public:
	std::vector<FrequencyNode> Frequency;		
private:
	std::vector<HankelTransforms> Hankel;
	
	double meanconductivity;
	double meanlog10conductivity;
	
	//Geometry Stuff
	double Xunrotated, Yunrotated, X, Y, Z, H;
	double R, R2, R3, R4, R5;
	double BigR, BigR2, BigR3, BigR5, BigR7;
	double XonR, YonR;
	cVec Source_Orientation;
	double xyrotation, cosxyrotation, sinxyrotation;
	cVec xaxis;
	cVec yaxis;
	
	//Hankle Stuff  
	size_t NumIntegrands;
	double LowerFractionalWidth;
	double UpperFractionalWidth;
	size_t number_integrand_calls;
	cdouble trapezoid_result[3];
	cdouble integrand_result[3];

public:	
	std::vector<LayerNode>  Layer;
	double ModellingLoopRadius = 0.0; //dipole by default
	size_t NumAbscissa;
	size_t NumLayers;
	RZeroMethod rzerotype;
	IPType iptype;
	size_t derivative_layer;
	CalculationType calculation_type;
	ResponseField Fields;

	cLEM(){ initialise(); }

	~cLEM(){};

	void initialise()
	{
		xaxis = cVec(1.0, 0.0, 0.0);
		yaxis = cVec(0.0, 1.0, 0.0);

		NumLayers = 0;
		NumFrequencies = 0;
		NumIntegrands = 3;
		NumAbscissa = 17;

		LowerFractionalWidth = 4.44;
		UpperFractionalWidth = 1.84;
	};

	void setlog10conductivitylog10thickness(const size_t nlayers, const double* log10conductivity, const double* log10thickness)
	{
		NumLayers = nlayers;
		Layer.resize(NumLayers);
		for (size_t i = 0; i < NumLayers; i++)Layer[i].Conductivity = pow(10.0, log10conductivity[i]);
		for (size_t i = 0; i < NumLayers - 1; i++)Layer[i].Thickness = pow(10.0, log10thickness[i]);
		setmeanconductivity();
		setmeanlog10conductivity();
	};

	void printearth()
	{
		size_t i;
		for (i = 0; i < NumLayers - 1; i++){
			printf("Layer %02zu:\t%10lf mS/m\t%10lf m\n", i, 1000.0*Layer[i].Conductivity, Layer[i].Thickness);
		}
		printf("Layer %02zu:\t%10lf mS/m\n", i, 1000.0*Layer[i].Conductivity);
	}

	std::vector<double> getconductivity()
	{
		std::vector<double> Conductivity;
		for (size_t i = 0; i < NumLayers; i++){
			Conductivity.push_back(Layer[i].Conductivity);
		}
		return Conductivity;
	}

	std::vector<double> getthickness()
	{
		std::vector<double> Thickness;
		for (size_t i = 0; i < NumLayers - 1; i++){
			Thickness.push_back(Layer[i].Thickness);
		}
		return Thickness;
	}

	void setmeanlog10conductivity()
	{
		if (NumLayers == 1){
			meanlog10conductivity = Layer[0].Conductivity;
			return;
		}
		//Returns the thickness weighted mean (in linear space) conductivity (but calculated in 10g10 space)
		double sumc = 0.0;
		double sumt = 0.0;
		for (size_t i = 0; i < NumLayers; i++){
			if (i < (NumLayers - 1)){
				sumc += log10(Layer[i].Conductivity)*Layer[i].Thickness;
				sumt += Layer[i].Thickness;
			}
			else{
				sumc += log10(Layer[i].Conductivity)*sumt; //make basement layer as thick as sum of all overlying
				sumt += sumt;
			}
		}
		meanlog10conductivity = pow(10.0, sumc / sumt);
	}

	void setmeanconductivity()
	{
		if (NumLayers == 1){
			meanconductivity = Layer[0].Conductivity;
			return;
		}
		//Returns the thickness weighted mean (in linear space) conductivity (but calculated in linear space)
		double sumc = 0.0;
		double sumt = 0.0;
		for (size_t i = 0; i < NumLayers; i++){
			if (i < (NumLayers - 1)){
				sumc += Layer[i].Conductivity*Layer[i].Thickness;
				sumt += Layer[i].Thickness;
			}
			else{
				sumc += Layer[i].Conductivity*sumt; //make basement layer as thick as sum of all overlying
				sumt += sumt;
			}
		}
		meanconductivity = sumc / sumt;
	}

	void setproperties(
		const std::vector<double>& conductivity,
		const std::vector<double>& thickness,
		const std::vector<double>& chargeability = std::vector<double>(0),
		const std::vector<double>& timeconstant = std::vector<double>(0),
		const std::vector<double>& frequencydependence = std::vector<double>(0))
	{
		NumLayers = conductivity.size();
		Layer.resize(NumLayers);
		for (size_t i = 0; i < NumLayers; i++){
			if (i < NumLayers - 1) Layer[i].Thickness = thickness[i];
			Layer[i].Conductivity = conductivity[i];
			if (chargeability.size() > 0){
				Layer[i].Chargeability = chargeability[i];
				Layer[i].TimeConstant = timeconstant[i];
				Layer[i].FrequencyDependence = frequencydependence[i];
			}
			else{
			  Layer[i].Chargeability = 0.0;
			  Layer[i].TimeConstant = 0.0;
			  Layer[i].FrequencyDependence = 0.0;
		   }
		}
		setmeanconductivity();
		setmeanlog10conductivity();
	};

	void setproperties(const cEarth1D& E)
	{
		setproperties(E.conductivity, E.thickness, E.chargeability, E.timeconstant, E.frequencydependence);
	}

	void setconductivitythickness(const std::vector<double>& conductivity, const std::vector<double>& thickness)
	{		
		setproperties(conductivity, thickness);		
	};

	void setconductivitythickness(const size_t nlayers, const double* conductivity, const double* thickness)
	{
		std::vector<double> c(conductivity, conductivity + nlayers);
		std::vector<double> t(thickness, thickness + nlayers - 1);
		setproperties(c, t);
	};

	cVec pitchrolldipole(double pitch, double roll)
	{
		//X = +ve in flight direction
		//Y = +ve on left wing
		//Z = +ve vertical up
		//ie different to Fugro convention

		//Steer left is positive yaw     X->Y axis
		//Left wing up is positive roll  Y->Z axis
		//Nose down is positive pitch	 Z->X axis

		cVec orientation = cVec(0.0, 0.0, 1.0);
		if (pitch != 0.0) orientation = orientation.rotate(pitch, yaxis);
		if (roll != 0.0) orientation = orientation.rotate(roll, xaxis);
		return orientation;
	};

	void setxyrotation()
	{
		//xyrotation is the anticlockwise angle (in degrees) 
		//that the horizontal coordinate system has to be rotated
		//so the horizontal dipole is all Y directed. 
		if (Source_Orientation.x == 0.0 && Source_Orientation.y == 0.0){
			xyrotation = 0.0;
		}
		else{
			xyrotation = atan2(Source_Orientation.y, Source_Orientation.x);
			xyrotation = xyrotation*R2D - 90.0;
		}
		cosxyrotation = cos(xyrotation*D2R);
		sinxyrotation = sin(xyrotation*D2R);
	}

	void xyrotate(const double& xin, const double& yin, double* xout, double* yout)
	{
		*xout =  xin*cosxyrotation + yin*sinxyrotation;
		*yout = -xin*sinxyrotation + yin*cosxyrotation;
	}

	void unxyrotate(const double& xin, const double& yin, double* xout, double* yout)
	{
		*xout = xin*cosxyrotation - yin*sinxyrotation;
		*yout = xin*sinxyrotation + yin*cosxyrotation;
	}

	void unxyrotateandscale(RealField& f, double scalefactor)
	{
		double x, y;
		unxyrotate(f.x, f.y, &x, &y);
		f.x = x*scalefactor;
		f.y = y*scalefactor;
		f.z = f.z*scalefactor;
	}

	void unxyrotateandscale(ComplexField& f, double scalefactor)
	{
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

	void setR(double r)
	{
		R = r;
		R2 = R*R;
		R3 = R2*R;
		R4 = R3*R;
		R5 = R4*R;
	}

	void setgeometry(cVec source_orientation, double h, double x, double y, double z)
	{
		Xunrotated = x;
		Yunrotated = y;
		Source_Orientation = source_orientation;
		setxyrotation();
		xyrotate(Xunrotated, Yunrotated, &X, &Y);
		R = sqrt(X*X + Y*Y);
		R2 = R*R;
		R3 = R2*R;
		R4 = R3*R;
		R5 = R4*R;
		Z = z;
		H = h;

		BigR = sqrt((Z - H)*(Z - H) + R*R);
		BigR2 = BigR*BigR;
		BigR3 = BigR*BigR2;
		BigR5 = BigR3*BigR2;
		BigR7 = BigR5*BigR2;
		
		if (R == 0.0){
			XonR = 0.0;
			YonR = 0.0;
		}
		else {
			XonR = X / R;
			YonR = Y / R;
		}
	}
		
	static cdouble ip_colecole_conductivity(
		const double& conductivity,
		const double& chargeability,
		const double& timeconstant,
		const double& frequencydependence,
		const double& omega)
	{
		//c = c0 - c0*(N / (1 + (1 - N)*(j*omega*T) ^ K));
		if (chargeability == 0.0){
			return cdouble(conductivity, 0.0);
		}
		else{
			cdouble c = conductivity - conductivity * (chargeability / (1.0 + (1.0 - chargeability)*(std::pow(cdouble(0.0, omega*timeconstant), frequencydependence))));
			return  c;
		}

	}

	static cdouble ip_pelton_conductivity(
		const double& conductivity,
		const double& chargeability,
		const double& timeconstant,
		const double& frequencydependence,
		const double& omega)
	{
		//p = p0[1 - m*(1 - (1 - 1/(1 + (j*omega*T) ^ K));
		if (chargeability == 0.0){
			return cdouble(conductivity, 0.0);
		}
		else{
			double  rho0 = 1.0 / conductivity;
			cdouble rho = rho0*(1.0 - chargeability*(1.0 - (1.0 / (1.0 + std::pow(cdouble(0.0, omega*timeconstant), frequencydependence)))));
			return  1.0/rho;
		}

	}

	void init_frequencies(const std::vector<double>& frequencies)
	{
		NumFrequencies = frequencies.size();
		if (Frequency.size() != NumFrequencies)Frequency.resize(NumFrequencies);
		if (Hankel.size() != NumFrequencies)Hankel.resize(NumFrequencies);
		for (size_t fi = 0; fi < NumFrequencies; fi++) {
			double omega = TWOPI * frequencies[fi];
			double muzeroomega = MUZERO * omega;
			Frequency[fi].Frequency = frequencies[fi];
			Frequency[fi].Omega = omega;
			Frequency[fi].MuZeroOmega = muzeroomega;
			Frequency[fi].iMuZeroOmega = cdouble(0.0, muzeroomega);
		}
	};

	void init_frequency(const size_t& fi)
	{
		init_integration_nodes(fi);		
		for (size_t ai = 0; ai < NumAbscissa; ai++) {
			init_abscissa(fi, ai);
		}		
	};

	void init_abscissa(const size_t& fi, const size_t& ai)
	{
		Frequency[fi].Abscissa[ai].Layer.resize(NumLayers);
		for (size_t li = 0; li < NumLayers; li++) {
			cdouble u;
			if (Layer[li].Chargeability == 0.0) {
				double gamma2 = Layer[li].Conductivity * Frequency[fi].MuZeroOmega;
				u = sqrt(cdouble(Frequency[fi].Abscissa[ai].Lambda2, gamma2));
			}
			else {
				cdouble cip;
				if (iptype == IPType::COLECOLE) {
					cip = ip_colecole_conductivity(Layer[li].Conductivity, Layer[li].Chargeability, Layer[li].TimeConstant, Layer[li].FrequencyDependence, Frequency[fi].Omega);
				}
				else {
					cip = ip_pelton_conductivity(Layer[li].Conductivity, Layer[li].Chargeability, Layer[li].TimeConstant, Layer[li].FrequencyDependence, Frequency[fi].Omega);
				}
				cdouble gamma2 = cip * cdouble(0.0, Frequency[fi].MuZeroOmega);
				u = std::sqrt(Frequency[fi].Abscissa[ai].Lambda2 + gamma2);
			}
			Frequency[fi].Abscissa[ai].Layer[li].U = u;
			if (li < NumLayers - 1)Frequency[fi].Abscissa[ai].Layer[li].Exp2UT = exp(-2.0 * u * Layer[li].Thickness);
		}
		init_layer_matrices(fi, ai);
		init_pmatrix(fi, ai);
	};

	double approximatehalfspace(const size_t& fi)
	{
		//approximate halfspace for the frequency at index fi
		if (NumLayers == 1)return Layer[0].Conductivity;

		//peak lambda
		double  peaklambda = sqrt(Frequency[fi].MuZeroOmega * meanlog10conductivity / 4.0);
		cdouble rz  = rzero_recursive(fi, peaklambda);
				
		cdouble  v = (1.0 + rz);

		v = Frequency[fi].iMuZeroOmega*v*v;
		return (-4.0*peaklambda*peaklambda*rz / v).real();
	}
		
	inline cdouble rzero_recursive(const size_t& fi, const double& lambda) const 
	{
		//Wait's recursive formulation
		const double lambda2 = lambda * lambda;
		const double muzeroomega = Frequency[fi].MuZeroOmega;
		const cdouble imuzeroomega(0.0, muzeroomega);

		double gamma2 = muzeroomega * Layer[NumLayers - 1].Conductivity;
		cdouble u = std::sqrt(cdouble(lambda2, gamma2));
		cdouble y = u / imuzeroomega;

		int i = (int)NumLayers - 2;
		while (i >= 0) {
			gamma2 = muzeroomega * Layer[i].Conductivity;
			u = std::sqrt(cdouble(lambda2, gamma2));
			const cdouble Nn = u / imuzeroomega;
			const cdouble v = u * Layer[i].Thickness;

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

	inline cdouble rzero_propogationmatrix(const size_t& fi, const size_t& ai)
	{
		//Oldenberg's propogation matrix formulation
		init_layer_matrices(fi, ai);
		init_pmatrix(fi, ai);
		return Frequency[fi].Abscissa[ai].P21onP11;
	};

	inline void init_layer_matrices(const size_t& fi, const size_t& ai)
	{
		cdouble e, eh, e1, e2;

		AbscissaNode& A = Frequency[fi].Abscissa[ai];

		//M1  
		e = A.Layer[0].U / A.Lambda;
		eh = e / 2.0;
		e1 = 0.5 + eh;
		e2 = 0.5 - eh;

		A.Layer[0].LayerMatrix.e11 = e1;
		A.Layer[0].LayerMatrix.e12 = e2;
		A.Layer[0].LayerMatrix.e21 = e2;
		A.Layer[0].LayerMatrix.e22 = e1;

		for (size_t li = 1; li < NumLayers; li++){
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
		if (NumLayers > 1){
			A.Layer[1].LayerPreMatrix = A.Layer[0].LayerMatrix;
			for (size_t li = 2; li < NumLayers; li++){
				A.Layer[li].LayerPreMatrix = A.Layer[li - 1].LayerPreMatrix * A.Layer[li - 1].LayerMatrix;
			}
		}

		//Set Postmatrices - postmatrix for last layer does not apply  
		if (NumLayers > 1){
			A.Layer[NumLayers - 2].LayerPostMatrix = A.Layer[NumLayers - 1].LayerMatrix;
			if (NumLayers > 2){
				for (size_t li = NumLayers - 2; li-- > 0;){
					A.Layer[li].LayerPostMatrix = A.Layer[li + 1].LayerMatrix * A.Layer[li + 1].LayerPostMatrix;
				}
			}
		}		
	};

	inline void init_pmatrix(const size_t& fi, const size_t& ai)
	{
		AbscissaNode& A = Frequency[fi].Abscissa[ai];
		//Set Full matrix
		if (NumLayers == 1){
			A.P_Full = A.Layer[0].LayerMatrix;
		}
		else{
			A.P_Full = Frequency[fi].Abscissa[ai].Layer[NumLayers - 1].LayerPreMatrix * Frequency[fi].Abscissa[ai].Layer[NumLayers - 1].LayerMatrix;
		}
		A.P21onP11 = A.P_Full.e21 / A.P_Full.e11;
	}

	inline PropogationMatrix dMjdCj(const size_t& fi, const size_t& ai, const size_t& li)
	{
		PropogationMatrix m;

		if (li == 0){
			cdouble a = Frequency[fi].iMuZeroOmega / (4.0*Frequency[fi].Abscissa[ai].Lambda*Frequency[fi].Abscissa[ai].Layer[li].U);
			m.e11 = a;
			m.e12 = -a;
			m.e21 = -a;
			m.e22 = a;
			return m;
		}
		else{
			cdouble a = Frequency[fi].iMuZeroOmega / (4.0*Frequency[fi].Abscissa[ai].Layer[li - 1].U * Frequency[fi].Abscissa[ai].Layer[li].U);
			cdouble ae = a*Frequency[fi].Abscissa[ai].Layer[li - 1].Exp2UT;
			m.e11 = a;
			m.e12 = -a;
			m.e21 = -ae;
			m.e22 = ae;
			return m;
		}

	}
	
	inline PropogationMatrix dMjplus1dCj(const size_t& fi, const size_t& ai, const size_t& li)
	{
		cdouble duds = Frequency[fi].iMuZeroOmega / (2.0*Frequency[fi].Abscissa[ai].Layer[li].U);
		cdouble y = Frequency[fi].Abscissa[ai].Layer[li + 1].U / Frequency[fi].Abscissa[ai].Layer[li].U;

		cdouble dydu = -y / Frequency[fi].Abscissa[ai].Layer[li].U;
		cdouble dyds = dydu*duds;

		cdouble v = Frequency[fi].Abscissa[ai].Layer[li].Exp2UT;
		cdouble dvdu = -2.0*Layer[li].Thickness * Frequency[fi].Abscissa[ai].Layer[li].Exp2UT;

		cdouble dvds = dvdu*duds;

		cdouble ydvds = y*dvds;
		cdouble vdyds = v*dyds;


		PropogationMatrix M;
		M.e11 = 0.5*dyds;
		M.e12 = -M.e11;
		M.e21 = 0.5*(dvds - (ydvds + vdyds));
		M.e22 = 0.5*(dvds + (ydvds + vdyds));

		return M;

	}
	
	inline PropogationMatrix dPdCj(const size_t& fi, const size_t& ai, const size_t& li)
	{
		AbscissaNode& A = Frequency[fi].Abscissa[ai];

		PropogationMatrix tmp;
		//One layer case
		if (NumLayers == 1)return tmp = dMjdCj(fi, ai, li);

		//Not last layer
		if (li < NumLayers - 1){
			PropogationMatrix M =
				(dMjdCj(fi, ai, li)      * A.Layer[li + 1].LayerMatrix)
				+ (A.Layer[li].LayerMatrix * dMjplus1dCj(fi, ai, li));

			//First layer case
			if (li == 0){
				if (NumLayers == 2) return M;
				return M * A.Layer[li + 1].LayerPostMatrix;
			}
			else if (li == NumLayers - 2){
				return A.Layer[li].LayerPreMatrix * M;
			}
			else{
				return A.Layer[li].LayerPreMatrix * M * A.Layer[li + 1].LayerPostMatrix;
			}
		}
		//Last layer
		else{
			return A.Layer[li].LayerPreMatrix * dMjdCj(fi, ai, li);
		}

	}
	
	inline cdouble dP21onP11dCj(const size_t& fi, const size_t& ai, const size_t& li)
	{
		PropogationMatrix m = dPdCj(fi, ai, li);
		return m.e21 / Frequency[fi].Abscissa[ai].P_Full.e11 - m.e11*Frequency[fi].Abscissa[ai].P21onP11 / Frequency[fi].Abscissa[ai].P_Full.e11;
	}
	
	inline PropogationMatrix dMjplus1dTj(const size_t& fi, const size_t& ai, const size_t& li)
	{
		cdouble y = Frequency[fi].Abscissa[ai].Layer[li + 1].U / Frequency[fi].Abscissa[ai].Layer[li].U;
		cdouble dvdt = -2.0*Frequency[fi].Abscissa[ai].Layer[li].U*Frequency[fi].Abscissa[ai].Layer[li].Exp2UT;

		PropogationMatrix m;
		m.e11 = 0.0;
		m.e12 = 0.0;
		m.e21 = 0.5*(1.0 - y)*dvdt;
		m.e22 = 0.5*(1.0 + y)*dvdt;
		return m;

	}
	
	inline PropogationMatrix dPdTj(const size_t& fi, const size_t& ai, const size_t& li)
	{
		AbscissaNode& A = Frequency[fi].Abscissa[ai];

		PropogationMatrix tmp;
		tmp.e11 = tmp.e12 = tmp.e21 = tmp.e22 = 0;

		//One layer case
		if (NumLayers == 1) return tmp;

		PropogationMatrix M = A.Layer[li].LayerMatrix * dMjplus1dTj(fi, ai, li);

		//First layer case
		if (li == 0){
			if (NumLayers == 2)return M;
			return M * A.Layer[li + 1].LayerPostMatrix;
		}
		else if (li > 0 && li < NumLayers - 1){
			if (li == NumLayers - 2) return A.Layer[li].LayerPreMatrix * M;
			return A.Layer[li].LayerPreMatrix * M * A.Layer[li + 1].LayerPostMatrix;
		}
		//Last layer case
		else{
			glog.warningmsg(_SRC_,"Zero thickness derivative for halfspace layer\n");
			tmp.e11 = tmp.e12 = tmp.e21 = tmp.e22 = 0;
			return tmp;
		}
	}

	inline cdouble dP21onP11dTj(const size_t& fi, const size_t& ai, const size_t& li)
	{
		PropogationMatrix m = dPdTj(fi, ai, li);
		return m.e21 / Frequency[fi].Abscissa[ai].P_Full.e11 - m.e11*Frequency[fi].Abscissa[ai].P21onP11 / Frequency[fi].Abscissa[ai].P_Full.e11;
	}

	void init_integration_nodes(const size_t& fi)
	{
		FrequencyNode& F = Frequency[fi];;
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
		for (size_t ai = 0; ai < NumAbscissa; ai++){
			AbscissaNode& A = F.Abscissa[ai];
			double lambda = exp(loglambda);
			A.Lambda = lambda;
			A.Lambda2 = A.Lambda * lambda;
			A.Lambda3 = A.Lambda2 * lambda;
			A.Lambda4 = A.Lambda3 * lambda;
			A.LambdaR = A.Lambda * R;
			A.j0LambdaR = std::cyl_bessel_j(0,A.LambdaR);
			A.j1LambdaR = std::cyl_bessel_j(1,A.LambdaR);
			loglambda += F.AbscissaSpacing;
		}
	}

	void dointegrals(const size_t& fi)
	{
		dointegrals_trapezoid(fi);		
	}

	inline void dointegrals_trapezoid(const size_t& fi)
	{
		HankelTransforms& H = Hankel[fi];

		number_integrand_calls = 0;

		trapezoid(fi);//the results go into the variable trapezoid_result

		if (calculation_type == CalculationType::FORWARDMODEL){
			H.I0.FM = trapezoid_result[0];
			H.I1.FM = trapezoid_result[1];
			H.I2.FM = trapezoid_result[2];
		}
		else if (calculation_type == CalculationType::CONDUCTIVITYDERIVATIVE){
			H.I0.dC = trapezoid_result[0];
			H.I1.dC = trapezoid_result[1];
			H.I2.dC = trapezoid_result[2];
		}
		else if (calculation_type == CalculationType::THICKNESSDERIVATIVE){
			H.I0.dT = trapezoid_result[0];
			H.I1.dT = trapezoid_result[1];
			H.I2.dT = trapezoid_result[2];
		}
		else if (calculation_type == CalculationType::ZDERIVATIVE){
			H.I0.dZ = trapezoid_result[0];
			H.I1.dZ = trapezoid_result[1];
			H.I2.dZ = trapezoid_result[2];
		}
		else if (calculation_type == CalculationType::HDERIVATIVE){
			H.I0.dH = trapezoid_result[0];
			H.I1.dH = trapezoid_result[1];
			H.I2.dH = trapezoid_result[2];
		}
		else if (calculation_type == CalculationType::RDERIVATIVE || calculation_type == CalculationType::XDERIVATIVE || calculation_type == CalculationType::YDERIVATIVE){
			H.I0.dR = trapezoid_result[0];
			H.I1.dR = trapezoid_result[1];
			H.I2.dR = trapezoid_result[2];
		}
		else
		{
			glog.errormsg(_SRC_,"LE::dointegrals_trapezoid Calculation type %lu not yet implemented\n", calculation_type);
		}
	}
	inline void trapezoid(const size_t& fi)
	{
		std::vector<cdouble> integrand1(3);
		std::vector<cdouble> integrand2(3);
		std::vector<cdouble> integrand3(3);

		trapezoid_result[0] = cdouble(0.0, 0.0);
		trapezoid_result[1] = cdouble(0.0, 0.0);
		trapezoid_result[2] = cdouble(0.0, 0.0);

		//First and last abscissa
		integrand(fi, 0);
		for (size_t ii = 0; ii < NumIntegrands; ii++){
			trapezoid_result[ii] += integrand_result[ii];
		}

		integrand(fi, NumAbscissa - 1);
		for (size_t ii = 0; ii < NumIntegrands; ii++){
			trapezoid_result[ii] += integrand_result[ii];
		}

		for (size_t ii = 0; ii < NumIntegrands; ii++){
			trapezoid_result[ii] *= 0.5;
		}

		//Cenral Abscissas
		for (size_t ai = 1; ai < NumAbscissa - 1; ai++){
			integrand(fi, ai);
			for (size_t ii = 0; ii < NumIntegrands; ii++){
				trapezoid_result[ii] += integrand_result[ii];
			}
		}

		for (size_t ii = 0; ii < NumIntegrands; ii++){
			trapezoid_result[ii] *= Frequency[fi].AbscissaSpacing;
		}
	}
	inline void integrand(const size_t& fi, const size_t& ai)
	{
		AbscissaNode& A = Frequency[fi].Abscissa[ai];
		number_integrand_calls++;

		double loopfactor = 1.0;
		if (ModellingLoopRadius > 0.0){
			double lambda_a = A.Lambda*ModellingLoopRadius;
			loopfactor = 2.0 *  std::cyl_bessel_j(1,lambda_a) / lambda_a;
		}

		double& lambdar = A.LambdaR;
		double& j0 = A.j0LambdaR;
		double& j1 = A.j1LambdaR;
		const double e = exp(-(Z + H)*A.Lambda);
		const double l2e = A.Lambda2*e;
		const double l3e = A.Lambda3*e;
		const double l4e = A.Lambda4*e;

		cdouble k;
		switch (calculation_type){
		case CalculationType::FORWARDMODEL:
			k = loopfactor * A.P21onP11;
			integrand_result[0] = k*l3e*j0;
			integrand_result[1] = k*l3e*j1;
			integrand_result[2] = k*l2e*j1;
			break;
		case CalculationType::CONDUCTIVITYDERIVATIVE:
			k = loopfactor * dP21onP11dCj(fi, ai, derivative_layer);
			integrand_result[0] = k*l3e*j0;
			integrand_result[1] = k*l3e*j1;
			integrand_result[2] = k*l2e*j1;
			break;
		case CalculationType::THICKNESSDERIVATIVE:
			k = loopfactor * dP21onP11dTj(fi, ai, derivative_layer);
			integrand_result[0] = k*l3e*j0;
			integrand_result[1] = k*l3e*j1;
			integrand_result[2] = k*l2e*j1;
			break;
		case CalculationType::ZDERIVATIVE:
			k = loopfactor * A.P21onP11;
			integrand_result[0] = k * -l4e*j0;
			integrand_result[1] = k * -l4e*j1;
			integrand_result[2] = k * -l3e*j1;
			break;
		case CalculationType::HDERIVATIVE:
			k = loopfactor * A.P21onP11;;
			integrand_result[0] = k * -l4e*j0;
			integrand_result[1] = k * -l4e*j1;
			integrand_result[2] = k * -l3e*j1;
			break;
		case CalculationType::RDERIVATIVE:
		case CalculationType::XDERIVATIVE:
		case CalculationType::YDERIVATIVE:
			k = loopfactor * A.P21onP11;
			if (R != 0.0) {
				integrand_result[0] = k * (-l4e*j1);			
				integrand_result[1] = k * (l4e*(j0 - j1 / lambdar));
				integrand_result[2] = k * (l3e*(j0 - j1 / lambdar));
			}
			else {
				integrand_result[0] = 0.0;
				integrand_result[1] = 0.0;
				integrand_result[2] = 0.0;
			}
			break;
		default:
			glog.errormsg(_SRC_,"LE::integrands Calculation type %lu not yet implemented", calculation_type);
			break;
		}

	}
	void setverticaldipolefields(const size_t& fi)
	{
		setverticaldipoleprimaryfields();
		setverticaldipolesecondaryfields(fi);
	}
	void setverticaldipoleprimaryfields()
	{
		Fields.v.p.x = 0.0; Fields.v.p.y = 0.0; Fields.v.p.z = 0.0;

		if (BigR == 0)return;
		if (Source_Orientation.z == 0.0)return;//ie no vertical dipole contribution

		if (calculation_type == CalculationType::FORWARDMODEL){
			Fields.v.p.x = THREEONFOURPI*X*(Z - H) / BigR5;
			Fields.v.p.y = THREEONFOURPI*Y*(Z - H) / BigR5;
			Fields.v.p.z = THREEONFOURPI*(Z - H)*(Z - H) / BigR5 - ONEONFOURPI / BigR3;
		}
		else if (calculation_type == CalculationType::CONDUCTIVITYDERIVATIVE || calculation_type == CalculationType::THICKNESSDERIVATIVE){
			Fields.v.p.x = 0.0;
			Fields.v.p.y = 0.0;
			Fields.v.p.z = 0.0;
		}
		else if (calculation_type == CalculationType::HDERIVATIVE){
			Fields.v.p.x = 0.0;
			Fields.v.p.y = 0.0;
			Fields.v.p.z = 0.0;
		}
		else if (calculation_type == CalculationType::ZDERIVATIVE){
			Fields.v.p.x = THREEONFOURPI*X*(1.0 / BigR5 - 5.0*(Z - H)*(Z - H) / BigR7);
			Fields.v.p.y = THREEONFOURPI*Y*(1.0 / BigR5 - 5.0*(Z - H)*(Z - H) / BigR7);
			Fields.v.p.z = THREEONFOURPI*(3.0*(Z - H) / BigR5 - 5.0*(Z - H)*(Z - H)*(Z - H) / BigR7);
		}
		else if (calculation_type == CalculationType::XDERIVATIVE || calculation_type == CalculationType::YDERIVATIVE || calculation_type == CalculationType::RDERIVATIVE){
			double dxdX = THREEONFOURPI*(Z - H)*(1.0 / BigR5 - 5.0*X*X / BigR7);
			double dydX = THREEONFOURPI*Y*(Z - H)*-5.0*X / BigR7;
			double dzdX = THREEONFOURPI*(Z - H)*(Z - H)*-5.0*X / BigR7 - ONEONFOURPI*-3.0*X / BigR5;

			double dxdY = THREEONFOURPI*X*(Z - H)*-5.0*Y / BigR7;
			double dydY = THREEONFOURPI*(Z - H)*(1.0 / BigR5 - 5.0*Y*Y / BigR7);
			double dzdY = THREEONFOURPI*(Z - H)*(Z - H)*-5.0*Y / BigR7 - ONEONFOURPI*-3.0*Y / BigR5;

			if (calculation_type == CalculationType::XDERIVATIVE){
				double dXdXo = cosxyrotation;
				double dYdXo = -sinxyrotation;
				Fields.v.p.x = dxdX*dXdXo + dxdY*dYdXo;
				Fields.v.p.y = dydX*dXdXo + dydY*dYdXo;
				Fields.v.p.z = dzdX*dXdXo + dzdY*dYdXo;
			}
			else if (calculation_type == CalculationType::YDERIVATIVE){
				double dXdYo = sinxyrotation;
				double dYdYo = cosxyrotation;
				Fields.v.p.x = dxdX*dXdYo + dxdY*dYdYo;
				Fields.v.p.y = dydX*dXdYo + dydY*dYdYo;
				Fields.v.p.z = dzdX*dXdYo + dzdY*dYdYo;
			}
			else if (calculation_type == CalculationType::RDERIVATIVE){
				Fields.v.p.x = dxdX*XonR + dxdY*YonR;
				Fields.v.p.y = dydX*XonR + dydY*YonR;
				Fields.v.p.z = dzdX*XonR + dzdY*YonR;
			}
		}
		else {
			glog.errormsg(_SRC_,"LE::setverticaldipoleprimaryfields Calculation type %lu not yet implemented", calculation_type);
		}
		unxyrotateandscale(Fields.v.p, Source_Orientation.z);
	}
	void setverticaldipolesecondaryfields(const size_t& fi)
	{
		Fields.v.s.x = cdouble(0.0, 0.0); Fields.v.s.y = cdouble(0.0, 0.0); Fields.v.s.z = cdouble(0.0, 0.0);

		if (Source_Orientation.z == 0.0)return;//ie no vertical dipole contribution

		if (calculation_type == CalculationType::FORWARDMODEL){
			Fields.v.s.x = -ONEONFOURPI * XonR * Hankel[fi].I1.FM;
			Fields.v.s.y = -ONEONFOURPI * YonR * Hankel[fi].I1.FM;
			Fields.v.s.z = -ONEONFOURPI * Hankel[fi].I0.FM;
		}
		else if (calculation_type == CalculationType::CONDUCTIVITYDERIVATIVE){
			Fields.v.s.x = -ONEONFOURPI * XonR * Hankel[fi].I1.dC;
			Fields.v.s.y = -ONEONFOURPI * YonR * Hankel[fi].I1.dC;
			Fields.v.s.z = -ONEONFOURPI * Hankel[fi].I0.dC;
		}
		else if (calculation_type == CalculationType::THICKNESSDERIVATIVE){
			Fields.v.s.x = -ONEONFOURPI * XonR * Hankel[fi].I1.dT;
			Fields.v.s.y = -ONEONFOURPI * YonR * Hankel[fi].I1.dT;
			Fields.v.s.z = -ONEONFOURPI * Hankel[fi].I0.dT;
		}
		else if (calculation_type == CalculationType::HDERIVATIVE){
			//these are negative of d/dz derivatives
			Fields.v.s.x = -ONEONFOURPI * XonR * Hankel[fi].I1.dH;
			Fields.v.s.y = -ONEONFOURPI * YonR * Hankel[fi].I1.dH;
			Fields.v.s.z = -ONEONFOURPI * Hankel[fi].I0.dH;
		}
		else if (calculation_type == CalculationType::ZDERIVATIVE){
			Fields.v.s.x = -ONEONFOURPI * XonR * Hankel[fi].I1.dZ;
			Fields.v.s.y = -ONEONFOURPI * YonR * Hankel[fi].I1.dZ;
			Fields.v.s.z = -ONEONFOURPI * Hankel[fi].I0.dZ;
		}
		else if (calculation_type == CalculationType::XDERIVATIVE || calculation_type == CalculationType::YDERIVATIVE || calculation_type == CalculationType::RDERIVATIVE) {

			cdouble dxdX = 0.0; cdouble dydX = 0.0; cdouble dzdX = 0.0;
			cdouble dxdY = 0.0; cdouble dydY = 0.0; cdouble dzdY = 0.0;

			if(R != 0.0){
				dxdX = -ONEONFOURPI * (Hankel[fi].I1.FM*(1.0 / R - X * X / R3) + XonR * Hankel[fi].I1.dR * XonR);
				dydX = -ONEONFOURPI * Y * (Hankel[fi].I1.FM*(-X / R3) + (1.0 / R) * Hankel[fi].I1.dR * XonR);
				dzdX = -ONEONFOURPI * Hankel[fi].I0.dR * XonR;

				dxdY = -ONEONFOURPI * X * (Hankel[fi].I1.FM*(-Y / R3) + (1.0 / R) * Hankel[fi].I1.dR * YonR);
				dydY = -ONEONFOURPI * (Hankel[fi].I1.FM*(1.0 / R - Y * Y / R3) + YonR * Hankel[fi].I1.dR * YonR);
				dzdY = -ONEONFOURPI * Hankel[fi].I0.dR * YonR;
			}

			if (calculation_type == CalculationType::XDERIVATIVE){
				double dXdXo = cosxyrotation;
				double dYdXo = -sinxyrotation;
				Fields.v.s.x = dxdX*dXdXo + dxdY*dYdXo;
				Fields.v.s.y = dydX*dXdXo + dydY*dYdXo;
				Fields.v.s.z = dzdX*dXdXo + dzdY*dYdXo;
			}
			else if (calculation_type == CalculationType::YDERIVATIVE){
				double dXdYo = sinxyrotation;
				double dYdYo = cosxyrotation;
				Fields.v.s.x = dxdX*dXdYo + dxdY*dYdYo;
				Fields.v.s.y = dydX*dXdYo + dydY*dYdYo;
				Fields.v.s.z = dzdX*dXdYo + dzdY*dYdYo;
			}
			else if (calculation_type == CalculationType::RDERIVATIVE){
				Fields.v.s.x = dxdX*XonR + dxdY*YonR;
				Fields.v.s.y = dydX*XonR + dydY*YonR;
				Fields.v.s.z = dzdX*XonR + dzdY*YonR;
			}
		}
		else {
			glog.errormsg(_SRC_,"LE::setverticaldipolesecondaryfields Calculation type %lu not yet implemented", calculation_type);
		}

		unxyrotateandscale(Fields.v.s, Source_Orientation.z);

	}
	void sethorizontaldipolefields(const size_t& fi){
		sethorizontaldipoleprimaryfields();
		sethorizontaldipolesecondaryfields(fi);
	}
	void sethorizontaldipoleprimaryfields(){

		Fields.h.p.x = 0.0; Fields.h.p.y = 0.0; Fields.h.p.z = 0.0;

		if (BigR == 0)return;
		if (Source_Orientation.x == 0.0 && Source_Orientation.y == 0.0)return;//ie. not horizontal dipole contribution

		if (calculation_type == CalculationType::FORWARDMODEL){
			Fields.h.p.x = THREEONFOURPI*X*Y / BigR5;
			Fields.h.p.y = THREEONFOURPI*Y*Y / BigR5 - ONEONFOURPI / BigR3;
			Fields.h.p.z = THREEONFOURPI*Y*(Z - H) / BigR5;
		}
		else if (calculation_type == CalculationType::CONDUCTIVITYDERIVATIVE || calculation_type == CalculationType::THICKNESSDERIVATIVE){
			Fields.h.p.x = 0.0;
			Fields.h.p.y = 0.0;
			Fields.h.p.z = 0.0;
		}
		else if (calculation_type == CalculationType::HDERIVATIVE){
			Fields.h.p.x = 0.0;
			Fields.h.p.y = 0.0;
			Fields.h.p.z = 0.0;
		}
		else if (calculation_type == CalculationType::ZDERIVATIVE){
			Fields.h.p.x = THREEONFOURPI*X*Y*(-5.0*(Z - H) / BigR7);
			Fields.h.p.y = THREEONFOURPI*Y*Y*(-5.0*(Z - H) / BigR7) - ONEONFOURPI*(-3.0*(Z - H) / BigR5);
			Fields.h.p.z = THREEONFOURPI*Y*(1.0 / BigR5 - 5.0*(Z - H)*(Z - H) / BigR7);
		}
		else if (calculation_type == CalculationType::XDERIVATIVE || calculation_type == CalculationType::YDERIVATIVE || calculation_type == CalculationType::RDERIVATIVE){
			double dxdX = THREEONFOURPI*Y*(1.0 / BigR5 - 5.0*X*X / BigR7);
			double dydX = THREEONFOURPI*Y*Y*-5.0*X / BigR7 + ONEONFOURPI*3.0*X / BigR5;
			double dzdX = THREEONFOURPI*Y*(Z - H)*-5.0*X / BigR7;

			double dxdY = THREEONFOURPI*X*(1.0 / BigR5 - 5.0*Y*Y / BigR7);
			double dydY = THREEONFOURPI*(3.0*Y / BigR5 - 5.0*Y*Y*Y / BigR7);
			double dzdY = THREEONFOURPI*(Z - H)*(1.0 / BigR5 - 5.0*Y*Y / BigR7);

			if (calculation_type == CalculationType::XDERIVATIVE){
				double dXdXo = cosxyrotation;
				double dYdXo = -sinxyrotation;
				Fields.h.p.x = dxdX*dXdXo + dxdY*dYdXo;
				Fields.h.p.y = dydX*dXdXo + dydY*dYdXo;
				Fields.h.p.z = dzdX*dXdXo + dzdY*dYdXo;
			}
			else if (calculation_type == CalculationType::YDERIVATIVE){
				double dXdYo = sinxyrotation;
				double dYdYo = cosxyrotation;
				Fields.h.p.x = dxdX*dXdYo + dxdY*dYdYo;
				Fields.h.p.y = dydX*dXdYo + dydY*dYdYo;
				Fields.h.p.z = dzdX*dXdYo + dzdY*dYdYo;
			}
			else if (calculation_type == CalculationType::RDERIVATIVE){
				Fields.h.p.x = dxdX*XonR + dxdY*YonR;
				Fields.h.p.y = dydX*XonR + dydY*YonR;
				Fields.h.p.z = dzdX*XonR + dzdY*YonR;
			}			
		}
		else {
			glog.errormsg(_SRC_,"LE::sethorizontaldipoleprimaryfields Calculation type %lu not yet implemented", calculation_type);
		}

		double scalefactor = sqrt(Source_Orientation.x*Source_Orientation.x + Source_Orientation.y*Source_Orientation.y);
		unxyrotateandscale(Fields.h.p, scalefactor);

	}
	void sethorizontaldipolesecondaryfields(const size_t& fi){

		Fields.h.s.x = cdouble(0.0, 0.0); Fields.h.s.y = cdouble(0.0, 0.0); Fields.h.s.z = cdouble(0.0, 0.0);

		if (R == 0)return;
		if (Source_Orientation.x == 0.0 && Source_Orientation.y == 0.0)return;//ie. not horizontal dipole contribution

		if (calculation_type == CalculationType::FORWARDMODEL){
			Fields.h.s.x = ONEONFOURPI * (X*Y) / (R2)* (2.0*Hankel[fi].I2.FM / R - Hankel[fi].I0.FM);
			Fields.h.s.y = ONEONFOURPI * ((Y*Y - X*X)*Hankel[fi].I2.FM / R3 - Y*Y*Hankel[fi].I0.FM / R2);
			Fields.h.s.z = ONEONFOURPI * Y/R * Hankel[fi].I1.FM;
		}
		else if (calculation_type == CalculationType::CONDUCTIVITYDERIVATIVE){
			Fields.h.s.x = ONEONFOURPI * (X*Y) / (R2)* (2.0*Hankel[fi].I2.dC / R - Hankel[fi].I0.dC);
			Fields.h.s.y = ONEONFOURPI * ((Y*Y - X*X)*Hankel[fi].I2.dC / R3 - Y*Y*Hankel[fi].I0.dC / R2);
			Fields.h.s.z = ONEONFOURPI * Y/R * Hankel[fi].I1.dC;
		}
		else if (calculation_type == CalculationType::THICKNESSDERIVATIVE){
			Fields.h.s.x = ONEONFOURPI * (X*Y) / (R2)* (2.0*Hankel[fi].I2.dT / R - Hankel[fi].I0.dT);
			Fields.h.s.y = ONEONFOURPI * ((Y*Y - X*X)*Hankel[fi].I2.dT / R3 - Y*Y*Hankel[fi].I0.dT / R2);
			Fields.h.s.z = ONEONFOURPI * Y/R * Hankel[fi].I1.dT;
		}
		else if (calculation_type == CalculationType::HDERIVATIVE){
			Fields.h.s.x = ONEONFOURPI * (X*Y) / (R2)* (2.0*Hankel[fi].I2.dH / R - Hankel[fi].I0.dH);
			Fields.h.s.y = ONEONFOURPI * ((Y*Y - X*X)*Hankel[fi].I2.dH / R3 - Y*Y*Hankel[fi].I0.dH / R2);
			Fields.h.s.z = ONEONFOURPI * Y/R * Hankel[fi].I1.dH;
		}
		else if (calculation_type == CalculationType::ZDERIVATIVE){
			Fields.h.s.x = ONEONFOURPI * (X*Y) / (R2)* (2.0*Hankel[fi].I2.dZ / R - Hankel[fi].I0.dZ);
			Fields.h.s.y = ONEONFOURPI * ((Y*Y - X*X)*Hankel[fi].I2.dZ / R3 - Y*Y*Hankel[fi].I0.dZ / R2);
			Fields.h.s.z = ONEONFOURPI * Y/R * Hankel[fi].I1.dZ;
		}
		else if (calculation_type == CalculationType::XDERIVATIVE || calculation_type == CalculationType::YDERIVATIVE || calculation_type == CalculationType::RDERIVATIVE){
			cdouble a, c, d, e, f, h;
			cdouble dadx, dcdx, dddx, dedx, dfdx, dhdx;
			cdouble dady, dcdy, dddy, dedy, dfdy, dhdy;

			cdouble b, dbdx, dbdy, dbdr;
			cdouble g, dgdx, dgdy;

			//x = a*b;
			//y = (d-c)*Hankel[fi].I2 - f*Hankel[fi].I0 = e*Hankel[fi].I2 - f*Hankel[fi].I0 = g;

			a = X*Y / R2;
			dadx = (R2*Y - X*Y*2.0*X) / R4;
			dady = (R2*X - X*Y*2.0*Y) / R4;

			b = 2.0*Hankel[fi].I2.FM / R - Hankel[fi].I0.FM;
			dbdr = 2.0*(R*Hankel[fi].I2.dR - Hankel[fi].I2.FM) / R2 - Hankel[fi].I0.dR;
			dbdx = dbdr*XonR;
			dbdy = dbdr*YonR;

			c = X*X / R3;
			dcdx = (R2*2.0*X - X*X*3.0*X) / R5;
			dcdy = X*X*-3.0*Y / R5;

			d = Y*Y / R3;
			dddx = Y*Y*-3.0*X / R5;
			dddy = (R2*2.0*Y - Y*Y*3.0*Y) / R5;

			e = d - c;
			dedx = dddx - dcdx;
			dedy = dddy - dcdy;

			f = Y*Y / R2;
			dfdx = Y*Y*-2.0*X / R4;
			dfdy = (R2*2.0*Y - Y*Y*2.0*Y) / R4;

			g = e*Hankel[fi].I2.FM - f*Hankel[fi].I0.FM;
			dgdx = e*Hankel[fi].I2.dR*XonR + Hankel[fi].I2.FM*dedx - (f*Hankel[fi].I0.dR*XonR + Hankel[fi].I0.FM*dfdx);
			dgdy = e*Hankel[fi].I2.dR*YonR + Hankel[fi].I2.FM*dedy - (f*Hankel[fi].I0.dR*YonR + Hankel[fi].I0.FM*dfdy);

			cdouble dxdX = ONEONFOURPI*(a*dbdx + b*dadx);
			cdouble dxdY = ONEONFOURPI*(a*dbdy + b*dady);

			cdouble dydX = ONEONFOURPI*(dgdx);
			cdouble dydY = ONEONFOURPI*(dgdy);

			h = Y / R;
			dhdx = -X*Y / R3;
			dhdy = 1.0 / R - Y*Y / R3;
			cdouble dzdX = ONEONFOURPI*(Hankel[fi].I1.FM*dhdx + h*Hankel[fi].I1.dR*XonR);
			cdouble dzdY = ONEONFOURPI*(Hankel[fi].I1.FM*dhdy + h*Hankel[fi].I1.dR*YonR);

			if (calculation_type == CalculationType::XDERIVATIVE){
				double dXdXo = cosxyrotation;
				double dYdXo = -sinxyrotation;
				Fields.h.s.x = dxdX*dXdXo + dxdY*dYdXo;
				Fields.h.s.y = dydX*dXdXo + dydY*dYdXo;
				Fields.h.s.z = dzdX*dXdXo + dzdY*dYdXo;
			}
			else if (calculation_type == CalculationType::YDERIVATIVE){
				double dXdYo = sinxyrotation;
				double dYdYo = cosxyrotation;
				Fields.h.s.x = dxdX*dXdYo + dxdY*dYdYo;
				Fields.h.s.y = dydX*dXdYo + dydY*dYdYo;
				Fields.h.s.z = dzdX*dXdYo + dzdY*dYdYo;
			}
			else if (calculation_type == CalculationType::RDERIVATIVE){
				Fields.h.s.x = dxdX*XonR + dxdY*YonR;
				Fields.h.s.y = dydX*XonR + dydY*YonR;
				Fields.h.s.z = dzdX*XonR + dzdY*YonR;
			}
		}
		else {
			glog.errormsg(_SRC_,"LE::sethorizontaldipolesecondaryfields Calculation type %lu not yet implemented", calculation_type);
		}

		double scalefactor = sqrt(Source_Orientation.x*Source_Orientation.x + Source_Orientation.y*Source_Orientation.y);
		unxyrotateandscale(Fields.h.s, scalefactor);
	}
	void setprimaryfields()
	{
		sethorizontaldipoleprimaryfields();
		setverticaldipoleprimaryfields();
		Fields.t.p.x = Fields.v.p.x + Fields.h.p.x;
		Fields.t.p.y = Fields.v.p.y + Fields.h.p.y;
		Fields.t.p.z = Fields.v.p.z + Fields.h.p.z;
	}
	void setsecondaryfields(const size_t& fi)
	{
		sethorizontaldipolesecondaryfields(fi);
		setverticaldipolesecondaryfields(fi);
		Fields.t.s.x = Fields.v.s.x + Fields.h.s.x;
		Fields.t.s.y = Fields.v.s.y + Fields.h.s.y;
		Fields.t.s.z = Fields.v.s.z + Fields.h.s.z;
	}

	//Horizontal coplanar
	cdouble ppmHCP(const size_t& fi){
		return 1.0e6*(R3*Hankel[fi].I0.FM);
	}
	cdouble dppmHCPdC(const size_t& fi){
		return 1.0e6*(R3*Hankel[fi].I0.dC);
	}
	cdouble dppmHCPdT(const size_t& fi){
		return 1.0e6*(R3*Hankel[fi].I0.dT);
	}
	cdouble dppmHCPdZ(const size_t& fi){
		return 1.0e6*(R3*Hankel[fi].I0.dZ);
	}
	cdouble dppmHCPdH(const size_t& fi){
		return 1.0e6*(R3*Hankel[fi].I0.dH);
	}
	cdouble dppmHCPdR(const size_t& fi){
		return 1.0e6*(R3*Hankel[fi].I0.dR + 3.0*R2*Hankel[fi].I0.FM);
	}

	//Perpendicular
	cdouble ppmPER(const size_t& fi){
		return -1.0e6*(1.0 - 0.5*R3*Hankel[fi].I1.FM);
	}
	cdouble dppmPERdC(const size_t& fi){
		return 0.5e6*(R3*Hankel[fi].I1.dC);
	}
	cdouble dppmPERdT(const size_t& fi){
		return 0.5e6*(R3*Hankel[fi].I1.dT);
	}
	cdouble dppmPERdZ(const size_t& fi){
		return 0.5e6*(R3*Hankel[fi].I1.dZ);
	}
	cdouble dppmPERdH(const size_t& fi){
		return 0.5e6*(R3*Hankel[fi].I1.dH);
	}
	cdouble dppmPERdR(const size_t& fi){
		return 0.5e6*(R3*Hankel[fi].I1.dR + 3.0*R2*Hankel[fi].I1.FM);
	}

	//Vertical coaxial
	cdouble ppmVCX(const size_t& fi){
		return -0.5e6*(R2*Hankel[fi].I2.FM - R3*Hankel[fi].I0.FM);
	}
	cdouble dppmVCXdC(const size_t& fi){
		return -0.5e6*(R2*Hankel[fi].I2.dC - R3*Hankel[fi].I0.dC);
	}
	cdouble dppmVCXdT(const size_t& fi){
		return -0.5e6*(R2*Hankel[fi].I2.dT - R3*Hankel[fi].I0.dT);
	}
	cdouble dppmVCXdZ(const size_t& fi){
		return -0.5e6*(R2*Hankel[fi].I2.dZ - R3*Hankel[fi].I0.dZ);
	}
	cdouble dppmVCXdH(const size_t& fi){
		return -0.5e6*(R2*Hankel[fi].I2.dH - R3*Hankel[fi].I0.dH);
	}
	cdouble dppmVCXdR(const size_t& fi){
		return -0.5e6*(R2*Hankel[fi].I2.dR + 2.0*R*Hankel[fi].I2.FM - R3*Hankel[fi].I0.dR - 3.0*R2*Hankel[fi].I0.FM);
	}

	//Vertical coplanar
	cdouble ppmVCP(const size_t& fi){
		return 1.0e6*(R2*Hankel[fi].I2.FM);
	}
	cdouble dppmVCPdC(const size_t& fi){
		return 1.0e6*(R2*Hankel[fi].I2.dC);
	}
	cdouble dppmVCPdT(const size_t& fi){
		return 1.0e6*(R2*Hankel[fi].I2.dT);
	}
	cdouble dppmVCPdZ(const size_t& fi){
		return 1.0e6*(R2*Hankel[fi].I2.dZ);
	}
	cdouble dppmVCPdH(const size_t& fi){
		return 1.0e6*(R2*Hankel[fi].I2.dH);
	}
	cdouble dppmVCPdR(const size_t& fi){
		return 1.0e6*(R2*Hankel[fi].I2.dR + 2.0*R*Hankel[fi].I2.FM);
	}

};

#endif
