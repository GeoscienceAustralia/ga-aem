/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

//---------------------------------------------------------------------------
#ifndef _le_H
#define _le_H

#include <complex>
#include <vector>

#include "geometry3d.h"

typedef std::complex<double> cdouble;

//Formulation mainly from the book 
//Geo-Electromagnetism, Wait James, R. Academic Press 1982
  
struct HankelTransform{ 
  cdouble FM;  
  cdouble dC;  
  cdouble dT;
  cdouble dZ;
  cdouble dH;
  cdouble dR;    
};

struct HankelTransforms{  
  HankelTransform I0;
  HankelTransform I1;
  HankelTransform I2;  
};

struct RealField{  
  double x;
  double y;
  double z;  
};

struct ComplexField{  
  cdouble x;
  cdouble y;
  cdouble z;  
};

struct TotalField{
  RealField    p;//primary filed component
  ComplexField s;//secondary field component  
};
 
struct ResponseField{
  TotalField   v;//due to vertical dipole
  TotalField   h;//due to y directed horizontal dipole
  TotalField   t;//due to total dipole
};

struct PropogationMatrix{

  cdouble e11;
  cdouble e12;
  cdouble e21;
  cdouble e22;

  PropogationMatrix operator+(const PropogationMatrix& b)
  {
	  PropogationMatrix c=b;
	  c.e11 += e11;
	  c.e12 += e12;
	  c.e21 += e21;
	  c.e22 += e22;
	  return c;
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
  cdouble   U;
  cdouble   Exp2UT;
  PropogationMatrix LayerMatrix;
  PropogationMatrix LayerPreMatrix;
  PropogationMatrix LayerPostMatrix;
};

struct AbscissaNode{  
  double Lambda;
  double Lambda2;
  double Lambda3;
  double Lambda4;  
  double LambdaR;
  double j0LambdaR;
  double j1LambdaR;
  std::vector<AbscissaLayerNode> Layer;
  PropogationMatrix P_Full;
  cdouble P21onP11;  
};

struct FrequencyNode{
  double Frequency;
  double Omega;
  double MuZeroOmega; 
  cdouble iMuZeroOmega;
  double ApproximateHalfspace;
  double PeakLambda;
  double LowerBound;
  double UpperBound;
  double AbscissaSpacing;
  std::vector<AbscissaNode> Abscissa;    
};

struct LayerNode{
  double Conductivity;
  double Thickness;  
};

enum eRZeroMethod     { RZM_PROPOGATIONMATRIX, RZM_RECURSIVE };
enum eCalculationType { CT_FORWARDMODEL, CT_CONDUCTIVITYDERIVATIVE, CT_THICKNESSDERIVATIVE, CT_HDERIVATIVE, CT_RDERIVATIVE, CT_XDERIVATIVE,	CT_YDERIVATIVE, CT_ZDERIVATIVE };

class LE{

  public:
  
  double ModellingLoopRadius = 0.0;//dipole by default
  size_t NumFrequencies;
  std::vector<FrequencyNode> Frequency;  
  size_t NumAbscissa;
  size_t NumLayers;    
  std::vector<LayerNode>  Layer;  
  std::vector<HankelTransforms> Hankel;
  ResponseField Fields;

  ///////////////////////////////////////////
  //Initialisation
  LE();
  LE(size_t nlayers, double* conductivity, double* thickness);
  ~LE();
  

  /////////////////////////////////////////////
  //Utility
  void printearth();
  void initialise();
  void setconductivitythickness(const size_t nlayers, const double* conductivity, const double* thickness);
  void setconductivitythickness(const std::vector<double>& conductivity, const std::vector<double>& thickness);  
  void setlog10conductivitylog10thickness(const size_t nlayers, const double* log10conductivity, const double* log10thickness);

  void setfrequencies(const std::vector<double>& frequencies);    
  double meanconductivity;
  double meanlog10conductivity;
  void setmeanconductivity();
  void setmeanlog10conductivity();

  std::vector<double> getconductivity();
  std::vector<double> getthickness();

  ////////////////////////////////  
  //Geometry Stuff
  double Xunrotated,Yunrotated,X,Y,Z,H;
  double R,R2,R3,R4,R5;
  double BigR,BigR2,BigR3,BigR5,BigR7;
  double XonR,YonR;
  cVec Source_Orientation;     
  cVec pitchrolldipole(double pitch, double roll);
  void setR(double r);
  void setgeometry(cVec source_orientation, double h, double x, double y, double z);
  double xyrotation,cosxyrotation,sinxyrotation;  
  void setxyrotation();
  void xyrotate(double xin, double yin, double* xout, double* yout);
  void unxyrotate(double xin, double yin, double* xout, double* yout);
  void unxyrotateandscale(RealField* f, double scalefactor);
  void unxyrotateandscale(ComplexField* f, double scalefactor);
  cVec xaxis;
  cVec yaxis;  
    
  /////////////////////////////////////////////////////
  //Hankle Stuff  
  size_t NumIntegrands;  
  double LowerFractionalWidth;
  double UpperFractionalWidth;  
  size_t number_integrand_calls;
  void setintegrationnodes(const size_t& fi);
  void dointegrals(const size_t& fi);  
  void dointegrals_trapezoid(const size_t& fi);  
    
  void trapezoid(const size_t& fi);  
  cdouble trapezoid_result[3];

  void integrand(const size_t& fi, const size_t& ai);        
  cdouble integrand_result[3];
  
      
  //////////////////////////////////////////////////
  //Computations  
  size_t derivative_layer;
  eCalculationType calculation_type;  
  eRZeroMethod rzerotype;
  double approximatehalfspace(const size_t& fi);    
  void setfrequencyabscissalayers(const size_t& fi);
  void setlayermatrices(const size_t& fi, const size_t& ai);
  void setpmatrix(const size_t& fi, const size_t& ai);
  cdouble rzero(const size_t& fi, const double& lambda);
  cdouble rzero_recursive(const size_t& fi, const double& lambda);
  cdouble rzero_propogationmatrix(const size_t& fi, const size_t& ai);
    
  //////////////////////////////////////////////////    
  //Propogtaion Matrix Stuff      
  PropogationMatrix dMjdCj(const size_t& fi, const size_t& ai, const size_t& li);
  PropogationMatrix dMjplus1dCj(const size_t& fi, const size_t& ai, const size_t& li);
  PropogationMatrix dPdCj(const size_t& fi, const size_t& ai, const size_t& li);
  PropogationMatrix dMjplus1dTj(const size_t& fi, const size_t& ai, const size_t& li);
  PropogationMatrix dPdTj(const size_t& fi, const size_t& ai, const size_t& li);
  cdouble dP21onP11dCj(const size_t& fi, const size_t& ai, const size_t& li);
  cdouble dP21onP11dTj(const size_t& fi, const size_t& ai, const size_t& li);

  ///////////////////////////////////////////////
  //Field Responses  
  void setverticaldipoleprimaryfields();
  void setverticaldipolesecondaryfields(const size_t& fi);  
  void sethorizontaldipoleprimaryfields();
  void sethorizontaldipolesecondaryfields(const size_t& fi);    

  void sethorizontaldipolefields(const size_t& fi);
  void setverticaldipolefields(const size_t& fi);
  void setprimaryfields();
  void setsecondaryfields(const size_t& fi);

  ///////////////////////////////////////////////
  //Frequency Domain Responses
  cdouble ppmHCP(const size_t& fi);
  cdouble dppmHCPdC(const size_t& fi);
  cdouble dppmHCPdT(const size_t& fi);
  cdouble dppmHCPdZ(const size_t& fi);
  cdouble dppmHCPdH(const size_t& fi);
  cdouble dppmHCPdR(const size_t& fi);

  cdouble ppmPER(const size_t& fi);
  cdouble dppmPERdC(const size_t& fi);
  cdouble dppmPERdT(const size_t& fi);
  cdouble dppmPERdZ(const size_t& fi);
  cdouble dppmPERdH(const size_t& fi);
  cdouble dppmPERdR(const size_t& fi);

  cdouble ppmVCX(const size_t& fi);
  cdouble dppmVCXdC(const size_t& fi);
  cdouble dppmVCXdT(const size_t& fi);
  cdouble dppmVCXdZ(const size_t& fi);
  cdouble dppmVCXdH(const size_t& fi);
  cdouble dppmVCXdR(const size_t& fi);

  cdouble ppmVCP(const size_t& fi);
  cdouble dppmVCPdC(const size_t& fi);
  cdouble dppmVCPdT(const size_t& fi);
  cdouble dppmVCPdZ(const size_t& fi);
  cdouble dppmVCPdH(const size_t& fi);
  cdouble dppmVCPdR(const size_t& fi);
  
};
////////////////////////////////////////////////////////////////////////////

#endif
