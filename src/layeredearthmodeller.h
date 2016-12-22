/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#ifndef layeredearthmodeller_H
#define layeredearthmodeller_H

#include <complex>
#include <vector>
#include "matrixvector23.h"

//Formulation mainly from the book 
//Geo-Electromagnetism, Wait James, R. Academic Press 1982
  
enum CALCULATIONTYPE
{
	eFM,eDX,eDY,eDZ,eDH,eDC,eDT,eDB
};

enum COILSETTYPE
{
	eHCP,eVCX,eVCP,ePER	
};

struct sHankelTransform{ 
  cdouble FM;  
  cdouble dC;    
  cdouble dT;
  cdouble dX;    
  cdouble dY;      
  cdouble dZ;
  cdouble dH;  
};


struct sPropogationMatrix{
  cdouble e11;
  cdouble e12;
  cdouble e21;
  cdouble e22;
};

struct sAbscissaLayerNode{      
  cdouble   U;
  cdouble   Exp2UT;
  sPropogationMatrix LayerMatrix;
  sPropogationMatrix LayerPreMatrix;
  sPropogationMatrix LayerPostMatrix;
};


struct sAbscissaNode{  
  double Lambda;
  double Lambda2;
  double Lambda3;
  double Lambda4;  
  double Lambda_r;
  double j0Lambda_r;
  double j1Lambda_r;  
  std::vector<sAbscissaLayerNode> Layer;
  sPropogationMatrix P_Full;
  cdouble P21onP11;
};

struct sLayerNode{
  double Conductivity;
  double Thickness;  
};

class LayeredEarthModeller{

  private:

  bool geometrychanged;
  bool earthchanged;

  size_t NumAbscissa;
  std::vector<sAbscissaNode> Abscissa;
  
  size_t NumLayers;    
  std::vector<double> Conductivity;
  std::vector<double> Thickness;  
  double  X,Y,Z,H; 

  double X2,Y2,X4,Y4,ZH,ZH2;
  double r,r2,r3,r5,r4,R,R2,R5,R7;
  

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
    
  ///////////////////////////////////////////
  //Initialisation
  LayeredEarthModeller();  
  ~LayeredEarthModeller();
  

  /////////////////////////////////////////////
  //Utility
  void printearth();
  void initialise();
  void setnumlayers(size_t nlayers);
  void setconductivitythickness(const size_t nlayers, const double* conductivity, const double* thickness);
  void setconductivitythickness(const std::vector<double>& conductivity, const std::vector<double>& thickness);  
  void setfrequency(const double frequency);      
  void setmeanconductivity();
  void setmeanlog10conductivity();
  ////////////////////////////////  
  //Geometry Stuff  
  void setx(double x);     
  void sety(double y);     
  void setz(double z);     
  void seth(double h);     
  void setr();     
  void setR();     
  void setgeometry(double h, double x, double y, double z);        
  /////////////////////////////////////////////////////     
  void setupcomputations();
  void dointegrals();        
  void trapezoid();  
  void compute_integrand(size_t ai);              
  //////////////////////////////////////////////////
  //Computations  
  CALCULATIONTYPE CalculationType;
  size_t DerivativeLayer;
  int  rzerotype;
  double approximatehalfspace();    
  cdouble rzero(double lambda);
  cdouble rzero_recursive(double lambda);
  cdouble rzero_propogationmatrix(size_t ai);

  void setnumabscissa(size_t nabscissa);
  void setintegrationnodes();
  void setfrequencyabscissalayers();
  void setlayermatrices(size_t ai);
  void setpmatrix(size_t ai);          
  //////////////////////////////////////////////////    
  //Propogtaion Matrix Stuff  
  //PropogationMatrix multiplymatrices(PropogationMatrix& m1, PropogationMatrix& m2);
  void multiplymatrices(const sPropogationMatrix& m1, const sPropogationMatrix& m2, sPropogationMatrix& m);
  void addmatrices(const sPropogationMatrix& m1, const sPropogationMatrix& m2, sPropogationMatrix& m);      
  sPropogationMatrix dMjdCj(size_t ai, size_t li);  
  sPropogationMatrix dMjplus1dCj(size_t ai, size_t li);    
  sPropogationMatrix dPdCj(size_t ai, size_t li);      
  sPropogationMatrix dMjplus1dTj(size_t ai, size_t li);
  sPropogationMatrix dPdTj(size_t ai, size_t li);  
  cdouble dP21onP11dCj(size_t ai, size_t li);  
  cdouble dP21onP11dTj(size_t ai, size_t li);  

  /////////////////////////////////////////////////////////
  void PTFM(Matrix33<double>& T);
  void dPTdX(Matrix33<double>& T);
  void dPTdY(Matrix33<double>& T);
  void dPTdZ(Matrix33<double>& T);
  void dPTdH(Matrix33<double>& T);  

  void STFM(Matrix33<cdouble>& T);
  void dSTdC(Matrix33<cdouble>& T);
  void dSTdT(Matrix33<cdouble>& T);
  void dSTdX(Matrix33<cdouble>& T);
  void dSTdY(Matrix33<cdouble>& T);
  void dSTdZ(Matrix33<cdouble>& T);
  void dSTdH(Matrix33<cdouble>& T);  


  void PT(const CALCULATIONTYPE& calculationtype, Matrix33<double>& T);
  void ST(const CALCULATIONTYPE& calculationtype, const size_t& derivativelayer, Matrix33<cdouble>& T);

  Vector3<double> pfield_if(const CALCULATIONTYPE& calculationtype, const Vector3<double>& sourcedirection);
  Vector3<cdouble> sfield_if(const CALCULATIONTYPE& calculationtype, const size_t& derivativelayer, const Vector3<double>& sourcedirection);
  
  double  p(const Vector3<double>& txdir, const Vector3<double>& rxdir);
  cdouble s(const Vector3<double>& txdir, const Vector3<double>& rxdir);

  double  dp(const CALCULATIONTYPE& calculationtype, const Vector3<double>& txdir,  const Vector3<double>& rxdir);
  cdouble ds(const CALCULATIONTYPE& calculationtype, const size_t& derivativelayer, const Vector3<double>& txdir,  const Vector3<double>& rxdir);
  
  cdouble dsdx(const Vector3<double>& txdir, const Vector3<double>& rxdir);
  cdouble dsdy(const Vector3<double>& txdir, const Vector3<double>& rxdir);
  cdouble dsdz(const Vector3<double>& txdir, const Vector3<double>& rxdir);
  cdouble dsdh(const Vector3<double>& txdir, const Vector3<double>& rxdir);
  cdouble dsdc(const size_t dlayer, const Vector3<double>& txdir, const Vector3<double>& rxdir);
  cdouble dsdt(const size_t dlayer, const Vector3<double>& txdir, const Vector3<double>& rxdir);  
  
  cdouble ppm(const Vector3<double>& txdir, const Vector3<double>& rxdir);  
  cdouble dppm(const CALCULATIONTYPE& calculationtype, const size_t& derivativelayer, const Vector3<double>& txdir, const Vector3<double>& rxdir);
};
////////////////////////////////////////////////////////////////////////////

#endif
