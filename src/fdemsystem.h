/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/


#ifndef _fdemsystemclass_H
#define _fdemsystemclass_H

#include <cstring>
#include <vector>
#include "blocklanguage.h"
#include "earth1d.h"
#include "layeredearthmodeller.h"

struct sLoop{
	Vector3<double>    npos;
	Vector3<double>    pos;
	Vector3<double>    naxis;     
	Vector3<double>    axis;     
	Vector3<double>    rpy;	
	Matrix33<double>   IR;
};


struct sFDEmCoilSet{
  LayeredEarthModeller L;
  double Frequency;
  double Separation;
  COILSETTYPE Orientation;  
  
  sLoop tx;
  sLoop rx;    
  Vector3<double> vtxrx;  
};

struct sFDEmGeometry{  
	double birdheight;
	double birdroll;
	double birdpitch;
	double birdyaw;
};

struct sFDEmData{	
	std::vector<double> inphase;
	std::vector<double> quadrature;	
};

struct sFDEmInversionOptions{		
	bool solve_conductivity;
	bool solve_thickness;
	bool solve_birdheight;	
	double lambda;
	double alpha;
	double minimumphid;
	double minimumimprovement;
	int maxiterations;	
};

struct sFDEmOutputOptions{		
	bool predictedbirdheight;

	bool conductivity;
	bool thickness;
	bool depth;	

	bool observeddata;	
	bool predicteddata;	
};

struct sFDEmResponse{  
	int     NC; //number of coilsets
	double* SIP;//secondary inphase field
	double* SQ;	//secondary quadrature field

	double* EIP;//std dev inphase    noise estimates	
	double* EQ; //std dev quadrature noise estimates	
};

enum AEMNOISETYPE{MAXPERCENTFLOOR,COMBINEPERCENTFLOOR};

struct sFDEmNoiseModel{  
  AEMNOISETYPE  type;
  cvector floor;
  cvector percentage;  
};





//---------------------------------------------------------------------------
class cFDEmSystem{

  public:
    
  std::string SystemName;   
  cBlock STM;
  void readsystemdescriptorfile(std::string systemdescriptorfile);        
  std::vector<sFDEmCoilSet> CoilSets;
  
  size_t NumberOfCoilSets(){ return CoilSets.size(); }
  
  void initialisecoilset(size_t cndex, double frequency, COILSETTYPE orientation, double separation);
  void setgeometry(const sFDEmGeometry& g);
  void setheight(const double height);
  void setrollpitchyaw(const double roll, const double pitch, const double yaw);

  //void setearth(const size_t nlayers, const double* conductivity, const double* thickness);
  void setearth(const std::vector<double>& conductivity, const std::vector<double>& thickness);
  void setearth(const cEarth1D& e);

  void setupcomputations();  
  std::vector<double> p();
  cvector s();  
  cvector ppms();
  cvector dppms(const CALCULATIONTYPE& calculationtype, const size_t& derivativelayer);  
  cvector noiseestimates(const cvector& response, const sFDEmNoiseModel& noisemodel);
  std::vector<double> cv2dv(const cvector& cv);
  cvector dv2cv(std::vector<double>& dv);
    
  cFDEmSystem();
  ~cFDEmSystem();  
};
////////////////////////////////////////////////////////////////////////////
#endif



