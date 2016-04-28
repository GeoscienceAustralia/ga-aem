/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

//---------------------------------------------------------------------------
#ifndef _tdemsystem_H
#define _tdemsystem_H

#include <complex>
#include "fftw3.h"
#include "general_utils.h"
#include "geometry3d.h"
#include "blocklanguage.h"
#include "le.h"

#define WAVEFORMTYPE_CURRENT  0  
#define WAVEFORMTYPE_RECEIVER 1  

struct sTDEmNoiseModelComponent{	
	double MultiplicativeNoise;
	std::vector<double> AdditiveNoise;	
};

struct sTDEmNoiseModel{	
	sTDEmNoiseModelComponent xcomponent;
	sTDEmNoiseModelComponent ycomponent;
	sTDEmNoiseModelComponent zcomponent;
};

struct sTDEmComponent{  
	double Primary;
	std::vector<double> Secondary;
};

struct sTDEmData{	
	sTDEmComponent xcomponent;
	sTDEmComponent ycomponent;
	sTDEmComponent zcomponent;
};

#ifndef _sTDEmGeometry_
#define _sTDEmGeometry_
struct sTDEmGeometry{  
	double tx_height;	
	double tx_roll;
	double tx_pitch;
	double tx_yaw;
	double txrx_dx;
	double txrx_dy;
	double txrx_dz;			
	double rx_roll;
	double rx_pitch;
	double rx_yaw;
};
#endif

#ifndef _sTDEmResponse_
#define _sTDEmResponse_
struct sTDEmResponse{  	
	double  PX;	
	double  PY;	
	double  PZ;	
	std::vector<double> SX;
	std::vector<double> SY;
	std::vector<double> SZ;
};
#endif

struct WindowSpecification{
	size_t SampleLow;
	size_t SampleHigh;
	size_t NumberOfSamples;
    
	double TimeLow;
	double TimeHigh;
	double TimeWidth;	

	std::vector<size_t> Sample;
	std::vector<double> Weight;
};


//---------------------------------------------------------------------------
class cTDEmSystem{

  public:


  cTDEmSystem();
  cTDEmSystem(std::string systemdescriptorfile);
  ~cTDEmSystem();
  void initialise();

  std::string SystemName;
  std::string SystemType;  
  LE Earth;
  cBlock STM;
  bool SaveDiagnosticFiles;
      
  cVec xaxis;
  cVec yaxis;
  cVec zaxis;

  std::string OutputType;
  std::string Normalisation;
  
  size_t NumberOfWindows;
  std::string WindowWeightingScheme;

  double BaseFrequency;
  double BasePeriod;
  double SampleFrequency;  
  double SampleInterval;  
  size_t SamplesPerWaveform;
  size_t NumberOfFFTFrequencies;
  size_t NumberOfSplinedFrequencies;

  int WaveformType;
  std::vector<double>  WaveformTime;
  std::vector<double>  WaveformCurrent;
  std::vector<double>  WaveformReceived;

  std::vector<double>  T_Waveform;
  std::vector<cdouble> F_Waveform;
  std::vector<cdouble> Transfer;
  std::vector<cdouble> FFTWork;
  std::vector<double>  fft_frequency;    
  fftw_plan       fftwplan_backward;  
	  
  size_t FrequenciesPerDecade;
  size_t NumberOfDiscreteFrequencies;  
  std::vector<double> DiscreteFrequencies;
  std::vector<double> DiscreteFrequenciesLog10;  
  double DiscreteFrequencyLow;
  double DiscreteFrequencyHigh;
  double FrequencyLog10Spacing;  
  std::vector<double> LowPassFilterCutoffFrequency;
  std::vector<double> LowPassFilterOrder;
  std::vector<double> SplinedFrequencieslog10;;

  double XOutputScaling;
  double YOutputScaling;
  double ZOutputScaling;

  double TX_LoopArea;
  double TX_NumberOfTurns;
  double TX_PeakCurrent;

  double TX_height;
  double TX_pitch;
  double TX_roll;
  cVec TX_orientation;

  double RX_height;
  double RX_pitch;
  double RX_roll;  

  cVec TX_RX_separation;
      
  void readsystemdescriptorfile(std::string systemdescriptorfile);
  void systeminitialise();
  void forwardmodel(const std::vector<double>& conductivity, const std::vector<double>& thickness, const sTDEmGeometry& geometry);
  void setconductivitythickness(const std::vector<double>& conductivity, const std::vector<double>& thickness);

  void createwaveform();
  void setupdiscretefrequencies();
  double calculate_fft_frequency(size_t index);
  void setgeometry(const sTDEmGeometry& G);
  void setgeometry(const double& tx_height, const double& tx_pitch, const double& tx_roll, const double& rx_pitch, const double& rx_roll, const double& txrx_dx, const double& txrx_dy, const double& txrx_dz);  
  void setupcomputations();					
  void setprimaryfields();
  void setsecondaryfields();
  void inversefft(){fftw_execute(fftwplan_backward);}
      
  void initialise_windows();
  void initialise_windows_area();
  void initialise_windows_boxcar();
  void initialise_windows_lineartaper();

  void computewindow(double* timeseries, std::vector<double>& W);
  std::vector<WindowSpecification> WinSpec;

  void printwindows();
  void write_windows(const std::string& path);
  void write_timedomainwaveform(const std::string& path);
  void write_frequencydomainwaveform(const std::string& path);
  void write_discretefrequencies(const std::string& path);
  void write_splinedfrequencies(const std::string& path);
  void write_frequencyseries(const std::string& path);
  void write_timesseries(const std::string& path);
      
  std::vector<double> X; //Secondary X field
  std::vector<double> Y; //Secondary Y field 
  std::vector<double> Z; //Secondary Z field 
  double PrimaryX;  //Primary X field
  double PrimaryY;  //Primary Y field
  double PrimaryZ;  //Primary Z field

  double PrimaryXforward;  //Primary X field (forwrad model only) for ppm normalization
  double PrimaryYforward;  //Primary Y field (forwrad model only)
  double PrimaryZforward;  //Primary Z field (forwrad model only)

  double primary(const size_t component){
	  if (component == 0) return PrimaryX;
	  else if (component == 1) return PrimaryY;
	  else if (component == 2) return PrimaryZ;
	  else return 0;
  }

  double secondary(const size_t component, const size_t window){
	  if (component == 0) return X[window];
	  else if (component == 1) return Y[window];
	  else if (component == 2) return Z[window];
	  else return 0;
  }
      
  std::vector<double> HxR;
  std::vector<double> HxI;
  std::vector<double> HyR;
  std::vector<double> HyI;
  std::vector<double> HzR;
  std::vector<double> HzI;
  std::vector<double> HxR_spline;
  std::vector<double> HxI_spline;
  std::vector<double> HyR_spline;
  std::vector<double> HyI_spline;
  std::vector<double> HzR_spline;
  std::vector<double> HzI_spline;

  std::vector<double> h2_spline;
  std::vector<double> a_spline;
  std::vector<double> b_spline;
  std::vector<size_t> klo_spline;
  std::vector<size_t> khi_spline;
  std::vector<double> a3ma_spline;
  std::vector<double> b3mb_spline;

  std::vector<cdouble> X_splined;
  std::vector<cdouble> Y_splined;
  std::vector<cdouble> Z_splined;
 
  void setup_splines();
  void spline(const std::vector<double>& x, const std::vector<double>& y, double yp1, double ypn, std::vector<double>& y2);
  void setup_splineinterp(const std::vector<double>& xn, const std::vector<double>& xi);  
  void splineinterp();
  
  
  cVec rotatetoreceiverorientation(cVec v);
  
  dmatrix readwaveformfile(const std::string& filename);
  void digitisewaveform(const dmatrix& w, std::vector<double>& t, std::vector<double>& v);
	  
  void drx_pitch(const double xb, const double zb, const double p, double& dxbdp, double& dzbdp);
  void drx_roll(const double yb, const double zb, const double r, double& dybdr, double& dzbdr);

  void drx_pitch(const std::vector<double>& xb, const std::vector<double>& zb, const double p, std::vector<double>& dxbdp, std::vector<double>& dzbdp);
  void drx_roll( const std::vector<double>& yb, const std::vector<double>& zb, const double r, std::vector<double>& dybdr, std::vector<double>& dzbdr);


};
////////////////////////////////////////////////////////////////////////////

#endif
