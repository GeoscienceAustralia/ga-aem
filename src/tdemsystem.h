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

enum eWaveFormType { WT_TX, WT_RX };
enum eOutputType { OT_BFIELD, OT_DBDT};
enum eNormalizationType { NT_NONE, NT_PPM, NT_PPM_PEAKTOPEAK};

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

class cTDEmResponse{

public:
	double PX;
	double PY;
	double PZ;
	std::vector<double> SX;
	std::vector<double> SY;
	std::vector<double> SZ;
};

class cTDEmGeometry{

public:
	double tx_height = 0.0;
	double tx_roll   = 0.0;
	double tx_pitch  = 0.0;
	double tx_yaw    = 0.0;
	double txrx_dx   = 0.0;
	double txrx_dy   = 0.0;
	double txrx_dz   = 0.0;
	double rx_roll   = 0.0;
	double rx_pitch  = 0.0;
	double rx_yaw    = 0.0;

	cTDEmGeometry(){};

	cTDEmGeometry(const double& _tx_height, const double& _tx_roll, const double& _tx_pitch, const double& _tx_yaw, const double& _txrx_dx, const double& _txrx_dy, const double& _txrx_dz, const double& _rx_roll, const double& _rx_pitch, const double& _rx_yaw){
		tx_height = _tx_height;
		tx_roll   = _tx_roll; tx_pitch = _tx_pitch; tx_yaw = _tx_yaw;
		txrx_dx   = _txrx_dx; txrx_dy  = _txrx_dy; txrx_dz = _txrx_dz;
		rx_roll   = _rx_roll; rx_pitch = _rx_pitch; rx_yaw = _rx_yaw;
	}

	cTDEmGeometry(const cBlock& b){
		set_zero();
		b.getvalue("tx_height", tx_height);
		b.getvalue("tx_roll", tx_roll);
		b.getvalue("tx_pitch", tx_pitch);
		b.getvalue("tx_yaw", tx_yaw);
		b.getvalue("txrx_dx", txrx_dx);
		b.getvalue("txrx_dy", txrx_dy);
		b.getvalue("txrx_dz", txrx_dz);
		b.getvalue("rx_roll", rx_roll);
		b.getvalue("rx_pitch", rx_pitch);
		b.getvalue("rx_yaw", rx_yaw);
	}

	void set_zero(){
		tx_height = 0.0;
		tx_roll = 0.0; tx_pitch = 0.0; tx_yaw = 0.0;
		txrx_dx = 0.0; txrx_dy = 0.0;  txrx_dz = 0.0;
		rx_roll = 0.0; rx_pitch = 0.0; rx_yaw = 0.0;
	}

	void fillundefined(const cTDEmGeometry& g)
	{
		if (g.tx_height != cBlock::ud_double())tx_height = g.tx_height;
		if (g.tx_roll != cBlock::ud_double())tx_roll = g.tx_roll;
		if (g.tx_pitch != cBlock::ud_double())tx_pitch = g.tx_pitch;
		if (g.tx_yaw != cBlock::ud_double())tx_yaw = g.tx_yaw;
		if (g.txrx_dx != cBlock::ud_double())txrx_dx = g.txrx_dx;
		if (g.txrx_dy != cBlock::ud_double())txrx_dy = g.txrx_dy;
		if (g.txrx_dz != cBlock::ud_double())txrx_dz = g.txrx_dz;
		if (g.rx_roll != cBlock::ud_double())rx_roll = g.rx_roll;
		if (g.rx_pitch != cBlock::ud_double())rx_pitch = g.rx_pitch;
		if (g.rx_yaw != cBlock::ud_double())rx_yaw = g.rx_yaw;
	}
};

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

private:
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

  eOutputType OutputType;
  eNormalizationType Normalisation;
  
  size_t NumberOfWindows;
  std::string WindowWeightingScheme;

  double BaseFrequency;
  double BasePeriod;
  double SampleFrequency;  
  double SampleInterval;  
  size_t SamplesPerWaveform;
  size_t NumberOfFFTFrequencies;
  size_t NumberOfSplinedFrequencies;

  eWaveFormType WaveformType;
  std::vector<double>  WaveformTime;
  std::vector<double>  WaveformCurrent;
  std::vector<double>  WaveformReceived;

  std::vector<double>  T_Waveform;
  std::vector<cdouble> F_Waveform;
  std::vector<cdouble> Transfer;
  std::vector<cdouble> FFTWork;
  std::vector<double>  fft_frequency;    
  fftw_plan fftwplan_backward;  
	  
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

  double TX_LoopArea;
  double TX_NumberOfTurns;
  double TX_PeakCurrent;
  double TX_PeakdIdT;
  double TX_height;
  double TX_pitch;
  double TX_roll;
  cVec TX_orientation;
  double RX_height;
  double RX_pitch;
  double RX_roll;  
  cVec TX_RX_separation;

  cTDEmGeometry NormalizationGeometry;
  double XScale;
  double YScale;
  double ZScale;

  std::vector<double> X; //Secondary X field
  std::vector<double> Y; //Secondary Y field 
  std::vector<double> Z; //Secondary Z field 
  double PrimaryX;  //Primary X field
  double PrimaryY;  //Primary Y field
  double PrimaryZ;  //Primary Z field
      
  void readsystemdescriptorfile(std::string systemdescriptorfile);
  void systeminitialise();
  void forwardmodel(const std::vector<double>& conductivity, const std::vector<double>& thickness, const cTDEmGeometry& geometry);
  void setconductivitythickness(const std::vector<double>& conductivity, const std::vector<double>& thickness);

  double compute_peak_didt();
  void setup_scaling();
  void createwaveform();
  void setupdiscretefrequencies();
  double calculate_fft_frequency(size_t index);
  void setgeometry(const cTDEmGeometry& G);
  void setgeometry(const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw);
  void setupcomputations();					
  void setprimaryfields();
  void setsecondaryfields();
  void inversefft(){fftw_execute(fftwplan_backward);}
      
  void initialise_windows();
  void initialise_windows_area();
  void initialise_windows_boxcar();
  void initialise_windows_lineartaper();

  void computewindow(const double* timeseries, std::vector<double>& W);
  std::vector<WindowSpecification> WinSpec;

  void printwindows();
  void write_windows(const std::string& path);
  void write_timedomainwaveform(const std::string& path);
  void write_frequencydomainwaveform(const std::string& path);
  void write_discretefrequencies(const std::string& path);
  void write_splinedfrequencies(const std::string& path);
  void write_frequencyseries(const std::string& path);
  void write_timesseries(const std::string& path);
      
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
     
  void setup_splines();
  void spline(const std::vector<double>& x, const std::vector<double>& y, double yp1, double ypn, std::vector<double>& y2);
  void setup_splineinterp(const std::vector<double>& xn, const std::vector<double>& xi);  
  void spline_interp();  
  cVec rotatetoreceiverorientation(cVec v);  
  dmatrix readwaveformfile(const std::string& filename);
  void digitisewaveform(const dmatrix& w, std::vector<double>& t, std::vector<double>& v);	  
  void drx_pitch(const double xb, const double zb, const double p, double& dxbdp, double& dzbdp);
  void drx_roll(const double yb, const double zb, const double r, double& dybdr, double& dzbdr);
  void drx_pitch(const std::vector<double>& xb, const std::vector<double>& zb, const double p, std::vector<double>& dxbdp, std::vector<double>& dzbdp);
  void drx_roll( const std::vector<double>& yb, const std::vector<double>& zb, const double r, std::vector<double>& dybdr, std::vector<double>& dzbdr);
  void forwardmodel(const cTDEmGeometry& G, const cEarth1D& E, cTDEmResponse& R);
};
////////////////////////////////////////////////////////////////////////////

#endif
