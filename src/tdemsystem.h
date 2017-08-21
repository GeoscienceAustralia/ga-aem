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
#include "earth1d.h"
#include "lem.h"

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

enum eGeometryElementType { 
	GE_TX_HEIGHT,
	GE_TX_ROLL,
	GE_TX_PITCH,
	GE_TX_YAW,
	GE_TXRX_DX,
	GE_TXRX_DY,
	GE_TXRX_DZ,
	GE_RX_ROLL,
	GE_RX_PITCH,
	GE_RX_YAW
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

	cTDEmGeometry(const std::vector<double> gvector){
		for (size_t i = 0; i < size(); i++){
			(*this)[i] = gvector[i];
		}		
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

	inline static size_t size(){
		return 10;
	}

	double& operator[](const size_t& index){

		switch (index) {
		case 0: return tx_height; break;
		case 1: return tx_roll; break;
		case 2: return tx_pitch; break;
		case 3: return tx_yaw; break;
		case 4: return txrx_dx; break;
		case 5: return txrx_dy; break;
		case 6: return txrx_dz; break;
		case 7: return rx_roll; break;
		case 8: return rx_pitch; break;
		case 9: return rx_yaw; break;
		default:
			rootmessage("Geometry index %llu out of range\n", index);
			std::string e = strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(e));
			break;
		}
	}

	double operator[](const size_t& index) const {
		return (*this)[index];
	};

	void set_zero(){
		for (size_t i = 0; i < size(); i++){
			(*this)[i] = 0.0;
		}		
	}

	void fillundefined(const cTDEmGeometry& g)
	{				
		for (size_t i = 0; i < size(); i++){			
			if ((*this)[i] == cBlock::ud_double()){
				(*this)[i] = g[i];				 
			}
		}		
	}	
	
	static std::string fname(const size_t& index){

		switch (index) {
		case 0: return "tx_height"; break;
		case 1: return "tx_roll"; break;
		case 2: return "tx_pitch"; break;
		case 3: return "tx_yaw"; break;
		case 4: return "txrx_dx"; break;
		case 5: return "txrx_dy"; break;
		case 6: return "txrx_dz"; break;
		case 7: return "rx_roll"; break;
		case 8: return "rx_pitch"; break;
		case 9: return "rx_yaw"; break;
		default:
			rootmessage("Geometry index %llu out of range\n", index);
			std::string e = strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(e));
			break;
		}
	};

	static size_t findex(const std::string& name){

		for (size_t i = 0; i < size(); i++){
			if (strcasecmp(name, fname(i)) == 0) return i;
		}
				
		rootmessage("Geometry field name %s is bad\n", name);
		std::string e = strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
		throw(std::runtime_error(e));		

		return 0;
	};

	static std::string units(const size_t& index){

		switch (index) {
		case 0: return "m"; break;
		case 1: return "degrees"; break;
		case 2: return "degrees"; break;
		case 3: return "degrees"; break;
		case 4: return "m"; break;
		case 5: return "m"; break;
		case 6: return "m"; break;
		case 7: return "degrees"; break;
		case 8: return "degrees"; break;
		case 9: return "degrees"; break;
		default:
			rootmessage("Geometry index %llu out of range\n", index);
			std::string e = strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(e));
			break;
		}
	};

	static std::string description(const size_t& index){

		switch (index) {
		case 0: return "Tx height above ground level"; break;
		case 1: return "Tx roll - left side up + ve";   break;
		case 2: return "Tx pitch - nose down + ve";  break;
		case 3: return "Tx yaw - turn left + ve";    break;
		case 4: return "Tx - Rx horizonatl inline separation";   break;
		case 5: return "Tx - Rx horizonatl transverse separation";   break;
		case 6: return "Tx - Rx vertical separation";   break;
		case 7: return "Rx roll - left side up + ve";   break;
		case 8: return "Rx pitch - nose down + ve";  break;
		case 9: return "Rx yaw - turn left + ve";    break;
		default:
			rootmessage("Geometry index %llu out of range\n", index);
			std::string e = strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(e));
			break;
		}
	};

	static eGeometryElementType elementtype(const size_t& index){		
		switch (index) {
		case 0: return GE_TX_HEIGHT; break;
		case 1: return GE_TX_ROLL;   break;
		case 2: return GE_TX_PITCH;  break;
		case 3: return GE_TX_YAW;    break;		
		case 4: return GE_TXRX_DX;   break;
		case 5: return GE_TXRX_DY;   break;
		case 6: return GE_TXRX_DZ;   break;
		case 7: return GE_RX_ROLL;   break;
		case 8: return GE_RX_PITCH;  break;
		case 9: return GE_RX_YAW;    break;
		default:
			rootmessage("Geometry index %llu out of range\n", index);
			std::string e = strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(e));
			break;
		}
	}

	static eCalculationType derivativetype(const size_t& index){		
		switch (index) {
		case 0: return CT_HDERIVATIVE; break;
		case 1: return CT_NONE; break;
		case 2: return CT_NONE; break;
		case 3: return CT_NONE; break;
		case 4: return CT_XDERIVATIVE; break;
		case 5: return CT_YDERIVATIVE; break;
		case 6: return CT_ZDERIVATIVE; break;
		case 7: return CT_NONE; break;		
		case 8: return CT_NONE; break;
		case 9: return CT_NONE; break;
		default:
			rootmessage("Geometry index %llu out of range\n", index);
			std::string e = strprint("Error: exception throw from %s (%d) %s\n", __FILE__, __LINE__, __FUNCTION__);
			throw(std::runtime_error(e));
			break;
		}
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
  cLEM LEM;
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
  
  void setearthproperties(const cEarth1D& E);
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

  std::vector<double> secondary(const size_t component){
	  if (component == 0) return X;
	  else if (component == 1) return Y;
	  else if (component == 2) return Z;
	  else return std::vector<double>(0);
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
