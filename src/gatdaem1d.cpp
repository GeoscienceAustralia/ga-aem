/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <stdio.h>
#include "tdemsystem.h"

#define  EXPORT_FCNS
#include "shrhelp.h"
#include "earth1d.h"
#include "gatdaem1d.h"

void* createhandle(const char* systemfile)
{			
	cTDEmSystem* T = new cTDEmSystem;
	T->readsystemdescriptorfile(std::string(systemfile));
	return (void*)T;
}

void deletehandle(void* hS)
{		
	cTDEmSystem* T = (cTDEmSystem*)hS;
	delete T;		
}

int nsamplesperwaveform(void* hS)
{		
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	return (int)T.SamplesPerWaveform;	
}

void waveform(void* hS, double* time, double* currentwaveform, double* voltagewaveform)
{			
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	for(size_t i=0; i<T.SamplesPerWaveform; i++){
		time[i] = T.WaveformTime[i];
		if(T.WaveformType == WT_TX){
			currentwaveform[i]  = T.WaveformCurrent[i];
		}

		if(T.WaveformType == WT_RX){
			voltagewaveform[i]  = T.WaveformReceived[i];
		}
	}		
}

int nwindows(void* hS)
{		
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	return (int)T.NumberOfWindows;	
}

int nlayers(void* hS)
{		
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	return (int)T.LEM.NumLayers;	
}

void windowtimes(void* hS, double* low, double* high)
{		
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	for(size_t i=0; i<T.NumberOfWindows; i++){
		low[i]   = T.WinSpec[i].TimeLow;
		high[i]  = T.WinSpec[i].TimeHigh;
	}
}

void setgeometry(void* hS, const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	T.setgeometry(tx_height,tx_roll,tx_pitch,tx_yaw,txrx_dx,txrx_dy,txrx_dz,rx_roll,rx_pitch,rx_yaw);
}

void setearth(void* hS, int nlayers, double* conductivity, double* thickness)
{		
	cTDEmSystem& T = *(cTDEmSystem*)hS;			
	T.LEM.setconductivitythickness(nlayers,conductivity,thickness);
}

void forwardmodel(void* hS,
	const double tx_height,
	const double tx_roll,
	const double tx_pitch,
	const double tx_yaw,
	const double txrx_dx,
	const double txrx_dy,
	const double txrx_dz,
	const double rx_roll,
	const double rx_pitch,
	const double rx_yaw,
	const int nlayers,
	const double* conductivity,
	const double* thickness,	
	double* PX,
	double* PY,
	double* PZ,
	double* SX,
	double* SY,
	double* SZ)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	T.setgeometry(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
	cEarth1D E(nlayers, conductivity, thickness);
	T.LEM.setproperties(E);
	T.setupcomputations();
	T.LEM.calculation_type = CT_FORWARDMODEL;
	T.LEM.derivative_layer = -1;
	T.setprimaryfields();
	T.setsecondaryfields();

	size_t nw = T.NumberOfWindows;
	size_t sz = sizeof(double)*nw;

	*PX = T.PrimaryX;
	*PY = T.PrimaryY;
	*PZ = T.PrimaryZ;
	memcpy(SX, T.X.data(), sz);
	memcpy(SY, T.Y.data(), sz);
	memcpy(SZ, T.Z.data(), sz);
}

void forwardmodel_ip(void* hS,
	const double tx_height,
	const double tx_roll,
	const double tx_pitch,
	const double tx_yaw,
	const double txrx_dx,
	const double txrx_dy,
	const double txrx_dz,
	const double rx_roll,
	const double rx_pitch,
	const double rx_yaw,
	const int nlayers,
	const double* conductivity,
	const double* thickness,
	const int iptype,
	const double* chargeability,
	const double* timeconstant,
	const double* frequencydependence,
	double* PX,
	double* PY,
	double* PZ,
	double* SX,
	double* SY,
	double* SZ)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	T.setgeometry(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
	cEarth1D E(nlayers, conductivity, thickness, chargeability, timeconstant, frequencydependence);
	T.LEM.iptype = (eIPType)iptype;
	T.LEM.setproperties(E);
	T.setupcomputations();
	T.LEM.calculation_type = CT_FORWARDMODEL;
	T.LEM.derivative_layer = -1;
	T.setprimaryfields();
	T.setsecondaryfields();

	size_t nw = T.NumberOfWindows;
	size_t sz = sizeof(double)*nw;

	*PX = T.PrimaryX;
	*PY = T.PrimaryY;
	*PZ = T.PrimaryZ;
	memcpy(SX, T.X.data(), sz);
	memcpy(SY, T.Y.data(), sz);
	memcpy(SZ, T.Z.data(), sz);
}

void derivative(void* hS, int dtype, int dlayer,
	double* PX, double* PY, double* PZ, double* SX, double* SY, double* SZ)
{		
	cTDEmSystem& T = *(cTDEmSystem*)hS;		
	T.LEM.calculation_type = (eCalculationType)dtype;
	T.LEM.derivative_layer = dlayer-1;	//subtract one from the layer number for zero based indexing
	T.setprimaryfields();	
	T.setsecondaryfields();

	size_t nw=T.NumberOfWindows;
	size_t sz=sizeof(double)*nw;
	
	*PX = T.PrimaryX;
	*PY = T.PrimaryY;
	*PZ = T.PrimaryZ;
	memcpy(SX, T.X.data(), sz);
	memcpy(SY, T.Y.data(), sz);
	memcpy(SZ, T.Z.data(), sz);
}

void fm_dlogc(void* hS,
	const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw,
	const int nlayers, const double* conductivity, const double* thickness,
	double* R)
{				
	cTDEmSystem& T = *(cTDEmSystem*)hS;				
	T.setgeometry(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
	T.LEM.setconductivitythickness(nlayers,conductivity,thickness);		
	T.setupcomputations();
	T.LEM.calculation_type = CT_FORWARDMODEL;		
	T.LEM.derivative_layer = -1;	
	T.setprimaryfields();	
	T.setsecondaryfields();

	size_t nw=T.NumberOfWindows;
	size_t sz=sizeof(double)*nw;
	double* p = R;
	*p=T.PrimaryX; p++;
	memcpy(p,T.X.data(),sz); p+=nw;
	*p=T.PrimaryY; p++;
	memcpy(p,T.Y.data(),sz); p+=nw;
	*p=T.PrimaryZ; p++;			
	memcpy(p,T.Z.data(),sz); p+=nw;	

	for(size_t k=0;k<(size_t)nlayers;k++){
		T.LEM.calculation_type = CT_CONDUCTIVITYDERIVATIVE;		
		T.LEM.derivative_layer = k;
		T.setprimaryfields();	
		T.setsecondaryfields();
		
		double c=T.LEM.Layer[k].Conductivity;
		*p=T.PrimaryX*c; p++;		
		for(size_t w=0;w<nw;w++){
			*p = T.X[w] * c;
			p++;			
		}

		*p=T.PrimaryY*c; p++;		
		for(size_t w=0;w<nw;w++){
			*p = T.Y[w] * c;
			p++;			
		}

		*p=T.PrimaryZ*c; p++;		
		for(size_t w=0;w<nw;w++){
			*p = T.Z[w] * c;
			p++;			
		}	

	}
}

