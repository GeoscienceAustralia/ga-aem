/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <stdio.h>
#include "tdemsystem.h"

#ifndef TESTOUTSIDEMATLAB
#define  EXPORT_FCNS
#include "shrhelp.h"
#endif

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
	for(int i=0; i<T.SamplesPerWaveform; i++){
		time[i]             = T.WaveformTime[i];
		if(T.WaveformType == WAVEFORMTYPE_CURRENT){
			currentwaveform[i]  = T.WaveformCurrent[i];
		}

		if(T.WaveformType == WAVEFORMTYPE_RECEIVER){
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
	return (int)T.Earth.NumLayers;	
}

void windowtimes(void* hS, double* low, double* high)
{		
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	for(int i=0; i<T.NumberOfWindows; i++){
		low[i]   = T.WinSpec[i].TimeLow;
		high[i]  = T.WinSpec[i].TimeHigh;
	}
}

void setgeometry(void* hS, struct sTDEmGeometry G)
{		
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	T.setgeometry(G);	
}

void setearth(void* hS, int nlayers, double* conductivity, double* thickness)
{		
	cTDEmSystem& T = *(cTDEmSystem*)hS;	
	T.Earth.setconductivitythickness(nlayers,conductivity,thickness);
}

void forwardmodel(void* hS, struct sTDEmGeometry G, int nlayers, double* conductivity, double* thickness, struct sTDEmResponseML* pR)
{		
	cTDEmSystem& T = *(cTDEmSystem*)hS;	
	T.setgeometry(G);
	T.Earth.setconductivitythickness(nlayers,conductivity,thickness);
	T.setupcomputations();
	T.Earth.calculation_type = CT_FORWARDMODEL;
	T.Earth.derivative_layer = -1;			
	T.setprimaryfields();	
	T.setsecondaryfields();

	size_t nw=T.NumberOfWindows;
	size_t sz=sizeof(double)*nw;

	sTDEmResponseML& R=*pR;
	R.PX=T.PrimaryX;
	R.PY=T.PrimaryY;
	R.PZ=T.PrimaryZ;		
	memcpy(R.SX,&(T.X[0]),sz);
	memcpy(R.SY,&(T.Y[0]),sz);
	memcpy(R.SZ,&(T.Z[0]),sz);
}

void derivative(void* hS, int dtype, int dlayer, struct sTDEmResponseML* pR)
{		
	cTDEmSystem& T = *(cTDEmSystem*)hS;		
	T.Earth.calculation_type = (eCalculationType)dtype;
	T.Earth.derivative_layer = dlayer-1;	//subtract one from the layer number for zero based indexing
	T.setprimaryfields();	
	T.setsecondaryfields();

	size_t nw=T.NumberOfWindows;
	size_t sz=sizeof(double)*nw;

	sTDEmResponseML& R=*pR;
	R.PX=T.PrimaryX;
	R.PY=T.PrimaryY;
	R.PZ=T.PrimaryZ;
	memcpy(R.SX,&(T.X[0]),sz);
	memcpy(R.SY,&(T.Y[0]),sz);
	memcpy(R.SZ,&(T.Z[0]),sz);
}

void fm_dlogc(void* hS, struct sTDEmGeometry G, int nlayers, double* conductivity, double* thickness, double* R)
{				
	cTDEmSystem& T = *(cTDEmSystem*)hS;				
	T.setgeometry(G);
	T.Earth.setconductivitythickness(nlayers,conductivity,thickness);		
	T.setupcomputations();
	T.Earth.calculation_type = CT_FORWARDMODEL;		
	T.Earth.derivative_layer = -1;	
	T.setprimaryfields();	
	T.setsecondaryfields();

	size_t nw=T.NumberOfWindows;
	size_t sz=sizeof(double)*nw;
	double* p = R;
	*p=T.PrimaryX; p++;
	memcpy(p,&(T.X[0]),sz); p+=nw;
	*p=T.PrimaryY; p++;
	memcpy(p,&(T.Y[0]),sz); p+=nw;
	*p=T.PrimaryZ; p++;			
	memcpy(p,&(T.Z[0]),sz); p+=nw;	

	for(size_t k=0;k<nlayers;k++){
		T.Earth.calculation_type = CT_CONDUCTIVITYDERIVATIVE;		
		T.Earth.derivative_layer = (int)k;
		T.setprimaryfields();	
		T.setsecondaryfields();
		
		double c=T.Earth.Layer[k].Conductivity;
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











