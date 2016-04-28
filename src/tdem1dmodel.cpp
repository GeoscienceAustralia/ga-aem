/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include "tdem1dmodel.h"
#include "tdemsystem.h"

void* tdem_create(const char* systemdescriptorfile)
{		
	cTDEmSystem* T;			
	T = new cTDEmSystem;	
	T->readsystemdescriptorfile(systemdescriptorfile); 
	return (void*) T;		
}

void tdem_delete(void* S)
{
	cTDEmSystem* T = (cTDEmSystem*)S;
	delete T;
}

int tdem_nwindows(void* S)
{    
	cTDEmSystem* T = (cTDEmSystem*)S;
	return (int)T->NumberOfWindows;
}

void tdem_1Dmodel(void* S, sTDEmGeometry& G, cEarth1D& E, sTDEmResponse& R)
{		
	cTDEmSystem* T = (cTDEmSystem*)S;				
	
	T->setconductivitythickness(E.conductivity,E.thickness);
	T->setgeometry(G.tx_height,G.tx_pitch,G.tx_roll,G.rx_pitch,G.rx_roll,G.txrx_dx,G.txrx_dy,G.txrx_dz);		
	T->setupcomputations();	
	T->setprimaryfields();	
	T->setsecondaryfields();

	R.PX=T->PrimaryX;
	R.PY=T->PrimaryY;
	R.PZ=T->PrimaryZ;

	R.SX.resize(T->NumberOfWindows);
	R.SY.resize(T->NumberOfWindows);
	R.SZ.resize(T->NumberOfWindows);
	for(size_t i=0; i<T->NumberOfWindows; i++){
		R.SX[i] = T->X[i];
		R.SY[i] = T->Y[i];
		R.SZ[i] = T->Z[i];						
	}	
}

void tdem_1Dderivative(void* S, eCalculationType ctype, size_t derivate_layerindex, sTDEmResponse& R)
{		
	cTDEmSystem* T = (cTDEmSystem*)S;	
	
	T->Earth.calculation_type = ctype;
	T->Earth.derivative_layer = derivate_layerindex;
	T->setprimaryfields();
	T->setsecondaryfields();

	R.PX=T->PrimaryX;
	R.PY=T->PrimaryY;
	R.PZ=T->PrimaryZ;

	R.SX.resize(T->NumberOfWindows);
	R.SY.resize(T->NumberOfWindows);
	R.SZ.resize(T->NumberOfWindows);
	for(size_t i=0; i<T->NumberOfWindows; i++){
		R.SX[i] = T->X[i];
		R.SY[i] = T->Y[i];
		R.SZ[i] = T->Z[i];				
	}
}
