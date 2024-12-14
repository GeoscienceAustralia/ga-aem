/*
This source code file is licensed under the GNU GPL Version 2.0 Licence by the following copyright holder:
Crown Copyright Commonwealth of Australia (Geoscience Australia) 2015.
The GNU GPL 2.0 licence is available at: http://www.gnu.org/licenses/gpl-2.0.html. If you require a paper copy of the GNU GPL 2.0 Licence, please write to Free Software Foundation, Inc. 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.

Author: Ross C. Brodie, Geoscience Australia.
*/

#include <stdio.h>
#include "tdemsystem.hpp"

#define  EXPORT_FCNS
#include "shrhelp.hpp"
#include "earth1d.hpp"
#include "gatdaem1d.h"


class cLogger glog; //The global instance of the log file manager
using namespace AEM;

void* createhandle(const char* systemfile)
{
	cTDEmSystem* T = new cTDEmSystem(std::string(systemfile));
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
	return (int)T.waveform().NumSamples;
}

void waveform(void* hS, double* time, double* currentwaveform, double* voltagewaveform)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	for (size_t i = 0; i < T.waveform().NumSamples; i++) {
		time[i] = T.waveform().Time[i];
		if (T.waveform().Type == Waveform::Type::TX) {
			currentwaveform[i] = T.waveform().TD_Waveform[i];
		}
		if (T.waveform().Type == Waveform::Type::RX) {
			voltagewaveform[i] = T.waveform().TD_Waveform[i];
		}
	}
}

int nwindows(void* hS)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	return (int)T.nwindows();
}

int nturns(void* hS)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	return (int)T.transmitter().NumberOfTurns;
}

double peakcurrent(void* hS)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	return (double)T.transmitter().PeakCurrent;
}

double looparea(void* hS)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	return (double)T.transmitter().LoopArea;
}

double basefrequency(void* hS)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	return (double)T.waveform().BaseFrequency;
}

int nlayers(void* hS)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	return (int)T.lem().nLayers();
}

void windowtimes(void* hS, double* low, double* high)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	for (size_t i = 0; i < T.nwindows(); i++) {
		low[i] = T.window(i).TimeLow;
		high[i] = T.window(i).TimeHigh;
	}
}

void setgeometry(void* hS, const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	const cTDEmGeometry G(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
	T.setgeometry(G);
}

void setearth(void* hS, int nlayers, double* conductivity, double* thickness)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	Earth1D E(nlayers, conductivity, thickness);
	T.lem().set_earth(E);
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
	cTDEmGeometry G(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
	T.setgeometry(G);
	Earth1D E(nlayers, conductivity, thickness);
	T.lem().set_earth(E);
	T.setup_computations();
	T.lem().set_calculationtype(CMode::FM);
	T.setprimaryfields();
	T.setsecondaryfields();

	size_t nw = T.nwindows();
	size_t sz = sizeof(double) * nw;

	*PX = T.PX();
	*PY = T.PY();
	*PZ = T.PZ();
	memcpy(SX, T.XS().data(), sz);
	memcpy(SY, T.YS().data(), sz);
	memcpy(SZ, T.ZS().data(), sz);
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
	cTDEmGeometry G(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
	T.setgeometry(G);
	Earth1D E(nlayers, conductivity, thickness, chargeability, timeconstant, frequencydependence);
	T.lem().set_iptype((AEM::IPType)iptype);
	T.lem().set_earth(E);
	T.setup_computations();
	T.lem().set_calculationtype(CMode::FM);
	T.setprimaryfields();
	T.setsecondaryfields();

	size_t nw = T.nwindows();
	size_t sz = sizeof(double) * nw;

	*PX = T.PX();
	*PY = T.PY();
	*PZ = T.PZ();
	memcpy(SX, T.XS().data(), sz);
	memcpy(SY, T.YS().data(), sz);
	memcpy(SZ, T.ZS().data(), sz);
}

void derivative(void* hS, int dtype, int dlayer, double* PX, double* PY, double* PZ, double* SX, double* SY, double* SZ) {
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	CMode mode = CalculationType::lookup_mode((size_t)dtype);
	CalculationType calc(mode, dlayer - 1);
	T.lem().set_calculationtype(calc);//subtract one from the layer number for zero based indexing
	T.setprimaryfields();
	T.setsecondaryfields();

	size_t nw = T.nwindows();
	size_t sz = sizeof(double) * nw;

	*PX = T.PX();
	*PY = T.PY();
	*PZ = T.PZ();
	memcpy(SX, T.XS().data(), sz);
	memcpy(SY, T.YS().data(), sz);
	memcpy(SZ, T.ZS().data(), sz);
}

void fm_dlogc(void* hS,
	const double tx_height, 
	const double tx_roll, const double tx_pitch, const double tx_yaw,
	const double txrx_dx, const double txrx_dy, const double txrx_dz, 
	const double rx_roll, const double rx_pitch, const double rx_yaw,
	const int nlayers, 
	const double* conductivity, 
	const double* thickness,
	double* R)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	cTDEmGeometry G(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
	T.setgeometry(G);
	Earth1D E(nlayers, conductivity, thickness);
	T.lem().set_earth(E);
	T.setup_computations();
	T.lem().set_calculationtype(CMode::FM);
	T.setprimaryfields();
	T.setsecondaryfields();

	size_t nw = T.nwindows();
	size_t sz = sizeof(double) * nw;
	double* p = R;
	*p = T.PX(); p++;
	memcpy(p, T.XS().data(), sz); p += nw;
	*p = T.PY(); p++;
	memcpy(p, T.YS().data(), sz); p += nw;
	*p = T.PZ(); p++;
	memcpy(p, T.ZS().data(), sz); p += nw;

	for (size_t k = 0; k < (size_t)nlayers; k++) {
		T.lem().set_calculationtype(CalculationType(CMode::DC, k));
		T.setprimaryfields();
		T.setsecondaryfields();

		//const double& c = T.lem().layers()[k].Conductivity;
		const double& c = conductivity[k];
		*p = T.PX() * c; p++;
		for (size_t w = 0; w < nw; w++) {
			*p = T.XS()[w] * c;
			p++;
		}

		*p = T.PY() * c; p++;
		for (size_t w = 0; w < nw; w++) {
			*p = T.YS()[w] * c;
			p++;
		}

		*p = T.PZ() * c; p++;
		for (size_t w = 0; w < nw; w++) {
			*p = T.ZS()[w] * c;
			p++;
		}

	}
}

void derivative_rx_pitch(void* hS, int n, double rx_pitch, double* xb, double* zb, double* dxbdp, double* dzbdp)
{
	cTDEmSystem& T = *(cTDEmSystem*)hS;
	//T.lem().calculation_type = cLEM::CalculationType::FM;
	//T.lem().derivative_layer = -1;
	//T.set_primaryfields();
	//T.set_secondaryfields();

	size_t sz = sizeof(double) * n;

	std::vector<double> dxbdpvec(n);
	std::vector<double> dzbdpvec(n);
	std::vector<double> xbvec(xb, xb + n);
	std::vector<double> zbvec(zb, zb + n);

	T.drx_pitch(xbvec, zbvec, rx_pitch, dxbdpvec, dzbdpvec);

	memcpy(dxbdp, dxbdpvec.data(), sz);
	memcpy(dzbdp, dzbdpvec.data(), sz);
}


