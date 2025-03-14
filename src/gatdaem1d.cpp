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
	TDEmSystem* T = new TDEmSystem(std::string(systemfile));
	return (void*)T;
}

void deletehandle(void* hS)
{
	TDEmSystem* T = (TDEmSystem*)hS;
	delete T;
}

int nsamplesperwaveform(void* hS)
{
	TDEmSystem& T = *(TDEmSystem*)hS;
	return (int)T.waveform().NumSamples;
}

void waveform(void* hS, double* time, double* currentwaveform, double* voltagewaveform)
{
	TDEmSystem& T = *(TDEmSystem*)hS;
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
	TDEmSystem& T = *(TDEmSystem*)hS;
	return (int)T.nWindows();
}

int nturns(void* hS)
{
	TDEmSystem& T = *(TDEmSystem*)hS;
	return (int)T.transmitter().nTurns;
}

double peakcurrent(void* hS)
{
	TDEmSystem& T = *(TDEmSystem*)hS;
	return (double)T.transmitter().PeakCurrent;
}

double looparea(void* hS)
{
	TDEmSystem& T = *(TDEmSystem*)hS;
	return (double)T.transmitter().LoopArea;
}

double basefrequency(void* hS)
{
	TDEmSystem& T = *(TDEmSystem*)hS;
	return (double)T.waveform().BaseFrequency;
}

int nlayers(void* hS)
{
	TDEmSystem& T = *(TDEmSystem*)hS;
	return (int)T.lem().nLayers();
}

void windowtimes(void* hS, double* low, double* high)
{
	TDEmSystem& T = *(TDEmSystem*)hS;
	for (size_t i = 0; i < T.nWindows(); i++) {
		low[i] = T.window(i).Low;
		high[i] = T.window(i).High;
	}
}

void setgeometry(void* hS, const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw)
{
	TDEmSystem& T = *(TDEmSystem*)hS;
	const TDEmGeometry G(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
	T.set_geometry(G);
}

void setearth(void* hS, int nlayers, double* conductivity, double* thickness)
{
	TDEmSystem& T = *(TDEmSystem*)hS;
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
	TDEmSystem& T = *(TDEmSystem*)hS;
	TDEmGeometry G(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
	Earth1D E(nlayers, conductivity, thickness);

	TDEmResponse R = T.forward_model(E, G);	
	size_t nw = T.nWindows();
	size_t sz = sizeof(double) * nw;

	memcpy(PX, R.primary(XCOMP).data(), sz);
	memcpy(PY, R.primary(YCOMP).data(), sz);
	memcpy(PZ, R.primary(ZCOMP).data(), sz);
	memcpy(SX, R.secondary(XCOMP).data(), sz);
	memcpy(SY, R.secondary(YCOMP).data(), sz);
	memcpy(SZ, R.secondary(ZCOMP).data(), sz);
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
	TDEmSystem& T = *(TDEmSystem*)hS;
	TDEmGeometry G(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
	Earth1D E(nlayers, conductivity, thickness, chargeability, timeconstant, frequencydependence);
	E.set_iptype((AEM::IPType)iptype);
	
	TDEmResponse R = T.forward_model(E, G);

	size_t nw = T.nWindows();
	size_t sz = sizeof(double) * nw;

	memcpy(PX, R.primary(XCOMP).data(), sz);
	memcpy(PY, R.primary(YCOMP).data(), sz);
	memcpy(PZ, R.primary(ZCOMP).data(), sz);

	memcpy(SX, R.secondary(XCOMP).data(), sz);
	memcpy(SY, R.secondary(YCOMP).data(), sz);
	memcpy(SZ, R.secondary(ZCOMP).data(), sz);	
}

void derivative(void* hS, int dtype, int dlayer, double* PX, double* PY, double* PZ, double* SX, double* SY, double* SZ) {
	TDEmSystem& T = *(TDEmSystem*)hS;
	CMode mode = CalculationType::lookup_mode((size_t)dtype);
	CalculationType calc(mode, dlayer - 1);//subtract one from the layer number for zero based indexing
	TDEmResponse R = T.derivative(calc);

	size_t nw = T.nWindows();
	size_t sz = sizeof(double) * nw;

	memcpy(PX, R.primary(XCOMP).data(), sz);
	memcpy(PY, R.primary(YCOMP).data(), sz);
	memcpy(PZ, R.primary(ZCOMP).data(), sz);
	memcpy(SX, R.secondary(XCOMP).data(), sz);
	memcpy(SY, R.secondary(YCOMP).data(), sz);
	memcpy(SZ, R.secondary(ZCOMP).data(), sz);
}

void fm_dlogc(void* hS,
	const double tx_height, 
	const double tx_roll, const double tx_pitch, const double tx_yaw,
	const double txrx_dx, const double txrx_dy, const double txrx_dz, 
	const double rx_roll, const double rx_pitch, const double rx_yaw,
	const int nlayers, 
	const double* conductivity, 
	const double* thickness,
	double* Response)
{
	TDEmSystem& T = *(TDEmSystem*)hS;
	TDEmGeometry G(tx_height, tx_roll, tx_pitch, tx_yaw, txrx_dx, txrx_dy, txrx_dz, rx_roll, rx_pitch, rx_yaw);
	T.set_geometry(G);
	Earth1D E(nlayers, conductivity, thickness);

	TDEmResponse R = T.forward_model(E, G);

	size_t nw = T.nWindows();
	size_t sz = sizeof(double) * nw;
	double* p = Response;
	memcpy(p, R.primary(XCOMP).data(), sz); p += nw;
	memcpy(p, R.secondary(XCOMP).data(), sz); p += nw;
	memcpy(p, R.primary(YCOMP).data(), sz); p += nw;
	memcpy(p, R.secondary(YCOMP).data(), sz); p += nw;
	memcpy(p, R.primary(ZCOMP).data(), sz); p += nw;
	memcpy(p, R.secondary(ZCOMP).data(), sz); p += nw;

	for (size_t k = 0; k < (size_t)nlayers; k++) {
		R = T.derivative(CalculationType(CMode::DC, k));
		const double& c = conductivity[k];
		R.P *= c; // Must scale by conductivity to get derivatibe w.r.t log(conductivity)
		R.S *= c;
		memcpy(p, R.primary(XCOMP).data(), sz); p += nw;
		memcpy(p, R.secondary(XCOMP).data(), sz); p += nw;
		memcpy(p, R.primary(YCOMP).data(), sz); p += nw;
		memcpy(p, R.secondary(YCOMP).data(), sz); p += nw;
		memcpy(p, R.primary(ZCOMP).data(), sz); p += nw;
		memcpy(p, R.secondary(ZCOMP).data(), sz); p += nw;
	}
}

/*
void derivative_rx_pitch(void* hS, int n, double rx_pitch, double* xb, double* zb, double* dxbdp, double* dzbdp) {
	TDEmSystem& T = *(TDEmSystem*)hS;
	size_t sz = sizeof(double) * n;

	std::vector<double> dxbdpvec(n);
	std::vector<double> dzbdpvec(n);
	std::vector<double> xbvec(xb, xb + n);
	std::vector<double> zbvec(zb, zb + n);

	//T.drx_pitch(xbvec, zbvec, rx_pitch, dxbdpvec, dzbdpvec);
	T.drx_pitch_new(xbvec, zbvec, rx_pitch, dxbdpvec, dzbdpvec);

	memcpy(dxbdp, dxbdpvec.data(), sz);
	memcpy(dzbdp, dzbdpvec.data(), sz);
}*/



