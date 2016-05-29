#ifndef gatdaem1d_H
#define gatdaem1d_H

#ifndef EXPORTED_FUNCTION
	#define EXPORTED_FUNCTION
#endif

EXPORTED_FUNCTION void* createhandle(const char* systemfile);
EXPORTED_FUNCTION void deletehandle(void* hS);
EXPORTED_FUNCTION int  nsamplesperwaveform(void* hS);
EXPORTED_FUNCTION void waveform(void* hS, double* time, double* currentwaveform, double* voltagewaveform);
EXPORTED_FUNCTION int  nwindows(void* hS);
EXPORTED_FUNCTION int  nlayers(void* hS);
EXPORTED_FUNCTION void windowtimes(void* hS, double* low, double* high);
EXPORTED_FUNCTION void setgeometry(void* hS, const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw);
EXPORTED_FUNCTION void setearth(void* hS, int nlayers, double* conductivity, double* thickness);
EXPORTED_FUNCTION void forwardmodel(void* hS, const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw, const int nlayers, const double* conductivity, const double* thickness, double* PX, double* PY, double* PZ, double* SX, double* SY, double* SZ);
EXPORTED_FUNCTION void derivative(void* hS, int dtype, int dlayer, double* PX, double* PY, double* PZ, double* SX, double* SY, double* SZ);
EXPORTED_FUNCTION void fm_dlogc(void* hS, const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw, const int nlayers, const double* conductivity, const double* thickness, double* R);

#endif
