#ifndef gatdaem1d_H
#define gatdaem1d_H

#ifndef EXPORTED_FUNCTION
	#define EXPORTED_FUNCTION
#endif

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

#ifndef _sTDEmResponseML_
#define _sTDEmResponseML_
struct sTDEmResponseML{  	
	double  PX;	
	double  PY;	
	double  PZ;	
	double* SX;
	double* SY;
	double* SZ;
};
#endif

EXPORTED_FUNCTION void* createhandle(const char* systemfile);
EXPORTED_FUNCTION void deletehandle(void* hS);
EXPORTED_FUNCTION int  nsamplesperwaveform(void* hS);
EXPORTED_FUNCTION void waveform(void* hS, double* time, double* currentwaveform, double* voltagewaveform);
EXPORTED_FUNCTION int  nwindows(void* hS);
EXPORTED_FUNCTION int  nlayers(void* hS);
EXPORTED_FUNCTION void windowtimes(void* hS, double* low, double* high);
EXPORTED_FUNCTION void setgeometry(void* hS, struct sTDEmGeometry G);
EXPORTED_FUNCTION void setearth(void* hS, int nlayers, double* conductivity, double* thickness);
EXPORTED_FUNCTION void forwardmodel(void* hS, struct sTDEmGeometry G, int nlayers, double* conductivity, double* thickness, struct sTDEmResponseML* pR);
EXPORTED_FUNCTION void derivative(void* hS, int dtype, int dlayer, struct sTDEmResponseML* pR);
EXPORTED_FUNCTION void fm_dlogc(void* hS, struct sTDEmGeometry G, int nlayers, double* conductivity, double* thickness, double* R);

#endif
