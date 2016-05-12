import os.path;
from ctypes import *;

#Function to load the shared library
def load_library():
    libname = os.path.dirname(os.path.realpath(__file__)) + "/gatdaem1d.so";
    print "Loading shared library ",libname;
    lib = CDLL(libname)
    return lib;

#Load the shared library
tdlib = load_library();

class WaveForm:
    """Current WaveForm class"""
    def __init__(self, handle):
        self.ns = tdlib.nsamplesperwaveform(handle);
        self.time    = (c_double * self.ns)();
        self.current = (c_double * self.ns)();
        tdlib.waveform(handle,byref(self.time),byref(self.current));

class Geometry(Structure):
    _fields_ = [("tx_height", c_double),
                ("tx_roll", c_double),
                ("tx_pitch", c_double),
                ("tx_yaw", c_double),
                ("txrx_dx", c_double),
                ("txrx_dy", c_double),
                ("txrx_dz", c_double),
                ("rx_roll", c_double),
                ("rx_pitch", c_double),
                ("rx_yaw", c_double)];

class Response(Structure):
    _fields_ = [("PX", c_double),
                ("PY", c_double),
                ("PZ", c_double),
                ("SX", POINTER(c_double)),
                ("SY", POINTER(c_double)),
                ("SZ", POINTER(c_double))];


#void* createhandle(const char* systemfile);
#void deletehandle(void* hS);

#int  nsamplesperwaveform(void* hS);
#void waveform(void* hS, double* time, double* currentwaveform, double* voltagewaveform);

#int  nwindows(void* hS);
#void windowtimes(void* hS, double* low, double* high);

#void setgeometry(void* hS, struct sTDEmGeometry G);
#void setearth(void* hS, int nlayers, double* conductivity, double* thickness);
#void forwardmodel(void* hS, struct sTDEmGeometry G, int nlayers, double* conductivity, double* thickness, struct sTDEmResponseML* pR);
#void derivative(void* hS, int dtype, int dlayer, struct sTDEmResponseML* pR);
#void fm_dlogc(void* hS, struct sTDEmGeometry G, int nlayers, double* conductivity, double* thickness, double* R);
#int nlayers(void* hS);
