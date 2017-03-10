import os.path;
import matplotlib.pyplot as plt;
import ctypes;
from ctypes import c_char_p;
from ctypes import c_void_p;
from ctypes import c_int;
from ctypes import c_double;
from ctypes import POINTER;
import numpy as np;
cptr = np.ctypeslib.as_ctypes;

#Function to load the shared library
def load_library():
    import platform;
    if(platform.system() == "Windows"):
        ext = '.dll'
    else:
        ext = '.so'
    libname = os.path.join(os.path.dirname(os.path.realpath(__file__)),"gatdaem1d"+ext)
    lib = ctypes.CDLL(libname)
    return lib;

#Load the shared library
tdlib = load_library();

class Earth:
    """Earth Class"""
    def __init__(self, conductivity, thickness):
        assert len(conductivity) == 1 + len(thickness);
        self.conductivity = np.array(conductivity,dtype=np.double,order='C');
        self.thickness    = np.array(thickness,dtype=np.double,order='C');

    def nlayers(self):
        return self.conductivity.size;

    def print(self):
        print("nlayers      =",self.nlayers())
        print("conductivity =",self.conductivity[:]);
        print("thickness    =",self.thickness[:]);

class Geometry:
    def __init__(self,tx_height=0,
                 tx_roll=0,tx_pitch=0,tx_yaw=0,
                 txrx_dx=0,txrx_dy=0,txrx_dz=0,
                 rx_roll=0,rx_pitch=0,rx_yaw=0):
        self.tx_height = tx_height;
        self.tx_roll   = tx_roll;
        self.tx_pitch  = tx_pitch;
        self.tx_yaw    = tx_yaw;
        self.txrx_dx   = txrx_dx;
        self.txrx_dy   = txrx_dy;
        self.txrx_dz   = txrx_dz;
        self.rx_roll   = rx_roll;
        self.rx_pitch  = rx_pitch;
        self.rx_yaw    = rx_yaw;
        
    def print(self):
        print("Tx height ",self.tx_height);
        print("Tx roll   ",self.tx_roll);
        print("Tx pitch  ",self.tx_pitch);
        print("Tx yaw    ",self.tx_yaw);
        print("Tx-Rx dx  ",self.txrx_dx);
        print("Tx-Rx dy  ",self.txrx_dy);
        print("Tx-Rx dz  ",self.txrx_dz);
        print("Rx roll   ",self.rx_roll);
        print("Rx pitch  ",self.rx_pitch);
        print("Rx yaw    ",self.rx_yaw);
        
class Response:
    """Response Class"""

    def __init__(self, nwindows):
        self.PX = np.double(0);
        self.PY = np.double(0);
        self.PZ = np.double(0);
        self.SX = np.zeros(nwindows,dtype=np.double,order='C');
        self.SY = np.zeros(nwindows,dtype=np.double,order='C');
        self.SZ = np.zeros(nwindows,dtype=np.double,order='C');

    def print(self):
        print(" Window       SX               SY               SZ");
        for i in range(0, self.SX.size):
            print('{0:5d} {1:16.6e} {2:16.6e} {3:16.6e}'.format(i+1,self.SX[i],self.SY[i],self.SZ[i]));

class Waveform:
    """Waveform Class"""
    def __init__(self, handle):
        n = tdlib.nsamplesperwaveform(handle);
        self.time      = np.zeros(n,dtype=np.double,order='C');
        self.current   = np.zeros(n,dtype=np.double,order='C');
        self.voltage   = np.zeros(n,dtype=np.double,order='C');
        tdlib.waveform(handle,cptr(self.time),cptr(self.current),cptr(self.voltage));

    def nsamples(self):
        return self.time.size;

    def print(self):
        print("Number of waveform samples = ",self.nsamples());
        for i in range(0, self.nsamples()):
            print('{0:5d} {1:10.8f} {2:10.8f}'.format(i+1,self.time[i],self.current[i]));

class Windows:
    """Receiver Windows Class"""
    def __init__(self, handle):
        n = tdlib.nwindows(handle);
        self.low  = np.zeros(n,dtype=np.double,order='C');
        self.high = np.zeros(n,dtype=np.double,order='C');
        tdlib.windowtimes(handle,cptr(self.low),cptr(self.high));
        self.centre = (self.low+self.high)/2.0;

    def nwindows(self):
        return self.low.size;

    def print(self):
        print("Number of windows = ",self.nwindows());
        for i in range(0, self.nwindows()):
            print("{0:5d} {1:10.8f} {2:10.8f}".format(i+1,self.low[i],self.high[i]));

###########################
            
#void* createhandle(const char* systemfile);
tdlib.createhandle.argtypes = [c_char_p];
tdlib.createhandle.restype  = c_void_p;

#void deletehandle(void* hS);
tdlib.deletehandle.argtypes = [c_void_p];
tdlib.deletehandle.restype  = None;

#int nsamplesperwaveform(void* hS);
tdlib.nsamplesperwaveform.argtypes = [c_void_p];
tdlib.nsamplesperwaveform.restype  = c_int;

#void waveform(void* hS, double* time, double* currentwaveform, double* voltagewaveform);
tdlib.waveform.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double), POINTER(c_double)];
tdlib.waveform.restype  = None;

#int nwindows(void* hS);
tdlib.nwindows.argtypes = [c_void_p];
tdlib.nwindows.restype  = c_int;

#void windowtimes(void* hS, double* low, double* high);
tdlib.windowtimes.argtypes = [c_void_p, POINTER(c_double), POINTER(c_double)];
tdlib.windowtimes.restype  = None;

#void forwardmodel(void* hS,
#	const double tx_height, const double tx_roll, const double tx_pitch, const double tx_yaw, const double txrx_dx, const double txrx_dy, const double txrx_dz, const double rx_roll, const double rx_pitch, const double rx_yaw,
#	const int nlayers, const double* conductivity, const double* thickness,
#	double* PX, double* PY, double* PZ, double* SX, double* SY, double* SZ)
#{
tdlib.forwardmodel.argtypes = [c_void_p,c_double,c_double,c_double,c_double,c_double,c_double,c_double,c_double,c_double,c_double,c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double),POINTER(c_double)];
tdlib.forwardmodel.restype  = None;

#void derivative(void* hS, int dtype, int dlayer, double* PX, double* PY, double* PZ, double* SX, double* SY, double* SZ)
tdlib.derivative.argtypes = [c_void_p, c_int, c_int, POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double), POINTER(c_double)];
tdlib.derivative.restype  = None;

class TDAEMSystem:
    """TDAEMSystem Class"""
    
    def __init__(self, stmfile):
        """Initialise the class with the STM file"""
        self.handle   = tdlib.createhandle(stmfile.encode("ascii"));
        self.waveform = Waveform(self.handle);
        self.windows  = Windows(self.handle);

        #Following are defioned in le.h #enum eCalculationType { CT_FORWARDMODEL, CT_CONDUCTIVITYDERIVATIVE, CT_THICKNESSDERIVATIVE, CT_HDERIVATIVE, CT_RDERIVATIVE, CT_XDERIVATIVE,	CT_YDERIVATIVE, CT_ZDERIVATIVE };
        self.FORWARDMODEL=0;
        self.CONDUCTIVITYDERIVATIVE=1;
        self.THICKNESSDERIVATIVE=2;
        self.HDERIVATIVE=3;
        self.RDERIVATIVE=4;
        self.XDERIVATIVE=5;
        self.YDERIVATIVE=6;
        self.ZDERIVATIVE=7;

    def __del__(self):
        """Delete the handle to free up internal resources"""
        tdlib.deletehandle(self.handle);
        
    def waveform_windows_plot(self,fig):
        ax1 = fig.add_subplot(1,1,1);
        ax1.plot(self.waveform.time,self.waveform.current,'-k');
        for i in range(0, self.windows.nwindows(), 2):
            x=[self.windows.low[i],self.windows.low[i],self.windows.high[i],self.windows.high[i]];
            y=[0,0.1,0.1,0];
            ax1.plot(x,y,'-r');
        for i in range(1, self.windows.nwindows(), 2):
            x=[self.windows.low[i],self.windows.low[i],self.windows.high[i],self.windows.high[i]];
            y=[0,-0.1,-0.1,0];
            ax1.plot(x,y,'-b');
        plt.ylabel('Time (s)');
        plt.ylabel('Normalized Current (A)');

    def nwindows(self):
        return self.windows.nwindows();

    def forwardmodel(self,G,E):
        R = Response(self.nwindows());
        tdlib.forwardmodel(self.handle,G.tx_height,G.tx_roll,G.tx_pitch,G.tx_yaw,G.txrx_dx,G.txrx_dy,G.txrx_dz,G.rx_roll,G.rx_pitch,G.rx_yaw,
        E.nlayers(),cptr(E.conductivity),cptr(E.thickness),cptr(R.PX),cptr(R.PY),cptr(R.PZ),cptr(R.SX),cptr(R.SY),cptr(R.SZ));
        return R;

    def derivative(self,dtype,dlayer):
        R = Response(self.nwindows());
        tdlib.derivative(self.handle,dtype,dlayer,cptr(R.PX),cptr(R.PY),cptr(R.PZ),cptr(R.SX),cptr(R.SY),cptr(R.SZ));
        return R;

