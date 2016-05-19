import os.path;
import ctypes;
from ctypes import c_char_p;
from ctypes import c_void_p;
from ctypes import c_int;
from ctypes import c_double;
from ctypes import Structure;
from ctypes import POINTER;

import numpy as np;
import matplotlib.pyplot as plt;

#Function to load the shared library
def load_library():
    #libname = os.path.dirname(os.path.realpath(__file__)) + "/gatdaem1d.so";
    libname = os.path.dirname(os.path.realpath(__file__)) + "\gatdaem1d.dll";
    print("Loading shared library ",libname);
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
    
    def __init__(self,tx_height=0,
                 tx_roll=0,tx_pitch=0,tx_yaw=0,
                 txrx_dx=0,txrx_dy=0,txrx_dz=0,
                 rx_roll=0,rx_pitch=0,rx_yaw=0):
        self.tx_height = c_double(tx_height);
        self.tx_roll   = c_double(tx_roll);
        self.tx_pitch  = c_double(tx_pitch);
        self.tx_yaw    = c_double(tx_yaw);
        self.txrx_dx   = c_double(txrx_dx);
        self.txrx_dy   = c_double(txrx_dy);
        self.txrx_dz   = c_double(txrx_dz);        
        self.rx_roll   = c_double(rx_roll);
        self.rx_pitch  = c_double(rx_pitch);
        self.rx_yaw    = c_double(rx_yaw);
        
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
        
class Response(Structure):
    """Response Class"""
    _fields_ = [("PX", c_double),
                ("PY", c_double),
                ("PZ", c_double),
                ("pSX", POINTER(c_double)),
                ("pSY", POINTER(c_double)),
                ("pSZ", POINTER(c_double))];

    def __init__(self, nwindows):
        self.PX = c_double(0);
        self.PY = c_double(0);
        self.PZ = c_double(0);
        self.SX = np.zeros(nwindows,dtype=np.double,order='C');
        self.SY = np.zeros(nwindows,dtype=np.double,order='C');
        self.SZ = np.zeros(nwindows,dtype=np.double,order='C');
        self.pSX = np.ctypeslib.as_ctypes(self.SX);
        self.pSY = np.ctypeslib.as_ctypes(self.SY);
        self.pSZ = np.ctypeslib.as_ctypes(self.SZ);

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
        tdlib.waveform(handle,
            np.ctypeslib.as_ctypes(self.time),
            np.ctypeslib.as_ctypes(self.current),
            np.ctypeslib.as_ctypes(self.voltage));

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
        self.low      = np.zeros(n,dtype=np.double,order='C');
        self.high     = np.zeros(n,dtype=np.double,order='C');
        tdlib.windowtimes(handle,
            np.ctypeslib.as_ctypes(self.low),
            np.ctypeslib.as_ctypes(self.high));
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

#void forwardmodel(void* hS, struct sTDEmGeometry G, int nlayers, double* conductivity, double* thickness, struct sTDEmResponseML* pR);
tdlib.forwardmodel.argtypes = [c_void_p, Geometry, c_int, POINTER(c_double), POINTER(c_double), POINTER(Response)];
tdlib.forwardmodel.restype  = None;

#void derivative(void* hS, int dtype, int dlayer, struct sTDEmResponseML* pR);
tdlib.derivative.argtypes = [c_void_p, c_int, c_int, POINTER(Response)];
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
        tdlib.forwardmodel(self.handle,G,E.nlayers(),np.ctypeslib.as_ctypes(E.conductivity),np.ctypeslib.as_ctypes(E.thickness),R);
        return R;

    def derivative(self,dtype,dlayer):
        R = Response(self.nwindows());
        tdlib.derivative(self.handle,dtype,dlayer,R);
        return R;

#void setgeometry(void* hS, struct sTDEmGeometry G);
#void setearth(void* hS, int nlayers, double* conductivity, double* thickness);
#void forwardmodel(void* hS, struct sTDEmGeometry G, int nlayers, double* conductivity, double* thickness, struct sTDEmResponseML* pR);
#void derivative(void* hS, int dtype, int dlayer, struct sTDEmResponseML* pR);
#void fm_dlogc(void* hS, struct sTDEmGeometry G, int nlayers, double* conductivity, double* thickness, double* R);
#int nlayers(void* hS);
