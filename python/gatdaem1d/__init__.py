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
    files = os.listdir(os.path.dirname(os.path.realpath(__file__)))
    libname = [file for file in files if 'gatdaem1d' in file][0]
    libname = os.path.join(os.path.dirname(os.path.realpath(__file__)),libname)
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
        self._PX = np.zeros(1,dtype=np.double,order='C');
        self._PY = np.zeros(1,dtype=np.double,order='C');
        self._PZ = np.zeros(1,dtype=np.double,order='C');
        self.SX = np.zeros(nwindows,dtype=np.double,order='C');
        self.SY = np.zeros(nwindows,dtype=np.double,order='C');
        self.SZ = np.zeros(nwindows,dtype=np.double,order='C');

    @property
    def PX(self):
        return self._PX[0]

    @property
    def PY(self):
        return self._PY[0]

    @property
    def PZ(self):
        return self._PZ[0]

    def print(self):
        print(" Primary       PX               PY               PZ");
        print('         {0:16.6e} {1:16.6e} {2:16.6e}'.format(self.PX,self.PY,self.PZ));

        print(" Windows       SX               SY               SZ");
        for i in range(0, self.SX.size):
            print('{0:8d} {1:16.6e} {2:16.6e} {3:16.6e}'.format(i+1,self.SX[i],self.SY[i],self.SZ[i]));

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

#int nturns(void* hS);
tdlib.nturns.argtypes = [c_void_p];
tdlib.nturns.restype  = c_int;

#double looparea(void* hS);
tdlib.looparea.argtypes = [c_void_p];
tdlib.looparea.restype  = c_double;

#double basefrequency(void* hS);
tdlib.basefrequency.argtypes = [c_void_p];
tdlib.basefrequency.restype  = c_double;

#double peakcurrent(void* hS);
tdlib.peakcurrent.argtypes = [c_void_p];
tdlib.peakcurrent.restype  = c_double;

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

tdlib.fm_dlogc.argtypes = [c_void_p,c_double,
                           c_double,c_double,c_double,
                           c_double,c_double,c_double,
                           c_double,c_double,c_double,
                           c_int, POINTER(c_double), POINTER(c_double),
                           POINTER(c_double)];
tdlib.fm_dlogc.restype  = None;
class TDAEMSystem:
    """TDAEMSystem Class"""

    def __init__(self, stmfile):
        """Initialise the class with the STM file"""
        self.handle   = tdlib.createhandle(stmfile.encode("ascii"));
        self.waveform = Waveform(self.handle);
        self.windows  = Windows(self.handle);

        #Following are defined in le.h #enum eCalculationType { CT_FORWARDMODEL, CT_CONDUCTIVITYDERIVATIVE, CT_THICKNESSDERIVATIVE, CT_HDERIVATIVE, CT_RDERIVATIVE, CT_XDERIVATIVE,	CT_YDERIVATIVE, CT_ZDERIVATIVE };
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

    def nTurns(self):
        return tdlib.nturns(self.handle);

    def loopRadius(self):
        return np.sqrt(self.loopArea() / np.pi)

    def peakCurrent(self):
        return tdlib.peakcurrent(self.handle)

    def baseFrequency(self):
        return tdlib.basefrequency(self.handle)

    def loopArea(self):
        return tdlib.looparea(self.handle);

    def forwardmodel(self, G, E):
        R = Response(self.nwindows());
        tdlib.forwardmodel(self.handle,G.tx_height,
                           G.tx_roll,G.tx_pitch,G.tx_yaw,
                           G.txrx_dx,G.txrx_dy,G.txrx_dz,
                           G.rx_roll,G.rx_pitch,G.rx_yaw,
                           E.nlayers(),cptr(E.conductivity),cptr(E.thickness),
                           cptr(R._PX),cptr(R._PY),cptr(R._PZ),
                           cptr(R.SX),cptr(R.SY),cptr(R.SZ));
        return R;

    def derivative(self, dtype, dlayer):
        R = Response(self.nwindows());
        tdlib.derivative(self.handle, dtype, dlayer,
                         cptr(R._PX),cptr(R._PY),cptr(R._PZ),
                         cptr(R.SX),cptr(R.SY),cptr(R.SZ));
        return R;

    def fm_dlogc(self, G, E):

        nchan = (1+self.nwindows());
        ncomp = 3;
        ncalc = (1+E.nlayers());
        len = nchan * ncomp * ncalc;

        tmp = np.zeros(len, dtype=np.double, order='C')


        tdlib.fm_dlogc(self.handle, G.tx_height,
                        G.tx_roll,G.tx_pitch,G.tx_yaw,
                        G.txrx_dx,G.txrx_dy,G.txrx_dz,
                        G.rx_roll,G.rx_pitch,G.rx_yaw,
                        E.nlayers(),cptr(E.conductivity),cptr(E.thickness),
                        cptr(tmp)
                       )

        tmp = np.reshape(tmp,(ncalc, ncomp, nchan));

        R = Response(self.nwindows())

        R._PX[0] = tmp[0, 0, 0]
        R._PY[0] = tmp[0, 1, 0]
        R._PZ[0] = tmp[0, 2, 0]

        R._SX = tmp[0, 0, 1:]
        R._SY = tmp[0, 1, 1:]
        R._SZ = tmp[0, 2, 1:]

        Jx = tmp[1:, 0, 1:]
        Jy = tmp[1:, 1, 1:]
        Jz = tmp[1:, 2, 1:]

        return R, Jx, Jy, Jz