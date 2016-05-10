from ctypes import *;

from gatdaem1d import tdlib;
from gatdaem1d import WaveForm;
from gatdaem1d import Geometry;
from gatdaem1d import Response;


#Get the AEM system's handle
hs  = tdlib.createhandle("/short/public/rcb547/apps/ga-aem/examples/bhmar-skytem/stmfiles/Skytem-LM.stm");
wfm = WaveForm(hs)
print wfm.__dict__


#Get the number of receiver windows and the start/end times
nw     = tdlib.nwindows(hs);
wtlow  = (c_double * nw)();
wthigh = (c_double * nw)();
tdlib.windowtimes(hs,byref(wtlow),byref(wthigh));

#print "Number of waveform samples = %d" % ns;
#for i in range(0, ns):
#	print '{0:5d} {1:10.8f} {2:10.8f}'.format(i+1,wftime[i],wfcurrent[i]);

#print "Number of windows =",nw;
#for i in range(0, nw):
#    print '{0:5d} {1:10.8f} {2:10.8f}'.format(i+1,wtlow[i],wthigh[i]);

G = Geometry();
G.tx_height = 30;
G.tx_roll   = 0;       G.tx_pitch  = 0; G.tx_yaw    = 0;
G.txrx_dx   = -12.62;  G.txrx_dy   = 0; G.txrx_dz   = +2.16;
G.rx_roll   = 0;       G.rx_pitch  = 0; G.rx_yaw    = 0;    

#Set the conductivity and thicknesses
nlayers = 3;
conductivity    = (c_double * nlayers)();
thickness       = (c_double * (nlayers-1))();
conductivity[0] = 0.01;
conductivity[1] = 0.1;
conductivity[2] = 0.001;
thickness[0]    = 40;
thickness[1]    = 20;

#Set up the response struct
R = Response();
R.SX = (c_double * nw)();
R.SY = (c_double * nw)();
R.SZ = (c_double * nw)();

#Run a forward model
tdlib.forwardmodel(hs,G,nlayers,byref(conductivity),byref(thickness),pointer(R));

#After finishing modelling delete the system handle to free up resources
tdlib.deletehandle(hs);

print "geometry     =",G.tx_height, G.tx_roll, G.tx_pitch, G.tx_yaw, G.txrx_dx, G.txrx_dy, G.txrx_dz, G.rx_roll, G.rx_pitch, G.rx_yaw
print "nlayers      =",nlayers
print "conductivity =",conductivity[:]
print "thickness    =",thickness[:]
print " Window       SX               SY               SZ"
for i in range(0, nw):
    print '{0:5d} {1:16.6e} {2:16.6e} {3:16.6e}'.format(i+1,R.SX[i],R.SY[i],R.SZ[i]);

