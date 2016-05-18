import os;
import time;
import random;

import numpy;
import matplotlib.pyplot as plt;

from gatdaem1d import tdlib;
from gatdaem1d import TDAEMSystem;
from gatdaem1d import Earth;
from gatdaem1d import Geometry;
from gatdaem1d import Response;

#Construct the AEM system class instance
#stmfile = "../../examples/bhmar-skytem/stmfiles/Skytem-HM.stm";
stmfile  = "..\\..\\examples\\bhmar-skytem\\stmfiles\\Skytem-LM.stm";

S = TDAEMSystem(stmfile);
S.windows.print();

if False:
    fig1 = plt.figure(1);
    S.waveform_windows_plot(fig1);
    plt.show(fig1);

#Set the conductivity and thicknesses
conductivity = [0.01, 0.1, 0.001];
thickness    = [40, 20];
E = Earth(conductivity,thickness);
E.print();

#Set the system geometry
G = Geometry(tx_height=35, txrx_dx = -12.62, txrx_dz = +2.16);
G.print();

#Run a forward model
R = Response(S.windows.nwindows);

t1=time.clock()
for i in range(1):
    conductivity[0] = 10**(random.uniform(-3, 0));
    conductivity[1] = 10**(random.uniform(-3, 0));
    conductivity[2] = 10**(random.uniform(-3, 0));
    thickness       = [40, 20];
    E = Earth(conductivity,thickness);
    S.forwardmodel(G,E,R);
t2=time.clock()
print(t2-t1);

conductivity = [0.05, 0.01, 0.1];
thickness    = [10, 10];
E = Earth(conductivity,thickness);
E.print();

print("\nForward model");
S.forwardmodel(G,E,R); R.print();
fm = numpy.ctypeslib.as_array(R.SZ,shape=(S.windows.nwindows,1)).copy();

print("\nLayer 1 conductivity derivative");
S.derivative(S.CONDUCTIVITYDERIVATIVE,1,R); R.print();
dl1c = numpy.ctypeslib.as_array(R.SZ,shape=(S.windows.nwindows,1)).copy();

print("\nLayer 2 conductivity derivative");
S.derivative(S.CONDUCTIVITYDERIVATIVE,2,R); R.print();
dl2c = numpy.ctypeslib.as_array(R.SZ,shape=(S.windows.nwindows,1)).copy();

print("\nLayer 3 conductivity derivative");
S.derivative(S.CONDUCTIVITYDERIVATIVE,3,R); R.print();
dl3c = numpy.ctypeslib.as_array(R.SZ,shape=(S.windows.nwindows,1)).copy();

print("\nLayer 1 thickness derivative");
S.derivative(S.THICKNESSDERIVATIVE,1,R); R.print();
dl1t = numpy.ctypeslib.as_array(R.SZ,shape=(S.windows.nwindows,1)).copy();

print("\nLayer 2 thickness derivative");
S.derivative(S.THICKNESSDERIVATIVE,2,R); R.print();
dl2t = numpy.ctypeslib.as_array(R.SZ,shape=(S.windows.nwindows,1)).copy();

print("\nHeight derivative");
S.derivative(S.HDERIVATIVE,-1,R); R.print();
print("\nHorizontal radial-distance derivative");
S.derivative(S.RDERIVATIVE,-1,R); R.print();
print("\nHorizontal dx distance derivative");
S.derivative(S.XDERIVATIVE,-1,R); R.print();
print("\nHorizontal dy distance derivative");
S.derivative(S.YDERIVATIVE,-1,R); R.print();
print("\nVertical dz distance derivative");
S.derivative(S.ZDERIVATIVE,-1,R); R.print();

fig2 = plt.figure(2);
ax1 = fig2.add_subplot(1,1,1);
ax1.semilogx(S.windows.centre,-fm,'-k',linewidth=2);
ax1.semilogx(S.windows.centre,-dl1c,'-r',linewidth=2);
ax1.semilogx(S.windows.centre,-dl2c,'-g',linewidth=2);
ax1.semilogx(S.windows.centre,-dl3c,'-b',linewidth=2);
ax1.semilogx(S.windows.centre,-dl1t,'-c',linewidth=2);
ax1.semilogx(S.windows.centre,-dl2t,'-m',linewidth=2);
plt.show(fig2);
quit();

