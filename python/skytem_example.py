#from ctypes import *;
import os;
import time;
import random;
from array import array;

from gatdaem1d import tdlib;
from gatdaem1d import TDAEMSystem;
from gatdaem1d import Earth;
from gatdaem1d import Geometry;
from gatdaem1d import Response;

#Construct the AEM system class instance
#stmfile = "../examples/bhmar-skytem/stmfiles/Skytem-HM.stm";
stmfile  = "..\\examples\\bhmar-skytem\\stmfiles\\Skytem-LM.stm";

S = TDAEMSystem(stmfile);
S.windows.print();
#S.plot_waveform_windows();

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

conductivity = [0.001, 0.1, 0.001];
thickness    = [1, 1];
E = Earth(conductivity,thickness);
E.print();

print("\nForward model");
S.derivative(S.FORWARDMODEL,-1,R); R.print();
print("\nLayer 1 conductivity derivative");
S.derivative(S.CONDUCTIVITYDERIVATIVE,1,R); R.print();
print("\nLayer 2 conductivity derivative");
S.derivative(S.CONDUCTIVITYDERIVATIVE,2,R); R.print();
print("\nLayer 3 conductivity derivative");
S.derivative(S.CONDUCTIVITYDERIVATIVE,3,R); R.print();

print("\nLayer 1 thickness derivative");
S.derivative(S.THICKNESSDERIVATIVE,1,R); R.print();

print("\nLayer 2 thickness derivative");
S.derivative(S.THICKNESSDERIVATIVE,2,R); R.print();

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


