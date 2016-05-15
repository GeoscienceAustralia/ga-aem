#from ctypes import *;
import os;

from gatdaem1d import tdlib;
from gatdaem1d import TDAEMSystem;
from gatdaem1d import Earth;
from gatdaem1d import Geometry;
from gatdaem1d import Response;

#Construct the AEM system class instance
#stmfile = "../examples/bhmar-skytem/stmfiles/Skytem-LM.stm";
stmfile  = "..\\examples\\bhmar-skytem\\stmfiles\\Skytem-LM.stm";

S = TDAEMSystem(stmfile);
S.windows.print();
#S.plot_waveform_windows();

#Set the conductivity and thicknesses
conductivity = [0.01,0.1,0.001];
thickness    = [40, 20];
E = Earth(conductivity,thickness);
E.print();

#Set the system geometry
G = Geometry(tx_height=35, txrx_dx=-12.62, txrx_dz = +2.16);
G.print();

#Run a forward model
R = Response(S.windows.nwindows);
S.forwardmodel(G,E,R);
R.print();



