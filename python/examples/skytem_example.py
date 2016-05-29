import time;
import random;
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

#Print the window times
S.windows.print();

if False:
    #Plot the waveform and window positions
    fig1 = plt.figure(1);
    #S.waveform.print(); #Too much printing
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

#Run a few random forward models in a loop
t1=time.clock()
for i in range(10):
    conductivity[0] = 10**(random.uniform(-3, 0));
    conductivity[1] = 10**(random.uniform(-3, 0));
    conductivity[2] = 10**(random.uniform(-3, 0));
    thickness       = [40, 20];
    E = Earth(conductivity,thickness);
    fm = S.forwardmodel(G,E);
t2=time.clock()
print("Time = ",t2-t1);

#Set another earth nodel
conductivity = [0.005, 0.2, 0.01];
thickness    = [20, 10];
E = Earth(conductivity,thickness);
print("\nEarth model");
E.print();

#Do a forward model
print("\nForward model");
fm = S.forwardmodel(G,E); fm.print();

#Then these derivative apply to the last forward model computed
print("\nLayer 1 conductivity derivative");
dl1c = S.derivative(S.CONDUCTIVITYDERIVATIVE,1); dl1c.print();
print("\nLayer 2 conductivity derivative");
dl2c = S.derivative(S.CONDUCTIVITYDERIVATIVE,2); dl2c.print();
print("\nLayer 3 conductivity derivative");
dl3c = S.derivative(S.CONDUCTIVITYDERIVATIVE,3); dl3c.print();
print("\nLayer 1 thickness derivative");
dl1t = S.derivative(S.THICKNESSDERIVATIVE,1); dl1t.print();
print("\nLayer 2 thickness derivative");
dl2t=S.derivative(S.THICKNESSDERIVATIVE,2); dl2t.print();

print("\nHeight derivative");
ddh=S.derivative(S.HDERIVATIVE,-1); ddh.print();
print("\nHorizontal radial-distance derivative");
ddr=S.derivative(S.RDERIVATIVE,-1); ddr.print();
print("\nHorizontal dx distance derivative");
ddx=S.derivative(S.XDERIVATIVE,-1); ddx.print();
print("\nHorizontal dy distance derivative");
ddy=S.derivative(S.YDERIVATIVE,-1); ddy.print();
print("\nVertical dz distance derivative");
ddz=S.derivative(S.ZDERIVATIVE,-1); ddz.print();


if False:
    #Plot the responses
    fig2 = plt.figure(2);
    ax1 = plt.subplot2grid((2,2), (0,0), rowspan=2);
    ax1.loglog(S.windows.centre,-fm.SZ,'-k',linewidth=2,label='Forward model');
    ax1.legend(fontsize=10);

    ax2 = fig2.add_subplot(2,2,2);
    ax2.semilogx(S.windows.centre,-dl1c.SZ,'-r',linewidth=2,label='dL1C');
    ax2.semilogx(S.windows.centre,-dl2c.SZ,'-g',linewidth=2,label='dL2C');
    ax2.semilogx(S.windows.centre,-dl3c.SZ,'-b',linewidth=2,label='dL3C');
    ax2.legend(fontsize=10);

    ax3 = fig2.add_subplot(2,2,4);
    ax3.semilogx(S.windows.centre,-dl1t.SZ,'-c',linewidth=2,label='dL1T');
    ax3.semilogx(S.windows.centre,-dl2t.SZ,'-m',linewidth=2,label='dL2T');
    ax3.legend(fontsize=10);

    plt.show(fig2);

quit();

