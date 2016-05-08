clc;
clear all;

%Add path to the gatdaem1d wrapper .m files and the shared library
addpath('..\bin\x64');
addpath('..\gatdaem1d_functions');

%Load the shared library
gatdaem1d_loadlibrary();

%Setup earth
E.conductivity = [0.05 0.1 0.001];
E.thickness    = [40   20];
%E.conductivity = [0.001];
%E.thickness    = [];


%get system handle
stmfile = '..\..\examples\thomson-vtem\stmfiles\VTEM-plus-7.3ms-pulse-southernthomson.stm';    
S.hS  = gatdaem1d_getsystemhandle(stmfile);
S.nw  = gatdaem1d_nwindows(S.hS);
S.wt  = gatdaem1d_windowtimes(S.hS);
S.wfm = gatdaem1d_waveform(S.hS);

%Setup geometry
G.tx_height = 30;	
G.tx_roll   = 0;    G.tx_pitch  = 0; G.tx_yaw    = 0;
G.txrx_dx   = 0;    G.txrx_dy   = 0; G.txrx_dz   = 0;			
G.rx_roll   = 0;    G.rx_pitch  = 0; G.rx_yaw    = 0;

R = gatdaem1d_forwardmodel(S.hS,G,E);

%Free the handle and unload the shared library
gatdaem1d_freesystemhandle(S.hS);
gatdaem1d_unloadlibrary();

disp('  Secondary');
disp('  X             Z');
disp(num2str([R.SX R.SZ]));

%%Plotting
figure;
hold on;
plot(S.wfm.time,S.wfm.current,'-r','linewidth',2);
t=S.wfm.time;dt=t(end)-t(1);xl = [t(1)-0.1*dt t(end)+0.1*dt];xlim(xl);ylim([-1.1 1.1]);
xlabel('Time (s)');
ylabel('Normalised current (A)');

figure
hold on;
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'xlim',[9e-6 20e-3]);
xlabel('Time (s)');
ylabel('Response (pV/Am^4)');
box on;
h1=plot(S.wt.centre,-R.SZ,'b','linewidth',2);
legend([h1],'(-vw) Z dB/dt');





