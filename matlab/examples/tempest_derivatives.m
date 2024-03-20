clc;
clear all;

% Add path to the gatdaem1d wrapper .m files and the shared library
% addpath('..\bin'); % Not needed if already in your Matlab path
% addpath('..\gatdaem1d_functions'); % Not needed if already in your Matlab path
% addpath('C:\fftw-3.3.5-dll64'); % Not needed if already in your Windows path

%Load the shared library
gatdaem1d_loadlibrary();

%Create a system object, get its handle, and some basic info
S.stmfile = '..\..\examples\Tempest-Frome-2010\stmfiles\Tempest-standard.stm';

S.hS  = gatdaem1d_getsystemhandle(S.stmfile);
S.nw  = gatdaem1d_nwindows(S.hS);
S.wt  = gatdaem1d_windowtimes(S.hS);
S.wfm = gatdaem1d_waveform(S.hS);

%Setup geometry
G.tx_height = 120;	
G.tx_roll   = 0;       G.tx_pitch  = 0; G.tx_yaw    = 0;
G.txrx_dx   = -120;    G.txrx_dy   = 0; G.txrx_dz   = -40;			
G.rx_roll   = 0;       G.rx_pitch  = 0; G.rx_yaw    = 0;

%Setup earth
E.conductivity = [0.05 0.1 0.001];
E.thickness    = [40   20];
%E.conductivity = [0.001];
%E.thickness    = [];


%Compute responses (put in loop change G and E as required)
dtype  = gatdaem1d_derivativestruct();    
dlayer = -1;

delta=0.01;
R1 = gatdaem1d_forwardmodel(S.hS,G,E);

%D  = gatdaem1d_derivative(S.hS,dtype.dH,dlayer);
%G.tx_height = G.tx_height+delta;

%D = gatdaem1d_derivative(S.hS,dtype.dX,dlayer);
%G.txrx_dx = G.txrx_dx+delta;

D = gatdaem1d_derivative(S.hS,dtype.dZ,dlayer);
G.txrx_dz = G.txrx_dz+delta;

R2 = gatdaem1d_forwardmodel(S.hS,G,E);

%Finished modelling so delete the system objects and unload the dll
gatdaem1d_freesystemhandle(S.hS);
gatdaem1d_unloadlibrary();


%%
figure
hold on;
set(gca,'xscale','log');
%set(gca,'yscale','log');
set(gca,'xlim',[9e-6 20e-3]);
xlabel('Time (s)');
ylabel('Response (fT)');
box on;
h1=plot(S.wt.centre, (R2.SX-R1.SX)/delta,'r+','linewidth',2);
h2=plot(S.wt.centre,-(R2.SZ-R1.SZ)/delta,'b+','linewidth',2);
h3=plot(S.wt.centre, (D.SX),'-r','linewidth',2);
h4=plot(S.wt.centre,-(D.SZ),'-b','linewidth',2);

legend([h1 h2 h3 h4],'X numerical','(-ve)Z numerical','X analytic','(-ve)Z analytic');



