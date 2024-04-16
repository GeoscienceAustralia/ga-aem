clc;
clear all;

% These examples are setup to be run as if your Matlab current directory is [ga-aem-install-dir]\matlab\examples
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
R = gatdaem1d_forwardmodel(S.hS,G,E);

[f1 A1] = gatdaem1d_amplitudespectrum(S.wfm.time,S.wfm.current);

%Finished modelling so delete the system objects and unload the dll
gatdaem1d_freesystemhandle(S.hS);
gatdaem1d_unloadlibrary();

disp('  Primary');
disp('  X             Z');
disp(num2str([R.PX R.PZ]));

disp('  Secondary');
disp('  X             Z');
disp(num2str([R.SX R.SZ]));

%%Plotting
figure;
hold on;box on;
plot(S.wfm.time,S.wfm.current,'-r.','linewidth',2);
t=S.wfm.time;dt=t(end)-t(1);xl = [t(1)-0.1*dt t(end)+0.1*dt];xlim(xl);ylim([-1.1 1.1]);
for w=1:1:S.nw    
    x=[S.wt.low(w) S.wt.high(w)];
    y=[0.1 0.1] * sign(mod(w,2)-0.5);
    plot(x,y,'g');
end
xlabel('Time (s)');
ylabel('Normalised current (A)');

%%
figure;
i1=find(A1>eps);
loglog(f1(i1),A1(i1),'-b');
hold on;
xlim([10 1e6]);
xlabel('Frequency (Hz)');
ylabel('Amplitude');

%%
figure
hold on;
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'xlim',[9e-6 20e-3]);
xlabel('Time (s)');
ylabel('Response (fT)');
box on;
h1=plot(S.wt.centre, R.SX,'-r.','linewidth',2);
h2=plot(S.wt.centre,-R.SZ,'-b.','linewidth',2);
legend([h1 h2],'X','(-ve)Z');



