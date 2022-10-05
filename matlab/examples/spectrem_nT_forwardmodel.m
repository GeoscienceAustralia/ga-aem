clc;
clear all;

%Add path to the gatdaem1d wrapper .m files and the shared library
addpath('..\bin');
addpath('..\gatdaem1d_functions');

%Load the shared library
gatdaem1d_loadlibrary();

%Create a system object, get its handle, and some basic info
S.stmfile = '..\..\examples\spectrem\stmfiles\Spectrem-nT.stm';
S.hS  = gatdaem1d_getsystemhandle(S.stmfile);
S.nw  = gatdaem1d_nwindows(S.hS);
S.wt  = gatdaem1d_windowtimes(S.hS);
S.wfm = gatdaem1d_waveform(S.hS);

%Setup geometry
G.tx_height = 90;	
G.tx_roll   = 0;       G.tx_pitch  = 0; G.tx_yaw    = 0;
G.txrx_dx   = -123;    G.txrx_dy   = 0; G.txrx_dz   = -36;			
G.rx_roll   = 0;       G.rx_pitch  = 0; G.rx_yaw    = 0;

%Setup earth
E.conductivity = [0.1  0.4 0.02];
E.thickness    = [10   20];

%Compute responses (put in loop change G and E as required)
R = gatdaem1d_forwardmodel(S.hS,G,E);

%Finished modelling so delete the system objects and unload the dll
gatdaem1d_freesystemhandle(S.hS);
gatdaem1d_unloadlibrary();

disp('  Primary');
disp('  X             Z');
disp(num2str([R.PX R.PZ]));

disp('  Secondary');
disp('  X             Z');
disp(num2str([R.SX R.SZ]));

%Plotting
figure;
hold on;box on;
plot(S.wfm.time,S.wfm.current,'r','linewidth',2);
t=S.wfm.time;dt=t(end)-t(1);xl = [t(1)-0.1*dt t(end)+0.1*dt];xlim(xl);ylim([-1.1 1.1]);
for w=1:1:S.nw    
    x=[S.wt.low(w) S.wt.high(w)];
    y=[0.1 0.1] * sign(mod(w,2)-0.5);
    plot(x,y,'g');
end
xlabel('Time (s)');
ylabel('Normalised current (A)');

%%
figure
hold on;box on;
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'xlim',[9e-6 20e-3]);
xlabel('Time (s)');
ylabel('Response (nT)');
box on;
h1=plot(S.wt.centre, R.SX,'r','linewidth',2);
h2=plot(S.wt.centre, -R.SZ,'b','linewidth',2);
legend([h1 h2],'X','(-ve)Z');





