clc;
clear all;


%Add path to the gatdaem1d wrapper .m files and the shared library
addpath('..\bin\x64');
addpath('..\gatdaem1d_functions');

%Load the shared library
gatdaem1d_loadlibrary();

%Create a system object, get its handle, and some basic info
S.stmfile = '..\..\examples\frome-tempest\stmfiles\Tempest-standard-diagnostics.stm';

S.hS  = gatdaem1d_getsystemhandle(S.stmfile);
S.nw  = gatdaem1d_nwindows(S.hS);
S.wt  = gatdaem1d_windowtimes(S.hS);
S.wfm = gatdaem1d_waveform(S.hS);

%Setup system geometry
G.tx_height = 120;
G.tx_roll   = 0;       G.tx_pitch  = 0; G.tx_yaw    = 0;
G.txrx_dx   = -120;    G.txrx_dy   = 0; G.txrx_dz   = -40;
G.rx_roll   = 0;       G.rx_pitch  = 0; G.rx_yaw    = 0;


%Setup earth
E.thickness           = [20     20];
E.conductivity        = [0.01   0.1    0.001];
E.chargeability       = [0.0    0.3     0.0];
E.timeconstant        = [0.0    0.001  0.0];
E.frequencydependence = [0.0    0.5    0.0];

%Compute responses
R = gatdaem1d_forwardmodel(S.hS,G,E);
R = gatdaem1d_forwardmodel(S.hS,G,E,'pelton');
%R = gatdaem1d_forwardmodel(S.hS,G,E,'colecole');


%Finished modelling so delete the system objects and unload the dll
gatdaem1d_freesystemhandle(S.hS);
gatdaem1d_unloadlibrary();

%%
%Load the intermediate modelling steps file
c = load('diag_discretefrequencies.txt');
df=c(:,1);
dxr=c(:,2);
dxi=c(:,3);
dyr=c(:,4);
dyi=c(:,5);
dzr=c(:,6);
dzi=c(:,7);

b = load('diag_splinedfrequencies.txt');
sf=b(:,1);
sxr=b(:,2);
sxi=b(:,3);
syr=b(:,4);
syi=b(:,5);
szr=b(:,6);
szi=b(:,7);

d = load('diag_frequencydomainwaveform.txt');
wf =d(:,1);
wfr=d(:,2);
wfi=d(:,3);
wtr=d(:,4);
wti=d(:,5);

xts = load('diag_xtimeseries.txt');
yts = load('diag_ytimeseries.txt');
zts = load('diag_ztimeseries.txt');
t=xts(:,1);
x=xts(:,2);
y=yts(:,2);
z=zts(:,2);

a = load('diag_windows.txt');
wtlo=a(:,2);
wthi=a(:,3);
wt = (wtlo+wthi)/2;
wx=a(:,4);
wy=a(:,5);
wz=a(:,6);

% dark_figure();
% maximize_figure();
% plot(wt,-wc,'w');
% hold on;
% plot(t,-z,'r');
% plot(l,0,'go');
% plot(h,0,'mo');
% grid on;





%% Plotting
dark_figure(1);
maximize_figure();
subplot(1,1,1);
h1=semilogx(df,dzr,'ro');
hold on;
h2=plot(df,dzi,'bo');
h3=plot(sf,szr,'r.');
h4=plot(sf,szi,'b.');

lh=legend([h1 h2 h3 h4],'Z real','Z imag','Z real splined','Z imag splined');
set(lh,'location','northwest');
ylabel('BZ (T)');
xlabel('Frequency (Hz)');
title("Frequency Domain Earth Response");

%%
dark_figure(2);
maximize_figure();
subplot(2,1,1);
semilogx(wf,wfr,'g');
hold on;
semilogx(wf,wfi,'m');
title("Frequency Domain Waveform");

subplot(2,1,2);
semilogx(wf,wtr,'g');
hold on;
semilogx(wf,wti,'m');
xlabel('Frequency (Hz)');
title("Frequency Domain Transfer Function");
%%

dark_figure(3); 
maximize_figure();
subplot(3,1,1);
plot(t,x);hold on;
hold on;
yl=get(gca,'ylim');
for k=1:1:length(wt)
    plot([wtlo(k) wtlo(k)],yl,':g')
    plot([wthi(k) wthi(k)],yl,':y')
end
xlabel('Time (s)');
ylabel('X Response');

subplot(3,1,2);
plot(t,y);hold on;
hold on;
yl=get(gca,'ylim');
for k=1:1:length(wt)
    plot([wtlo(k) wtlo(k)],yl,':g')
    plot([wthi(k) wthi(k)],yl,':y')
end
xlabel('Time (s)');
ylabel('Y Response');

subplot(3,1,3);
plot(t,z*(4*pi*1e-7));hold on;
hold on;
yl=get(gca,'ylim');
for k=1:1:length(wt)
    plot([wtlo(k) wtlo(k)],yl,':g')
    plot([wthi(k) wthi(k)],yl,':y')
end

xlabel('Time (s)');
ylabel('Z Response');

%%
dark_figure(4); 
maximize_figure();
subplot(1,1,1);
h1=loglog(wt,wx,'-ro');hold on;
%loglog(wt,wy);
h2=loglog(wt,-wz,'-bo');
xlabel('Centre Time (s)');
ylabel('Response (fT)');
lh=legend([h1 h2],'X','-Z');




