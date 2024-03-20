clc;
clear all;
clear mex;

%Add path to the gatdaem1d wrapper .m files and the shared library
addpath('..\bin');
addpath('..\gatdaem1d_functions');

%Load the shared library
gatdaem1d_loadlibrary();

%Create a LM system object, get its handle, and some basic info
LM.stmfile = '..\..\examples\SkyTEM-BHMAR-2009\stmfiles\Skytem-LM.stm';
LM.hS  = gatdaem1d_getsystemhandle(LM.stmfile);
LM.nw  = gatdaem1d_nwindows(LM.hS);
LM.wt  = gatdaem1d_windowtimes(LM.hS);
LM.wfm = gatdaem1d_waveform(LM.hS);

%Create a HM system object, get its handle, and some basic info
HM.stmfile = '..\..\examples\SkyTEM-BHMAR-2009\stmfiles\Skytem-HM.stm';
HM.hS  = gatdaem1d_getsystemhandle(HM.stmfile);
HM.nw  = gatdaem1d_nwindows(HM.hS);
HM.wt  = gatdaem1d_windowtimes(HM.hS);
HM.wfm = gatdaem1d_waveform(HM.hS);

dtype = gatdaem1d_derivativestruct();    
dlayer=1;
    

%Compute responses (put in loop change G and E as required)
for k=1:1:10
    %Setup geometry
    G.tx_height = 30;
    G.tx_roll   = 0;       G.tx_pitch  = 0; G.tx_yaw    = 0;
    G.txrx_dx   = -12.62;  G.txrx_dy   = 0; G.txrx_dz   = +2.16;
    G.rx_roll   = 0;       G.rx_pitch  = 0; G.rx_yaw    = 0;    
    
    %Setup earth
    E.conductivity = [5.05 0.1 0.05 0.001];        
    E.thickness    = [20   50  20];        
            
    %LR = gatdaem1d_forwardmodel(LM.hS,G,E);
    %HR = gatdaem1d_forwardmodel(HM.hS,G,E);                        
    %LD = gatdaem1d_derivative(LM.hS,dtype.dC,dlayer);
    %HD = gatdaem1d_derivative(HM.hS,dtype.dC,dlayer);       
    L = gatdaem1d_fm_dlogc(LM.hS,G,E);
    H = gatdaem1d_fm_dlogc(HM.hS,G,E);
end


%Finished modelling so delete the system objects and unload the dll
gatdaem1d_freesystemhandle(LM.hS);
gatdaem1d_freesystemhandle(HM.hS);
gatdaem1d_unloadlibrary();

%%
%[num2str(L.FM.SX,'%18.6e')  num2str(L.FM.SZ,'%18.6e')]
%[num2str(H.FM.SX,'%18.6e')  num2str(H.FM.SZ,'%18.6e')]

%%Plotting
figure
hold on;
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'xlim',[1e-5  1e-2]);
set(gca,'ylim',[1e-15 1e-7]);
xlabel('Time (s)');
ylabel('Response (V/A.m^4)');
box on;
h1=plot(LM.wt.centre,-L.FM.SZ,'-r.','linewidth',2);
h2=plot(HM.wt.centre,-H.FM.SZ,'-b.','linewidth',2);
legend([h1 h2],'LM Z','HM Z');

%%
figure
hold on;
set(gca,'xscale','log');
set(gca,'yscale','linear');
set(gca,'xlim',[9e-6 1e-2]);
xlabel('Time (s)');
ylabel('Response (V/A.m^4)');
box on;
for k=1:1:length(E.conductivity)
    h1=plot(LM.wt.centre,-L.dlogC(k).SZ,':r','linewidth',2);
    h2=plot(HM.wt.centre,-H.dlogC(k).SZ,':b','linewidth',2);        
end
legend([h1 h2],'d(LM Z)/dlogC','d(HM Z)/dC');

%%
figure;
subplot(2,1,1);
hold on;box on;
plot(LM.wfm.time,LM.wfm.current,'-r','linewidth',2);
t=LM.wfm.time;dt=t(end)-t(1);xl = [t(1)-0.1*dt t(end)+0.1*dt];xlim(xl);
for w=1:1:LM.nw    
    x=[LM.wt.low(w) LM.wt.high(w)];
    y=[0.1 0.1] * sign(mod(w,2)-0.5);
    plot(x,y,'g');
end
xlabel('Time (s)');
ylabel('Normalised current (A)');
title('Low momemnt');

subplot(2,1,2);
hold on;box on;
plot(HM.wfm.time,HM.wfm.current,'-b','linewidth',2);
t=HM.wfm.time;dt=t(end)-t(1);xl = [t(1)-0.1*dt t(end)+0.1*dt];xlim(xl);
for w=1:1:HM.nw    
    x=[HM.wt.low(w) HM.wt.high(w)];
    y=[0.1 0.1] * sign(mod(w,2)-0.5);
    plot(x,y,'g');
end
xlabel('Time (s)');
ylabel('Normalised current (A)');
title('High momemnt');

