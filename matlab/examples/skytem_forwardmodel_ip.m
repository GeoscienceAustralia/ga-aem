clc;
clear all;
clear mex;

%Add path to the gatdaem1d wrapper .m files and the shared library
addpath('..\bin\x64');
addpath('..\gatdaem1d_functions');

%Load the shared library
gatdaem1d_loadlibrary();

%Create a system object, get its handle, and some basic info
%S.stmfile = '..\..\examples\bhmar-skytem\stmfiles\Skytem-LM.stm';
S.stmfile = '..\..\examples\bhmar-skytem\stmfiles\Skytem-HM.stm';
S.hS  = gatdaem1d_getsystemhandle(S.stmfile);
S.nw  = gatdaem1d_nwindows(S.hS);
S.wt  = gatdaem1d_windowtimes(S.hS);
S.wfm = gatdaem1d_waveform(S.hS);
   
%Compute responses (put in loop change G and E as required)
for k=1:1:1
    %Setup geometry
    G.tx_height = 30;
    G.tx_roll   = 0;       G.tx_pitch  = 0; G.tx_yaw    = 0;
    G.txrx_dx   = -12.62;  G.txrx_dy   = 0; G.txrx_dz   = +2.16;
    G.rx_roll   = 0;       G.rx_pitch  = 0; G.rx_yaw    = 0;    
    
    %Setup earth
    E1.thickness           = [20     20];        
    E1.conductivity        = [0.01   0.1    0.001];        
    
    E2 = E1;
    E2.chargeability       = [0.0    0.1   0.0]; 
    E2.timeconstant        = [0.0    0.0005  0.0];  
    E2.frequencydependence = [0.0    1      0.0]; 
    
    R1 = gatdaem1d_forwardmodel(S.hS,G,E1);        
    R2 = gatdaem1d_forwardmodel(S.hS,G,E2,'colecole');
    R3 = gatdaem1d_forwardmodel(S.hS,G,E2,'pelton');
end
%Finished modelling so delete the system objects and unload the dll
gatdaem1d_freesystemhandle(S.hS);
gatdaem1d_unloadlibrary();

%% Plotting
%figure();
dark_figure();
maximize_figure();
hold on;
set(gca,'xscale','log');
set(gca,'yscale','log');
set(gca,'xlim',[min(S.wt.centre)/1.1 1.1*max(S.wt.centre)]);
%set(gca,'ylim',[1e-15 1e-7]);
xlabel('Time (s)');
ylabel('Response (V/A.m^4)');
box on;
h1=plot(S.wt.centre,-R1.SZ,'-w.','linewidth',2);
h2=plot(S.wt.centre,-R2.SZ,'-b.','linewidth',2);
plot(S.wt.centre,R2.SZ,'-b+','linewidth',2);
h3=plot(S.wt.centre,-R3.SZ,'-r.','linewidth',2);
plot(S.wt.centre,R3.SZ,'-r+','linewidth',2);
legend([h1 h2 h3],'Z','Z (pelton-IP)', 'Z (colecole-IP)');


