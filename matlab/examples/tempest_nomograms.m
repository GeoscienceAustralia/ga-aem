clc;
clear all;

%Add path to the tdaem1d_shrlib wrapper .m files and the shared library
addpath('..\bin\x64');
addpath('..\gatdaem1d_functions');

%Load the shared library
gatdaem1d_loadlibrary();

%Create a system object, get its handle, and some basic info
S.stmfile = '..\..\examples\frome-tempest\stmfiles\Tempest-standard.stm';
S.hS  = gatdaem1d_getsystemhandle(S.stmfile);
S.nw  = gatdaem1d_nwindows(S.hS);
S.wt  = gatdaem1d_windowtimes(S.hS);
S.wfm = gatdaem1d_waveform(S.hS);

%Setup geometry
G.tx_height = 120;	
G.tx_roll   = 0;       G.tx_pitch  = 0; G.tx_yaw    = 0;
G.txrx_dx   = -120;    G.txrx_dy   = 0; G.txrx_dz   = -40;			
G.rx_roll   = 0;       G.rx_pitch  = 0; G.rx_yaw    = 0;

c = logspace(-4,1,50);
%Setup earth
for i=1:1:length(c);    
    E.conductivity = c(i);
    E.thickness    = [];
    %Compute responses (put in loop change G and E as required)
    R(i) = gatdaem1d_forwardmodel(S.hS,G,E);
    X(i,:)=R(i).SX;
    Z(i,:)=R(i).SZ;
end
%Finished modelling so delete the system objects and unload the dll
gatdaem1d_freesystemhandle(S.hS);
gatdaem1d_unloadlibrary();
%%

h1=semilogx(c,X,'-b'); 
hold on;
h2=semilogx(c,-Z,'-r');        
xlabel('Halfspace Conductivity (S/m)');
ylabel('Response (fT)');
legend([h1(1) h2(1)],'X','Z');


