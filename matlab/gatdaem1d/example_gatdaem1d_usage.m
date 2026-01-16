clc;
clear all;

%help gatdaem1d_system_handle
%help gatdaem1d_system_handle.forwardmodel
%help gatdaem1d_system_handle.waveform

%libname = 'gatdaem1d'; if libisloaded(libname) ; disp("Unloading library"); unloadlibrary(libname); end; rehash;
%return;

dllpath = 'C:\Users\rossc\AppData\Local\GA-AEM-DEV\bin\gatdaem1d.dll';
dllpath = 'C:\Users\rossc\Work\code\repos\ga-aem\matlab\bin\gatdaem1d.dll';

%stmfile = 'C:\Users\rossc\AppData\Local\GA-AEM-DEV\examples\SkyTEM-BHMAR-2009\stmfiles\Skytem-LM.stm';
stmfile = 'C:\Users\rossc\AppData\Local\GA-AEM-DEV\examples\SkyTEM-BHMAR-2009\stmfiles\Skytem-HM.stm';
S  = gatdaem1d_system_handle(stmfile, "DllPath", dllpath);
%S  = gatdaem1d_system_handle(stmfile);
disp(sprintf('Library ABI version: %d', S.LibraryABIVersion));
disp(sprintf('Proto ABI version: %d', S.ProtoABIVersion));

nw = S.nwindows();
bf = S.basefrequency();
wt = S.windowtimes();

[wt_low,wt_high] = S.windowtimes();


[time, currentwaveform, voltagewaveform] = waveform(S);
%dark_figure(); plot(time,currentwaveform);

G.tx_height =     30;
G.tx_roll   =      0;  G.tx_pitch =  0; G.tx_yaw  = 0;
G.txrx_dx   = -12.62;  G.txrx_dy  = -5; G.txrx_dz = +2.16;
G.rx_roll   =      0;  G.rx_pitch =  0; G.rx_yaw  = 0;

% We only need to get these enum ints once (and same for all systems)
%Either via function call with stringstring
iptype_none     = S.iptype_from_string('NONE');
iptype_colecole = S.iptype_from_string('COLECOLE');
iptype_pelton   = S.iptype_from_string('PELTON');

% Or via the read-only enum values
iptype_none     = S.Enums.IPTYPE.GATDAEM1D_IPTYPE_NONE;
iptype_colecole = S.Enums.IPTYPE.GATDAEM1D_IPTYPE_COLECOLE;
iptype_pelton   = S.Enums.IPTYPE.GATDAEM1D_IPTYPE_PELTON;

%E.iptype = iptype_pelton;
%E.iptype = iptype_colecole;
E.iptype = iptype_none;
E.thickness           = [20   20];
E.conductivity        = [0.01   0.1 0.001];
E.frequencydependence = [0.0    0.5   0.0];

%E.iptype = iptype_colecole;
%E.iptype = iptype_pelton;
%E.chargeability       = [0.0    0.3   0.0];
%E.timeconstant        = [0.0  0.001   0.0];
%E.frequencydependence = [0.0    0.5   0.0];

R = S.forwardmodel(G,E);
DC1 = S.derivative(G, E, 'DC', 0);
DC2 = S.derivative(G, E, 'DC', 1);
DC3 = S.derivative(G, E, 'DC', 2);
DT1 = S.derivative(G, E, 'DT', 0);
DT2 = S.derivative(G, E, 'DT', 1);
A = S.fm_dlogc(G, E);
A;

%%
if(false)
    wt = (wt_low + wt_high)/2;
    figh = dark_figure();
    maximize_figure(figh);
    plot(wt,R.SX,'-ro'); hold on;
    plot(wt,R.SY,'-go');
    plot(wt,-R.SZ,'-bo');
    set(gca,'xscale','log');    
    set(gca,'xlim',[min(wt)*0.9 max(wt)/0.9]);
    
    set(gca,'yscale','log');
    %set(gca,'yscale','linear');
    %set(gca,'ylim',1e-12*[-0.1 1]);
end

%%
S.free_system_handle();
disp('Done...');
sound(1);

