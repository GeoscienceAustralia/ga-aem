clc;
clear all;
clear mex;

%Add path to the gatdaem1d wrapper .m files and the shared library
%addpath('..\bin\x64');
addpath('..\gatdaem1d_functions');

%Load the shared library
gatdaem1d_loadlibrary();

%Create a LM system object, get its handle, and some basic info
LM.stmfile = '..\..\examples\bhmar-skytem\stmfiles\Skytem-LM.stm';
LM.hS  = gatdaem1d_getsystemhandle(LM.stmfile);
LM.nw  = gatdaem1d_nwindows(LM.hS);
LM.wt  = gatdaem1d_windowtimes(LM.hS);
LM.wfm = gatdaem1d_waveform(LM.hS);

% Additive noise for BHMAR
LM.z_mnp = 3.6;
%LM.z_an  = [5.7761e-011 7.7154e-012 5.7849e-012 3.9164e-012 3.1502e-012 2.5105e-012 2.2912e-012  1.921e-012  1.733e-012  1.529e-012 1.2258e-012 9.6876e-013 9.0323e-013 8.2181e-013 7.4835e-013 6.2648e-013 6.2901e-013 5.7157e-013 5.1475e-013]';
LM.z_an  = [7.7154e-012 5.7849e-012 3.9164e-012 3.1502e-012 2.5105e-012 2.2912e-012  1.921e-012  1.733e-012  1.529e-012 1.2258e-012 9.6876e-013 9.0323e-013 8.2181e-013 7.4835e-013 6.2648e-013 6.2901e-013 5.7157e-013 5.1475e-013]';
HM.z_mnp = 3.6;
HM.z_an  = [2.5545e-013 2.0815e-013 1.9144e-013  1.592e-013 1.4598e-013 1.3402e-013 1.2712e-013 1.0844e-013 1.0214e-013 9.7184e-014 9.0881e-014 8.4579e-014 7.7776e-014 6.9864e-014 6.6747e-014 5.9365e-014   5.33e-014  4.843e-014 4.2199e-014 3.7096e-014  3.571e-014]';

%Create a HM system object, get its handle, and some basic info
HM.stmfile = '..\..\examples\bhmar-skytem\stmfiles\Skytem-HM.stm';
HM.hS  = gatdaem1d_getsystemhandle(HM.stmfile);
HM.nw  = gatdaem1d_nwindows(HM.hS);
HM.wt  = gatdaem1d_windowtimes(HM.hS);
HM.wfm = gatdaem1d_waveform(HM.hS);

filename = '..\..\examples\bhmar-skytem\data\synthetic.dat';
fp=fopen(filename,'w');
flight = 45;
line = 20010;
fid  = 0;
x    = 300000 - 25;
y    = 6200000;
elevation = 200;

N=100;
for k = 0:1:N
    
    x   = x+25;
    fid = fid+1;
    
    %Setup geometry
    G.tx_height = 30;
    G.tx_roll   = 0;       G.tx_pitch  = 0; G.tx_yaw    = 0;
    G.txrx_dx   = -12.62;  G.txrx_dy   = 0; G.txrx_dz   = +2.16;
    G.rx_roll   = 0;       G.rx_pitch  = 0; G.rx_yaw    = 0;    
    
    %Setup earth
    f=k/N;
    a = [0.010 20+f*20
        0.100  1+(1-f)*10
        0.030  50
        0.1    20+(1-f)*10
        0.001  Inf];
    E.conductivity = a(:,1);
    E.thickness    = a(1:size(a,1)-1,2);
            
    LR = gatdaem1d_forwardmodel(LM.hS,G,E);
    HR = gatdaem1d_forwardmodel(HM.hS,G,E);                        
    
    LR.EZ  = get_noise_instance(LR.SZ,LM.z_mnp,LM.z_an);%noise instance
    LR.SZN = LR.SZ + LR.EZ;%noise contaminated data

    HR.EZ  = get_noise_instance(HR.SZ,HM.z_mnp,HM.z_an);%noise instance
    HR.SZN = HR.SZ + HR.EZ;%noise contaminated data
    
    fprintf(fp,'%8d',flight);
    fprintf(fp,'%8d',line);
    fprintf(fp,'%8d',fid);    
    fprintf(fp,'%10.1f',x);
    fprintf(fp,'%10.1f',y);
    fprintf(fp,'%10.1f',elevation);
    
    fprintf(fp,'%8.3f',G.tx_height);
    fprintf(fp,'%8.3f',G.tx_roll);
    fprintf(fp,'%8.3f',G.tx_pitch);
    fprintf(fp,'%8.3f',G.tx_yaw);
    
    fprintf(fp,'%8.3f',G.txrx_dx);
    fprintf(fp,'%8.3f',G.txrx_dy);
    fprintf(fp,'%8.3f',G.txrx_dz);
    
    fprintf(fp,'%8.3f',G.rx_roll);
    fprintf(fp,'%8.3f',G.rx_pitch);
    fprintf(fp,'%8.3f',G.rx_yaw);    
    
    fprintf(fp,'%16.6e',-LR.SZ);
    fprintf(fp,'%16.6e',-LR.SZN);
    fprintf(fp,'%16.6e',-LR.EZ);
    
    fprintf(fp,'%16.6e',-HR.SZ);
    fprintf(fp,'%16.6e',-HR.SZN);
    fprintf(fp,'%16.6e',-HR.EZ);
    
    fprintf(fp,'%8d',length(E.conductivity));
    fprintf(fp,'%16.6e',E.conductivity);
    fprintf(fp,'%16.6e',E.thickness);
    
    fprintf(fp,'\n');
end


%Finished modelling so delete the system objects and unload the dll
gatdaem1d_freesystemhandle(LM.hS);
gatdaem1d_freesystemhandle(HM.hS);
gatdaem1d_unloadlibrary();
fclose(fp);
fclose all;
disp("done");
return

%%

LR.EZ  = get_noise_instance(LR.SZ,LM.z_mnp,LM.z_an);
LR.SZN = LR.SZ + LR.EZ;

HR.EZ  = get_noise_instance(HR.SZ,HM.z_mnp,HM.z_an);
HR.SZN = HR.SZ + HR.EZ;

dark_figure();
maximize_figure();

subplot(1,2,1);
semilogy(-LR.SZ,'-w');
hold on;
semilogy(-LR.SZN,'r+');
errorbar(-LR.SZN,LR.EZ);

subplot(1,2,2);
semilogy(-HR.SZ,'-w');
hold on;
semilogy(-HR.SZN,'r+');
errorbar(-HR.SZN,HR.EZ);


