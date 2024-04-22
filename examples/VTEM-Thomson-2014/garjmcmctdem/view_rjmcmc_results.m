clc;
clear all;
close all;

addpath("..\..\..\matlab\rjmcmctdem_functions");
addpath("..\..\..\matlab\utility_functions");

%%True synthetic model if applicable
%TM.c = [0.0500 0.2000 0.0100];
%TM.t = [40 30];
TM=[];%Otherwise nothing

basedir = 'output\';

ncdir   = [basedir 'pmaps\'];
plotdir = [basedir 'pmap_plots\'];
if(~exist(plotdir','dir'))
    mkdir(plotdir);
end

%Component info for plotting nicely
C(1).name   = 'Z';%name of the component
C(1).dindex = 1:45;%index in the data array for this component
C(1).wtime =  [2.05e-05     2.6e-05    3.15e-05    3.65e-05     4.2e-05     4.8e-05     5.5e-05    6.35e-05     7.3e-05     8.4e-05    9.65e-05   0.0001105    0.000127    0.000146   0.0001675   0.0001925    0.000221   0.0002535   0.0002915    0.000335   0.0003845   0.0004415   0.0005075    0.000583   0.0006695   0.0007695    0.000884   0.0010155   0.0011665   0.0013395    0.001539   0.0017685   0.0020315   0.0023335     0.00268   0.0030785   0.0035365    0.004061    0.004664    0.005358    0.006155   0.0070705   0.0081225    0.009331    0.010717];%window centre times

%Add other components as required
%C(2).name   = ;%name of the component
%C(2).dindex = ;%index in the data array for this component
%C(2).wtime =  ;%window centre times

subsam   = 10; % > 1 to plot a bit faster
select_nlayers = [];%select any number of layers with []
%select_nlayers = 3;%select nl=3 layer models only

%maxdepth = [];  % Use this for full depth range
maxdepth = 400; % Or set to particular max depth

F = dir([ncdir '*.nc']);
nsoundings=length(F);
for i=1:1:nsoundings
    close all;
    ncfile = [ncdir F(i).name];
    P = read_rjmcmc_pmap(ncfile);

    %Plot 1
    figure1_handle = dark_figure(1); maximize_figure(); figure_defaults();
    [fh,ah] = plot_rjmcmc_pmap(figure1_handle,P,maxdepth,TM);
    titlestr = sprintf('#%d Line %d %.0fmE %.0fmN',i,P.line,P.x,P.y);
    axes(ah(1)); title(titlestr,'Fontsize',14);

    set(gcf,'inverthardcopy','off');
    jpgfile = sprintf('%s\\sounding_%03d_pmap.jpg',plotdir,i);
    saveas(gcf,jpgfile);

    %Plot 2
    figure2_handle = dark_figure(2); maximize_figure(); figure_defaults();

    subplot(1,2,1);
    plot_rjmcmc_data_predicted(gca,P,C,subsam,titlestr);
    xlabel('Time (s)');
    ylabel('Response (pV/Am^4)');
    xlim([10e-6 2e-2]);
    ylim([5e-3 2e2]); % data min/max

    subplot(1,2,2);
    plot_rjmcmc_models(gca,P,subsam,maxdepth,[],select_nlayers);
    
    set(gcf,'inverthardcopy','off');
    jpgfile = sprintf('%s\\sounding_%03d_predicted_models.jpg',plotdir,i);
    saveas(gcf,jpgfile);
    
end


