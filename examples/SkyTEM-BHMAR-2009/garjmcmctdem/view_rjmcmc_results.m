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
C(1).name   = 'LM-Z';%name of the component
C(1).dindex = 1:18;%index in the data array for this component
C(1).wtime =  [1.7195e-05  2.1695e-05  2.7695e-05  3.5195e-05  4.4195e-05  5.5695e-05  7.0195e-05  8.8695e-05  0.00011219  0.00014119   0.0001782   0.0002247   0.0002827   0.0003562   0.0004487  0.00056519   0.0007117   0.0008962];%window centre times

%Add other components as required
C(2).name   = 'HM-Z';%name of the component
C(2).dindex = 19:39;%index in the data array for this component
C(2).wtime =  [8.5695e-05   0.0001092   0.0001382   0.0001752  0.00022169  0.00027969  0.00035319   0.0004457  0.00056219   0.0007087  0.00089319   0.0011257   0.0014182   0.0017862   0.0022497   0.0028332   0.0035677   0.0044927   0.0056572   0.0071227   0.0088392];%window centre times

subsam   = 10; % > 1 to plot a bit faster
select_nlayers = [];%select any number of layers with []
%select_nlayers = 3;%select nl=3 layer models only

maxdepth = [];  % Use this for full depth range
%maxdepth = 200; % Or set to particular max depth

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
    ylabel('Response (V/Am^4)');
    xlim([10e-6 2e-2]);
    ylim([1e-13 2e-8]); % data min/max

    subplot(1,2,2);
    plot_rjmcmc_models(gca,P,subsam,maxdepth,[],select_nlayers);
    
    set(gcf,'inverthardcopy','off');
    jpgfile = sprintf('%s\\sounding_%03d_predicted_models.jpg',plotdir,i);
    saveas(gcf,jpgfile);
    
end


