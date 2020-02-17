clc;
clear all;
close all;
addpath("..\..\..\matlab\rjmcmctdem_functions");

ncdir = 'output\pmaps\';
%ncdir = 'output_0.10\pmaps\';
%ncdir = 'output_0.20\pmaps\';
%ncdir = 'output_100k_noise_free\pmaps\';
%ncdir = 'output_100k_noise_contaminated\pmaps\';

f = dir([ncdir '*.nc']);
%ncfile  = 'output\pmaps\seq.00000001.7011.28534.000000.nc';
%ncfile  = 'output\pmaps\seq.00000002.7011.28539.000000.nc';

TM.c = [1.000000e-02    1.000000e-01    3.000000e-02    1.000000e-01    1.000000e-03];
TM.t = [2.000000e+01    1.100000e+01    5.000000e+01    3.000000e+01];

%for i=1:1:length(f)
for i=1:1:1
    ncfile = [ncdir f(i).name];
    P = read_rjmcmc_pmap(ncfile);
    [fh ah] = view_rjmcmc_pmap(P,TM);
end

%%
dark_figure()
plot(-P.observations);
hold on;
errorbar(-P.observations,P.errors);
set(gca,'yscale','log');

%%
ind = P.temperature == 1;
nmf = (P.misfit(ind)/P.ndata);
dark_figure()
%semilogy(nmf)
%mean(nmf(200:end))
plot(P.ar_valuechange(ind));
hold on;
plot(P.ar_move(ind));
ylim([0 100]);
