clc;
clear all;

addpath("..\..\..\matlab\rjmcmctdem_functions");
ncdir = 'output\pmaps\';
f = dir([ncdir '*.nc']);
%ncfile  = 'output\pmaps\seq.00000001.7011.28534.000000.nc';
%ncfile  = 'output\pmaps\seq.00000002.7011.28539.000000.nc';
close all;
%for i=1:1:length(f)
for i=1:1:1
    ncfile = [ncdir f(i).name];
    P = read_rjmcmc_pmap(ncfile);
    view_rjmcmc_pmap(P);
end
return
%%
dark_figure()
plot(-P.observations);
hold on;
errorbar(-P.observations,P.errors);
set(gca,'yscale','log');

%%
ind = P.temperature == 1;
a = (P.misfit(ind)/P.ndata);
dark_figure()
semilogy(a)
mean(a(600:end))

%%
dark_figure()
plot(P.nlhist)
hold on;
plot(double(max(P.nlhist))./(1:1:20));

