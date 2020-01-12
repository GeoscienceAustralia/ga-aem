clc;
clear all;

ncfile = 'output\pmaps\seq.00000001.7011.28534.000000.nc';
%ncfile = 'output\pmaps\seq.00000100.7011.29029.000000.nc';
%ncfile = 'output\pmaps\seq.00000200.7011.29529.000000.nc';

ncinfo(ncfile);
depth = ncread(ncfile,'depth');
value = ncread(ncfile,'value');
layer  = ncread(ncfile,'layer');

lchist = ncread(ncfile,'log10conductivity_histogram');
cphist = ncread(ncfile,'interface_depth_histogram');
nlhist = ncread(ncfile,'nlayers_histogram');
obs   = ncread(ncfile,'observations');
ndata  = length(obs);

chain = ncread(ncfile,'chain');
cvs = ncread(ncfile,'convergence_sample');
misfit = ncread(ncfile,'misfit')';
nlayers = ncread(ncfile,'nlayers')';
logppd = ncread(ncfile,'logppd')';
ar_valuechange = ncread(ncfile,'ar_valuechange')';

%%
mean_model = sum(double(lchist) .* value) ./ sum(lchist);

%%
dark_figure(1);
cmap = jet(256); cmap(1,:)=[0 0 0];
colormap(cmap);

imagesc(value,depth,lchist')
hold on;
plot(mean_model,depth,'-w');
%%
dark_figure(2);
bar(layer,nlhist)

%%
dark_figure(3);
plot(lchist(32,:))

%%
dark_figure(4);
subplot(3,1,1);
plot(cvs,nlayers)

subplot(3,1,2);
semilogy(cvs,misfit/ndata);
ylim([0 10]);

subplot(3,1,3);
plot(cvs,ar_valuechange);

