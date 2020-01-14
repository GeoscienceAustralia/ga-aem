clc;
clear all;

ncfile = 'output1.5p\pmaps\seq.00000001.7011.28534.000000.nc';
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
temperature = ncread(ncfile,'temperature')';
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

subplot(1,10,1:8)
imagesc(value,depth,lchist')
hold on;
plot(mean_model,depth,'-w');

subplot(1,10,10)
imagesc([0 1],depth,cphist)

box on;
%%
dark_figure(2);
bar(layer,nlhist)

%%
dark_figure(3);
plot(lchist(32,:))

%%
dark_figure(4);

subplot(3,1,1);
plot(cvs,temperature)

subplot(3,1,2);
plot(cvs,nlayers)

subplot(3,1,3);
semilogy(cvs,misfit/ndata);
ylim([0 10]);

subplot(3,1,4);
plot(cvs,ar_valuechange);

