using NetCDF

function read_rjmcmc_pmap(ncfilename)

#=                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
#        [ x1   x2   y1   y2 ];
p(1,:) = [0.05 0.50 0.775 0.980];
p(2,:) = [0.05 0.50 0.550 0.755];
p(3,:) = [0.05 0.50 0.325 0.530];
p(4,:) = [0.05 0.50 0.100 0.305];

p(5,:) = [0.55 0.99 0.82 0.98];
p(6,:) = [0.55 0.94 0.05 0.70];
p(7,:) = [0.95 0.99 0.05 0.70];
=#


    # ncinfo(ncfilename);

    P = Dict()

    P["depth"] = ncread(ncfilename,"depth");
    P["value"] = ncread(ncfilename,"value");
    P["layer"]  = ncread(ncfilename,"layer");

    P["lchist"] = ncread(ncfilename,"log10conductivity_histogram");
    P["cphist"] = ncread(ncfilename,"interface_depth_histogram");
    P["nlhist"] = ncread(ncfilename,"nlayers_histogram");

    P["observations"] = ncread(ncfilename,"observations");
    P["errors"]       = ncread(ncfilename,"errors");
    P["ndata"]  = length(P["observations"]);

    P["chain"] = ncread(ncfilename,"chain");
    P["nchains"] = size(P["chain"],1);
    P["cvs"] = ncread(ncfilename,"convergence_sample");
    P["temperature"] = ncread(ncfilename,"temperature");
    P["misfit"] = ncread(ncfilename,"misfit");
    P["nlayers"] = ncread(ncfilename,"nlayers");
    P["logppd"] = ncread(ncfilename,"logppd");

    P["mean_model"] = ncread(ncfilename,"mean_model");
    P["mode_model"] = ncread(ncfilename,"mode_model");
    P["p10_model"]  = ncread(ncfilename,"p10_model");
    P["p50_model"]  = ncread(ncfilename,"p50_model");
    P["p90_model"]  = ncread(ncfilename,"p90_model");

    P["ar_valuechange"] = ncread(ncfilename,"ar_valuechange");
    P["ar_move"]        = ncread(ncfilename,"ar_move");
    P["ar_birth"]       = ncread(ncfilename,"ar_birth");
    P["ar_death"]       = ncread(ncfilename,"ar_death");
    P["swap_histogram"] = ncread(ncfilename,"swap_histogram");

    P

end
#=
dark_figure();
maximize_figure();
for i=1:1:size(p,1);
    ap = [p(i,1) p(i,3) p(i,2)-p(i,1) p(i,4)-p(i,3)];
    ax(i) = axes('position',ap);
    box on;    
end
cmap = jet(256); cmap(1,:)=[0 0 0];
colormap(cmap);

%%
axes(ax(6));
imagesc(value,depth,lchist)
hold on;
plot(mean_model,depth,'-w');
%plot(mode_model,depth,'-y');
plot(p10_model,depth,':m');
plot(p50_model,depth,'-m');
plot(p90_model,depth,':m');

ylabel('Depth (m)');
xlabel('Conductivity (S/m)');
set(gca,'xaxislocation','top');
set(gca,'xtick',[-3:1:1]);
set(gca,'xticklabel',10.^[-3:1:1]);

axes(ax(7));
imagesc([0 1],depth,cphist)
box on;
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);

axes(ax(5));;
bar(layer,nlhist,'r');
xlabel('Number of layers');
ylim([0 max(nlhist)]);
%set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);


%%%
axes(ax(1));
semilogy(cvs,misfit/ndata);
%ylim([0 10]);
ylabel('Misfit');
set(gca,'xticklabel',[]);
set(gca,'layer','top');
ylim([0.09 100]);

axes(ax(2));
semilogy(cvs,temperature);
ylim([0.9 max(temperature(:))*1.1]);
ylabel('Temperature');
set(gca,'xticklabel',[]);
set(gca,'layer','top');


axes(ax(3));
plot(cvs,nlayers)
ylabel('#Layers');
set(gca,'xticklabel',[]);
set(gca,'layer','top');

axes(ax(4));
h(1)=plot(cvs,ar_valuechange(1,:),'-m');lab{1}='value';
hold on;
h(2)=plot(cvs,ar_move(1,:),'-b');lab{2}='move';
h(3)=plot(cvs,ar_birth(1,:),'-g');lab{3}='birth';
h(4)=plot(cvs,ar_death(1,:),'-r');lab{4}='death';

%plot(cvs,ar_valuechange,'-m');
%plot(cvs,ar_move,'-b');
%plot(cvs,ar_birth,'-g');
%plot(cvs,ar_death,'-r');

xlabel('Sample#');
set(gca,'layer','top');
ylabel('Accept rate');
ylim([0 100]);
set(get(gca,'xaxis'),'ExponentMode','manual');
set(get(gca,'xaxis'),'Exponent',0);
lh=legend(h,lab);
set(lh,'orientation','horizontal');
=#



