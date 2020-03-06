function [figure_handle ax] = view_rjmcmc_pmap(P,TM)

%        [ x1   x2   y1   y2 ];
p(1,:) = [0.05 0.50 0.775 0.980];
p(2,:) = [0.05 0.50 0.550 0.755];
p(3,:) = [0.05 0.50 0.325 0.530];
p(4,:) = [0.05 0.50 0.100 0.305];

p(5,:) = [0.55 0.99 0.82 0.98];
p(6,:) = [0.55 0.94 0.05 0.70];
p(7,:) = [0.95 0.99 0.05 0.70];

%%
figure_handle = dark_figure();
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
imagesc(P.value,P.depth,P.lchist)
hold on;
plot(P.mean_model,P.depth,'-w');
%plot(P.mode_model,P.depth,'-y');
plot(P.p10_model,P.depth,':m');
plot(P.p50_model,P.depth,'-m');
plot(P.p90_model,P.depth,':m');
if(nargin>1)
    TM.c = [1.000000e-02    1.000000e-01    3.000000e-02    1.000000e-01    1.000000e-03];
    TM.t = [2.000000e+01    1.100000e+01    5.000000e+01    3.000000e+01];
    [TM.tc,TM.td] = ct2cd(TM.c,TM.t,400);    
    plot(log10(TM.tc),TM.td,'-g');
end
ylabel('Depth (m)');
xlabel('Conductivity (S/m)');
set(gca,'xaxislocation','top');
set(gca,'xtick',[-3:1:1]);
set(gca,'xticklabel',10.^[-3:1:1]);

axes(ax(7));
imagesc([0 1],P.depth,P.cphist)
box on;
set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);

axes(ax(5));;
bar(P.layer,P.nlhist,'r');
xlabel('Number of layers');
ylim([0 max(P.nlhist)]);
%set(gca,'xticklabel',[]);
set(gca,'yticklabel',[]);


%%%
axes(ax(1));
semilogy(P.cvs,P.misfit/P.ndata);
%ylim([0 10]);
ylabel('Normalised misfit');
set(gca,'xticklabel',[]);
set(gca,'layer','top');
ylim([0.09 100]);

axes(ax(2));
semilogy(P.cvs,P.temperature);
ylim([0.9 max(P.temperature(:))*1.1]);
ylabel('Temperature');
set(gca,'xticklabel',[]);
set(gca,'layer','top');


axes(ax(3));
plot(P.cvs,P.nlayers)
ylabel('#Layers');
set(gca,'xticklabel',[]);
set(gca,'layer','top');

axes(ax(4));
h(1)=plot(P.cvs,P.ar_valuechange(1,:),'-m');lab{1}='value';
hold on;
h(2)=plot(P.cvs,P.ar_move(1,:),'-b');lab{2}='move';
h(3)=plot(P.cvs,P.ar_birth(1,:),'-g');lab{3}='birth';
h(4)=plot(P.cvs,P.ar_death(1,:),'-r');lab{4}='death';

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



