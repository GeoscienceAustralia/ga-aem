function plot_rjmcmc_models(axes_handle,P,subsam,maxdepth,titlestr,select_nlayers)

if(isempty(maxdepth))
   maxdepth = P.pmax;
end

modelcolor  = [0.5 0.5 0.5];
for si=1:subsam:length(P.cvs);
    if(P.cvs(si) < P.nburnin) continue; end;
    for ci=1:1:P.nchains
        if(P.temperature(ci,si) ~= 1.0) continue; end;
        nl = P.nlayers(ci,si);

        %For example only plot the "select_nlayers" layer models
        if(~isempty(select_nlayers))
            if(nl~=select_nlayers) continue; end;
        end

        clay = 10.^P.layer_value(1:nl,si,ci);
        dlay = P.layer_depth_top(1:nl,si,ci);
        [c,d] = cdtop2cd(clay,dlay,maxdepth);
        semilogx(c,d,'-','color',modelcolor,'linewidth',1);
        hold on;

        %disp(clay')
        %disp(diff(dlay'))
        %disp(dlay')
    end
end
set(gca,'xaxislocation','top');
set(gca,'ydir','reverse');
xlim(10.^[P.vmin P.vmax]);
ylim([0 maxdepth]);
xlabel('Conductivity (S/m)');
ylabel('Depth (m)');
if(~isempty(titlestr))
    title(titlestr);
end




