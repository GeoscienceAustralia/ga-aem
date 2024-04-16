function plot_rjmcmc_data_predicted(axes_handle,P,C,subsam,titlestr)

modelcolor  = [0.5 0.5 0.5];
for j=1:1:length(C);
    dind  = C(j).dindex;
    wtime = C(j).wtime;
    for si=1:subsam:length(P.cvs);
        if(P.cvs(si) < P.nburnin) continue; end;
        for ci=1:1:P.nchains
            if(P.temperature(ci,si) ~= 1.0) continue; end;
            loglog(wtime,-P.predicted(dind,si,ci),'-','color',modelcolor,'linewidth',1);
            hold on;
        end
    end
    errorbar(wtime,-P.observations(dind),-P.errors(dind),'r');
end
xlabel('Time (s)');
ylabel('Response (pV/Am^4)');
xlim([10e-6 2e-2]);
if(~isempty(titlestr))
    title(titlestr);
end




