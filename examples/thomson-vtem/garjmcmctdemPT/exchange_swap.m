clc;
clear all;

Pi = 1;
Pj = 10;
Ti = 1.0;

dark_figure(1);
for Tj=1:1:10
    
    logar = (1.0 / Ti - 1.0 / Tj) * (Pi - Pj);
    
    N=10000;
    n=0;
    for j=1:1:N
        logu = log(rand());
        if(logu < logar)
            n=n+1;
        end
    end
    aprob=n/N;
    plot(Tj,aprob,'ro');    
    hold on;
    ylim([0 1]);
end

