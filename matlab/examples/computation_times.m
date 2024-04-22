clc;
clear all;

% Timing result below computed in src/example_forward_model.cpp
% Compiler: Visual Studio 2019 compiler on Win10
% Computer: Microsoft Surface 7 Intel Core i71065G& CPU @1.30Ghz 1.50GHz CPU 
% Time to randomly popoulate c,t and forward model both Skytem HM and LM for
% examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-LM.stm	
% examples\\SkyTEM-BHMAR-2009\\stmfiles\\Skytem-HM.stm


%  nlayers and times below
%  n secs
a=[1 0.968
   2 1.243
   3 1.388
   4 1.581
   5 1.735
   6 1.846
   7 2.037
   8 2.225
   9 2.411
   10 2.616
   11 2.805
   12 3.004
   13 3.19
   14 3.441
   15 3.599
   16 3.793
   17 3.981
   18 4.175
   19 4.394
   20 4.722
   21 4.747
   22 4.937
   23 5.061
   24 5.259
   25 5.442
   26 5.814
   27 5.834
   28 6.014
   29 6.202
   30 6.395
   31 6.646
   32 6.791
   33 6.965
   34 7.195
   35 7.377
   36 7.808
   37 7.746
   38 7.947
   39 8.039
   40 8.225
   41 8.465
   42 8.666
   43 8.832
   44 9.06
   45 9.197
   46 9.402
   47 9.768
   48 9.799
   49 10.031
   50 10.185];


nlayers = a(:,1);
time    = a(:,2);

p=polyfit(nlayers,time,1)
y = p(1)*nlayers + p(2);

figure();
plot(nlayers,time,'-bo');
hold on;
plot(nlayers,y,'-r');
xlabel('#Layers');
ylabel('Time (ms) / forward model');


