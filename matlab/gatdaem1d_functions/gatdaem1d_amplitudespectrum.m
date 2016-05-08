function [freq A] = gatdaem1d_amplitudespectrum(time,current)

dt = (time(2)-time(1));
N=length(time);


freq = (0.5/dt)*linspace(0,1,N/2+1)';
F    = fft(current)/N;
A    = abs(F(1:N/2+1));
freq = freq(1:N/2+1);

