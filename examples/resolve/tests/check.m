clc;
clear all;

f = load('frequencies.dat');
fm = load('fm.dat');
suffix = '_dC_000';
%suffix = '_dB';
fm1 = load(['fm1'  suffix '.dat']);
fm2 = load(['fm2'  suffix '.dat']);
dba = load(['dfma' suffix '.dat']);
dbn = load(['dfmn' suffix '.dat']);

dark_figure(1);
maximize_figure();

subplot(2,2,[1,3]);
semilogx(f,fm(:,1),'-ro');
hold on;
semilogx(f,fm(:,2),'-bo');
semilogx(f,fm1(:,1),':r');
semilogx(f,fm1(:,2),':b');
semilogx(f,fm2(:,1),':r');
semilogx(f,fm2(:,2),':b');

subplot(2,2,2);
semilogx(f,dba(:,1),'-ro');
hold on;
semilogx(f,dbn(:,1),'-c');


subplot(2,2,4);
semilogx(f,dba(:,2),'-bo');
hold on;
semilogx(f,dbn(:,2),'-m');

