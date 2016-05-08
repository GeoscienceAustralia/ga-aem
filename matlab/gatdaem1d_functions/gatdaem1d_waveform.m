function w = gatdaem1d_waveform(hS)

ns = calllib(gatdaem1d_libname(),'nsamplesperwaveform',  hS);

t = zeros(ns,1);
c = zeros(ns,1);
v = zeros(ns,1);

pt = libpointer('doublePtr',t);
pc = libpointer('doublePtr',c);
pv = libpointer('doublePtr',v);

calllib(gatdaem1d_libname(),'waveform',hS,pt,pc,pv);

t = get(pt,'Value');
c = get(pc,'Value');
v = get(pv,'Value');

w.time    = t;
w.current = c;
w.voltage = v;