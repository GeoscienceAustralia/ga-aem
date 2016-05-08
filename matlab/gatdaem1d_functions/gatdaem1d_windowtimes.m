function wt = gatdaem1d_windowtimes(hS)

nw = calllib(gatdaem1d_libname(),'nwindows',  hS);

low  = zeros(nw,1);
high = zeros(nw,1);

plow  = libpointer('doublePtr',low);
phigh = libpointer('doublePtr',high);

calllib(gatdaem1d_libname(),'windowtimes',hS,plow,phigh);

low  = get(plow,'Value');
high = get(phigh,'Value');

wt.low  = low;
wt.high = high;
wt.centre = (wt.low+wt.high)/2;

