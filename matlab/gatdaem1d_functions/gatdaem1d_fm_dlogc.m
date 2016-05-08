function A = gatdaem1d_fm_dlogc(hS,G,E)

libname = gatdaem1d_libname();
nw = calllib(libname,'nwindows',hS);
nlayers = length(E.conductivity);

nchan = (1+nw);
ncomp = 3;
ncalc = (1+nlayers);
len = (nchan*ncomp*ncalc);

D = zeros(1,len);
pD = libpointer('doublePtr',D);
calllib(libname,'fm_dlogc',hS,G,nlayers,E.conductivity,E.thickness,pD);
D=get(pD,'Value');
delete(pD);
clear  pD;

D=reshape(D,[nchan,ncomp,ncalc]);
w=2:nw+1;
A.FM.PX = D(1,1,1);
A.FM.PY = D(1,2,1);
A.FM.PZ = D(1,3,1);
A.FM.SX = D(w,1,1);
A.FM.SY = D(w,2,1);
A.FM.SZ = D(w,3,1);


for k=1:1:nlayers  %must add 1
    calc=k+1;
    R.PX = D(1,1,calc);
    R.PY = D(1,2,calc);
    R.PZ = D(1,3,calc);
    R.SX = D(w,1,calc);
    R.SY = D(w,2,calc);
    R.SZ = D(w,3,calc);
    A.dlogC(k)=R;
end
