function R = gatdaem1d_derivative(hS,derivativetype,layernumber)
      
libname = gatdaem1d_libname();
nw = calllib(libname,'nwindows',hS);
R.PX    = 0;
R.PY    = 0;
R.PZ    = 0;
R.SX    = zeros(nw,1);
R.SY    = zeros(nw,1);	
R.SZ    = zeros(nw,1);	
pR = libpointer('sTDEmResponseML',R);
calllib(libname,'derivative',hS,derivativetype,layernumber,pR);
R = get(pR,'Value');
delete(pR); 
clear pR;

