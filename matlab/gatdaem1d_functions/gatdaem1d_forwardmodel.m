function R = gatdaem1d_forwardmodel(hS,G,E)

nlayers=length(E.conductivity);
if(nlayers ~= 1+length(E.thickness))
    error('Sorry the number of thicknesses must by one less than the number of conductivities');;
end
      
libname = gatdaem1d_libname();
nw = calllib(libname,'nwindows',hS);
R.PX    = 0;
R.PY    = 0;
R.PZ    = 0;
R.SX    = zeros(nw,1);
R.SY    = zeros(nw,1);	
R.SZ    = zeros(nw,1);	
pR = libpointer('sTDEmResponseML',R);
calllib(libname,'forwardmodel',hS,G,nlayers,E.conductivity,E.thickness,pR);
R = get(pR,'Value');
delete(pR); 
clear pR;




