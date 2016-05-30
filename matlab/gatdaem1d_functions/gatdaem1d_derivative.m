function R = gatdaem1d_derivative(hS,derivativetype,layernumber)
      
libname = gatdaem1d_libname();
nw = calllib(libname,'nwindows',hS);
R.PX    = 0;
R.PY    = 0;
R.PZ    = 0;
R.SX    = zeros(nw,1);
R.SY    = zeros(nw,1);	
R.SZ    = zeros(nw,1);	

ptr_px = libpointer('doublePtr',R.PX);
ptr_py = libpointer('doublePtr',R.PY);
ptr_pz = libpointer('doublePtr',R.PZ);
ptr_sx = libpointer('doublePtr',R.SX);
ptr_sy = libpointer('doublePtr',R.SY);
ptr_sz = libpointer('doublePtr',R.SZ);

calllib(libname,'derivative',hS,derivativetype,layernumber,ptr_px,ptr_py,ptr_pz,ptr_sx,ptr_sy,ptr_sz);

R.PX = get(ptr_px,'Value');delete(ptr_px); clear ptr_px;
R.PY = get(ptr_py,'Value');delete(ptr_py); clear ptr_py;
R.PZ = get(ptr_pz,'Value');delete(ptr_pz); clear ptr_pz;
R.SX = get(ptr_sx,'Value');delete(ptr_sx); clear ptr_sx;
R.SY = get(ptr_sy,'Value');delete(ptr_sy); clear ptr_sy;
R.SZ = get(ptr_sz,'Value');delete(ptr_sz); clear ptr_sz;
