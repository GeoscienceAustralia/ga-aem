function [dxbdp,dzbdp] = gatdaem1d_rx_pitch_derivative(hS, rx_pitch, xb, zb)

if(nargin ~= 4)
    error('Sorry must be 4 arguments');
end

libname = gatdaem1d_libname();
n = length(xb);
dxbdp = zeros(n,1);
dzbdp = zeros(n,1);

ptr_xb = libpointer('doublePtr',xb);
ptr_zb = libpointer('doublePtr',zb);
ptr_dxbdp = libpointer('doublePtr',dxbdp);
ptr_dzbdp = libpointer('doublePtr',dzbdp);

calllib(libname,'derivative_rx_pitch',hS,n,rx_pitch,ptr_xb,ptr_zb,ptr_dxbdp,ptr_dzbdp); 
dxbdp = get(ptr_dxbdp,'Value');delete(ptr_dxbdp); clear ptr_dxbdp;
dzbdp = get(ptr_dzbdp,'Value');delete(ptr_dzbdp); clear ptr_dzbdp;
delete(ptr_xb); clear ptr_xb;
delete(ptr_zb); clear ptr_zb;

