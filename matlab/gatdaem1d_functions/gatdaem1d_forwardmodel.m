function R = gatdaem1d_forwardmodel(hS,G,E,iptypestr)

if(nargin < 3)
    error('Sorry must be at least 3 arguments');
end

if(nargin > 4)
    error('Sorry must be no more than 4 arguments');
end


nlayers=length(E.conductivity);
if(nlayers ~= 1+length(E.thickness))
    error('Sorry the length of thicknesses must be one less than the length of conductivity');;
end

if(nargin == 4)
    if(strcmpi(iptypestr,'colecole'))
        iptype=1;
    elseif(strcmpi(iptypestr,'pelton'))
        iptype=2;
    else
        error('Sorry for IP modelling iptypestr must be either "colecole" or "pelton"');
    end

    
    if(isfield(E,'chargeability'))
        if(nlayers ~= length(E.chargeability))
            error('Sorry the length of chargeability must equal to the length of conductivity');;
        end
        
        if(isfield(E,'timeconstant'))
            if(nlayers ~= length(E.timeconstant))
                error('Sorry the length of timeconstant must be equal to the length of conductivity');;
            end
        else
            error('Sorry you must set the timeconstant');;
        end
        
        if(isfield(E,'frequencydependence'))
            if(nlayers ~= length(E.frequencydependence))
                error('Sorry the length of frequencydependence must equal to the length of conductivity');;
            end
        else
            error('Sorry you must set the frequencedependence');;
        end
    end    
end

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

if(nargin == 3)
    calllib(libname,'forwardmodel',hS,G.tx_height, G.tx_roll, G.tx_pitch, G.tx_yaw, G.txrx_dx, G.txrx_dy, G.txrx_dz, G.rx_roll, G.rx_pitch, G.rx_yaw, nlayers, E.conductivity,E.thickness,ptr_px,ptr_py,ptr_pz,ptr_sx,ptr_sy,ptr_sz);
else
    calllib(libname,'forwardmodel_ip',hS,G.tx_height, G.tx_roll, G.tx_pitch, G.tx_yaw, G.txrx_dx, G.txrx_dy, G.txrx_dz, G.rx_roll, G.rx_pitch, G.rx_yaw, nlayers, E.conductivity, E.thickness, iptype, E.chargeability, E.timeconstant, E.frequencydependence, ptr_px, ptr_py, ptr_pz, ptr_sx, ptr_sy, ptr_sz);    
end

R.PX = get(ptr_px,'Value');delete(ptr_px); clear ptr_px;
R.PY = get(ptr_py,'Value');delete(ptr_py); clear ptr_py;
R.PZ = get(ptr_pz,'Value');delete(ptr_pz); clear ptr_pz;
R.SX = get(ptr_sx,'Value');delete(ptr_sx); clear ptr_sx;
R.SY = get(ptr_sy,'Value');delete(ptr_sy); clear ptr_sy;
R.SZ = get(ptr_sz,'Value');delete(ptr_sz); clear ptr_sz;

