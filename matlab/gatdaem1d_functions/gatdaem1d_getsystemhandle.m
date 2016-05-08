function hS = gatdaem1d_getsystemhandle(stmfile)

if(exist(stmfile,'file'))
    hS = calllib(gatdaem1d_libname(),'createhandle', stmfile);
else    
    error('Sorry but the specified stmfile %s does not exist',stmfile);;
end

