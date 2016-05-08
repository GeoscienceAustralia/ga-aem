function gatdaem1d_loadlibrary()        
    libname = gatdaem1d_libname();
    headerfile = gatdaem1d_headerfile();
    if(~libisloaded(libname))
        [notfound,warnings] = loadlibrary(libname, headerfile);
        %libfunctions(libname,'-full');
    end
end



