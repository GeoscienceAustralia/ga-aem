
clc;
clear all;

%% This needs to be run on a machine that has the compiler installed and in the path

libname    = 'gatdaem1d';
aliasname  = libname;
dlldir     = 'C:\Users\rossc\AppData\Local\GA-AEM-DEV\bin';
headerdir  = 'C:\Users\rossc\AppData\Local\GA-AEM-DEV\include';

dllpath       = [dlldir    '\' libname '.dll'];
headerpath    = [headerdir '\' libname '.h'];
protopath     = [libname '_proto.m'];
thunkpath     = [libname '_thunk_pcwin64.dll'];
thunklibname  = [libname '_thunk_pcwin64'];

disp(dllpath);
disp(headerpath);
disp(protopath);
disp(thunkpath);
disp(thunklibname);


if libisloaded(libname)
     disp("Unloading library");
     unloadlibrary(libname);
end
rehash;
%disp('Done...'); sound(1); return;


if libisloaded(thunklibname)
     disp("Unloading thunk library");
     unloadlibrary(thunklibname);
end

[functions warnings]  = loadlibrary(dllpath, headerpath, 'mfilename', protopath, 'alias', aliasname);

% (Optional) See what MATLAB thinks it loaded
libfunctionsview('gatdaem1d');


if libisloaded(libname)
     disp("Unloading library");
     unloadlibrary(libname);
end

if libisloaded(thunklibname)
     disp("Unloading thunk library");
     unloadlibrary(thunklibname);
end

disp('Done...');
sound(1);