Run compile_julia_shared_lib.sh to make the shared GA-AEM library. 

Ensure submodules are pulled first:
```bash script
>> cd myrepos/ga-aem 
>> git submodule init 
>> git submodule update 
>> git submodule status 
```
Tested with Julia 1.3.1 using Cxx.jl
Please compile Cxx.jl with the same version of GCC used to compile GA-AEM
Before running for the first time ensure at BASH to export CXXJL_ROOTDIR=/path/to/gcc/version
