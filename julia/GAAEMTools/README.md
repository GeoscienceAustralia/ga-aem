## GAAEMTools

Package to enable integration of forward modelling functionality in GA-AEM with Julia, as well as plotting tools for visualising inversion outputs.

### Installation

This package requires Julia 1.3.x or earlier due to its dependence on [Cxx.jl](https://github.com/JuliaInterop/Cxx.jl/) which does not yet support Julia 1.4+.

Only Linux is supported currently. Installation was tested on CentOS 6.10 using Julia 1.3.1 and gcc 4.4.7. 

Ensure submodules are pulled first. Assuming you have cloned the ga-aem git repository to your home directory:
```bash script
>> cd ~/ga-aem 
>> git submodule init 
>> git submodule update 
>> git submodule status 
```

Then enter the `julia` subdirectory of the git repository and start julia. To install the package, enter the `pkg` REPL using `]` and type:
```julia
(v1.3) pkg> dev GAAEMTools/
```
The installation will compile a shared library that exposes the forward modelling functionality from GA-AEM to Julia. Because this compilation step depends on C++ source files from the parent repository, currently this Julia package can only be installed in development mode from its subdirectory within the GA-AEM repository.

#### For NCI VDI users
Before running for the first time ensure to set the environment variable CXXJL_ROOTDIR=/path/to/gcc/version

### Usage

After running the installation step above, you can load GAAEMTools in a Julia script or at the REPL by typing:
```julia
using GAAEMTools
```

#### Forward modelling
Forward modelling in GAAEMTools requires you to create a "geometry" object, which stores key parameters of the AEM system you are trying to model.