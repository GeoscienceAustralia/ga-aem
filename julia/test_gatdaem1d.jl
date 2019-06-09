using Libdl
using Cxx
#CXXJL_ROOTDIR should be set at least when first time installing Cxx
#export CXXJL_ROOTDIR=/apps/gcc/5.2.0

# Importing shared library and header file
const path_to_lib = pwd()
addHeaderDir(path_to_lib, kind=C_System)

# Julia compiler also needs to find headers in these dirs
addHeaderDir("/short/cr78/rcb547/code/repos/ga-aem-develop/src", kind=C_System)
addHeaderDir("/short/cr78/rcb547/code/repos/ga-aem-develop/submodules/cpp-utils/src", kind=C_System)
addHeaderDir("/apps/fftw3/3.3.7-gcc/include", kind=C_System)

cxxinclude("gatdaem1d_julia.h")
Libdl.dlopen(path_to_lib * "/gatdaem1d_julia.so", Libdl.RTLD_GLOBAL)

# Create the cTDEmSystem class object
# stmfile = "../examples/bhmar-skytem/stmfiles/Skytem-LM.stm";
# stmfile = "../examples/bhmar-skytem/stmfiles/Skytem-HM.stm";
stmfile = "../examples/thomson-vtem/stmfiles/VTEM-plus-7.3ms-pulse-southernthomson.stm";
T = @cxxnew cTDEmSystem(pointer(stmfile));

# Set the input arrays
nwindows = @cxx T->NumberOfWindows;
geometry     = [35.0,0.0,0.0,0.0,-12.2,0.0,+2.5,0.0,0.0,0.0];
conductivity = [0.01, 1.0, 0.001];
thickness    = [20.0, 40.0];
nlayers      = length(conductivity);

# Allocate reusable storage for outputs
PX = 0.0;
PY = 0.0;
PZ = 0.0;
SX = Array{Float64}(undef,nwindows);
SY = Array{Float64}(undef,nwindows);
SZ = Array{Float64}(undef,nwindows);

nloops = 1000;

# Forward model with separate calls
@time for i=1:nloops
	@cxx T->setgeometry(pointer(geometry));
	@cxx T->setconductivitythickness(nlayers,pointer(conductivity),pointer(thickness));
	@cxx T->setupcomputations();
	@cxx T->setprimaryfields();
	@cxx T->setsecondaryfields();
	@cxx T->getfields(Ref(PX),Ref(PY),Ref(PZ),pointer(SX),pointer(SY),pointer(SZ));
end

# Forward model with single call (maybe less overhead?)
@time for i=1:nloops
	@cxx T->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),pointer(geometry),Ref(PX),Ref(PY),Ref(PZ),pointer(SX),pointer(SY),pointer(SZ));
end


#println("nwindows: $nwindows \n")
#println("Primary X: PX: $PX \n")
#println("Primary Y: PY: $PY \n")
#println("Primary Z: PZ: $PZ \n")
#println("Secondary X: SX: $SX \n")
#println("Secondary Y: SY: $SY \n")
println("Secondary Z: SZ: $SZ \n")

