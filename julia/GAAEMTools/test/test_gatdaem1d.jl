using Libdl
using Cxx
#CXXJL_ROOTDIR should be set at least when first time installing Cxx
#export CXXJL_ROOTDIR=/apps/gcc/5.2.0

# Importing shared library and header file
const path_to_lib = pwd()
addHeaderDir(path_to_lib, kind=C_System)

# Julia compiler also needs to find headers in these dirs
addHeaderDir("../src", kind=C_System)
addHeaderDir("../submodules/cpp-utils/src", kind=C_System)
addHeaderDir("/apps/fftw3/3.3.7-gcc/include", kind=C_System)

cxxinclude("gatdaem1d_julia.h")
Libdl.dlopen(path_to_lib * "/gatdaem1d_julia.so", Libdl.RTLD_GLOBAL)

# Create the cTDEmSystem class object
stmfile = "../examples/bhmar-skytem/stmfiles/Skytem-LM.stm";
TLM = @cxxnew cTDEmSystem(pointer(stmfile));
stmfile = "../examples/bhmar-skytem/stmfiles/Skytem-HM.stm";
THM = @cxxnew cTDEmSystem(pointer(stmfile));
#stmfile = "../examples/thomson-vtem/stmfiles/VTEM-plus-7.3ms-pulse-southernthomson.stm";

# Set the input arrays
nwindowsLM = @cxx TLM->NumberOfWindows;
nwindowsHM = @cxx THM->NumberOfWindows;
geometryLM   = [35.0,0.0,0.0,0.0,-17.0,0.0,+2.0,0.0,0.0,0.0];
geometryHM   = [35.0,0.0,0.0,0.0,-17.0,0.0,+0.2,0.0,0.0,0.0];
conductivity = [0.01, 1.0, 0.001];
thickness    = [20.0, 40.0];
nlayers      = length(conductivity);

# Allocate reusable storage for outputs
PX = 0.0;
PY = 0.0;
PZ = 0.0;
SXLM = Array{Float64}(undef,nwindowsLM);
SYLM = Array{Float64}(undef,nwindowsLM);
SZLM = Array{Float64}(undef,nwindowsLM);

SXHM = Array{Float64}(undef,nwindowsHM);
SYHM = Array{Float64}(undef,nwindowsHM);
SZHM = Array{Float64}(undef,nwindowsHM);

nloops = 1000;

# Forward model with separate calls
#@time for i=1:nloops
#	@cxx T->setgeometry(pointer(geometry));
#	@cxx T->setconductivitythickness(nlayers,pointer(conductivity),pointer(thickness));
#	@cxx T->setupcomputations();
#	@cxx T->setprimaryfields();
#	@cxx T->setsecondaryfields();
#	@cxx T->getfields(Ref(PX),Ref(PY),Ref(PZ),pointer(SX),pointer(SY),pointer(SZ));
#end

#@cxx TLM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),pointer(geometry),Ref(PX),Ref(PY),Ref(PZ),pointer(SXLM),pointer(SYLM),pointer(SZLM));
#@cxx THM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),pointer(geometry),Ref(PX),Ref(PY),Ref(PZ),pointer(SXHM),pointer(SYHM),pointer(SZHM));
# Forward model with single call (maybe less overhead?)
@time for i=1:nloops
    @cxx TLM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),pointer(geometryLM),Ref(PX),Ref(PY),Ref(PZ),pointer(SXLM),pointer(SYLM),pointer(SZLM));
    @cxx THM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),pointer(geometryHM),Ref(PX),Ref(PY),Ref(PZ),pointer(SXHM),pointer(SYHM),pointer(SZHM));
end


#println("nwindows: $nwindows \n")
#println("Primary X: PX: $PX \n")
#println("Primary Y: PY: $PY \n")
#println("Primary Z: PZ: $PZ \n")
#println("Secondary X: SX: $SX \n")
#println("Secondary Y: SY: $SY \n")
#println("Secondary Z: SZ: $SZ \n")


