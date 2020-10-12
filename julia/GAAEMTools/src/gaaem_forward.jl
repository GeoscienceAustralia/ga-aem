module UseGA_AEM
using Libdl
using Cxx

mutable struct EMforward
    TLM  :: Cxx.CxxCore.CppPtr
    THM  :: Cxx.CxxCore.CppPtr
    PX   :: Float64
    PY   :: Float64
    PZ   :: Float64
    SXLM :: Array{Float64, 1}
    SYLM :: Array{Float64, 1}
    SZLM :: Array{Float64, 1}
    SXHM :: Array{Float64, 1}
    SYHM :: Array{Float64, 1}
    SZHM :: Array{Float64, 1}
    tLM  :: Array{Float64, 1} # only because I don't know how to get time from TLM C++ class
    tHM  :: Array{Float64, 1}
end


function init_GA_AEM(; stmfile_LM = (@__DIR__)*"/../examples/bhmar-skytem/stmfiles/Skytem-LM.stm",
                       stmfile_HM = (@__DIR__)*"/../examples/bhmar-skytem/stmfiles/Skytem-HM.stm")
    #CXXJL_ROOTDIR should be set at least when first time installing Cxx
    #export CXXJL_ROOTDIR=/apps/gcc/5.2.0

    # Importing shared library and header file
    path_to_lib = (@__DIR__)*"/"
    addHeaderDir(path_to_lib, kind=C_System)

    ga_aem_path = abspath(path_to_lib*"../")

    # Julia compiler also needs to find headers in these dirs
    addHeaderDir(ga_aem_path*"/src", kind=C_System)
    addHeaderDir(ga_aem_path*"/submodules/cpp-utils/src/", kind=C_System)
    addHeaderDir("/usr/local/include", kind=C_System)


    # needs fftw3 as well
    if "FFTW_DIR" in keys(ENV)
        addHeaderDir(ENV["FFTW_DIR"]*"/include/", kind=C_System)
        Libdl.dlopen(ENV["FFTW_DIR"]*"/lib/libfftw3.so",Libdl.RTLD_GLOBAL)
    else
        Libdl.dlopen("libfftw3.so", Libdl.RTLD_GLOBAL)
    end

    cxxinclude("gatdaem1d_julia.h")
    Libdl.dlopen(path_to_lib * "/gatdaem1d_julia.so", Libdl.RTLD_GLOBAL)


    # Create the cTDEmSystem class object
    TLM = @cxxnew cTDEmSystem(pointer(stmfile_LM));
    THM = @cxxnew cTDEmSystem(pointer(stmfile_HM));

    # Allocate reusable storage for outputs
    PX = 0.0
    PY = 0.0
    PZ = 0.0
    nwindowsLM = @cxx TLM->NumberOfWindows;
    nwindowsHM = @cxx THM->NumberOfWindows;
    SXLM = Array{Float64, 1}(undef,nwindowsLM)
    SYLM = Array{Float64, 1}(undef,nwindowsLM)
    SZLM = Array{Float64, 1}(undef,nwindowsLM)

    SXHM = Array{Float64, 1}(undef,nwindowsHM)
    SYHM = Array{Float64, 1}(undef,nwindowsHM)
    SZHM = Array{Float64, 1}(undef,nwindowsHM)

    # get the centres of the window times
    cxx"""
    #include <vector>
    #include <cmath>
    #ifndef WINDOW_CENTRES_H
    #define WINDOW_CENTRES_H
    using namespace std;

    vector<double> window_centres(const cTDEmSystem* sys) {
        vector<double> centre;
        for (size_t i = 0; i < sys->NumberOfWindows; i++) {
            centre.push_back(sqrt(sys->WinSpec[i].TimeHigh * sys->WinSpec[i].TimeLow));
        }
        return centre;
    }
    #endif
    """

    tLM = collect(@cxx window_centres(TLM));
    tHM = collect(@cxx window_centres(THM));

    # call it once bug
    # Set the input arrays
    geometryLM   = [35.0,0.0,0.0,0.0,-17.0,0.0,+2.0,0.0,0.0,0.0]
    geometryHM   = [35.0,0.0,0.0,0.0,-17.0,0.0,+0.2,0.0,0.0,0.0]
    conductivity = [0.01, 1.0, 0.001]
    thickness    = [20.0, 40.0]
    nlayers      = length(conductivity)

    @cxx TLM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),
                           pointer(geometryLM),Ref(PX),Ref(PY),Ref(PZ),
                           pointer(SXLM),pointer(SYLM),pointer(SZLM))
    @cxx THM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),
                           pointer(geometryHM),Ref(PX),Ref(PY),Ref(PZ),
                           pointer(SXHM),pointer(SYHM),pointer(SZHM))

    EMforward(TLM, THM, PX, PY, PZ, SXLM, SYLM, SZLM, SXHM, SYHM, SZHM, tLM, tHM)
end

function(em::EMforward)()

    # Set the input arrays
    geometryLM   = [35.0,0.0,0.0,0.0,-17.0,0.0,+2.0,0.0,0.0,0.0]
    geometryHM   = [35.0,0.0,0.0,0.0,-17.0,0.0,+0.2,0.0,0.0,0.0]
    conductivity = [0.01, 1.0, 0.001]
    thickness    = [20.0, 40.0]
    nlayers      = length(conductivity)

    @cxx em.TLM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),
                           pointer(geometryLM),Ref(em.PX),Ref(em.PY),Ref(em.PZ),
                           pointer(em.SXLM),pointer(em.SYLM),pointer(em.SZLM))
    @cxx em.THM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),
                           pointer(geometryHM),Ref(em.PX),Ref(em.PY),Ref(em.PZ),
                           pointer(em.SXHM),pointer(em.SYHM),pointer(em.SZHM))
end

mutable struct Geometry
    geometryLM  :: Array{Float64, 1}
    geometryHM  :: Array{Float64, 1}
end

function Geometry(;ztxLM  = 35.0,
                   ztxHM  = 35.0,
                   rrx  =  -17.0,
                   dzrxLM =  2.0,
                   dzrxHM =  0.2
                 )
    @assert ztxLM  > 0.0
    @assert ztxHM  > 0.0
    @assert rrx  < 0.0
    @assert dzrxLM > 0.0
    @assert dzrxHM > 0.0

    geometryLM   = [ztxLM,0.0,0.0,0.0,rrx,0.0,dzrxLM,0.0,0.0,0.0]
    geometryHM   = [ztxHM,0.0,0.0,0.0,rrx,0.0,dzrxHM,0.0,0.0,0.0]

    Geometry(geometryLM, geometryHM)
end


function (em::EMforward)(g::Geometry, conductivity::Array{Float64, 1}, thickness::Array{Float64, 1})
    nlayers      = length(conductivity)
    @cxx em.TLM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),
                           pointer(g.geometryLM),Ref(em.PX),Ref(em.PY),Ref(em.PZ),
                           pointer(em.SXLM),pointer(em.SYLM),pointer(em.SZLM))
    @cxx em.THM->forwardmodel(nlayers,pointer(conductivity),pointer(thickness),
                           pointer(g.geometryHM),Ref(em.PX),Ref(em.PY),Ref(em.PZ),
                           pointer(em.SXHM),pointer(em.SYHM),pointer(em.SZHM))
end

mutable struct EMoperator
    em :: EMforward
    g  :: Geometry
end

function (op::EMoperator)(ztxLM::Float64, ztxHM::Float64, conductivity::Array{Float64, 1}, thickness::Array{Float64,1 })
    op.g.geometryLM[1] = ztxLM
    op.g.geometryHM[1] = ztxHM
    op.em(op.g, conductivity, thickness)
end

mutable struct Sounding
    dataLM    :: Array{Float64, 1}
    dataHM    :: Array{Float64, 1}
    sdLM      :: Array{Float64, 1}
    sdHM      :: Array{Float64, 1}
    x         :: Array{Float64, 1}
    ztx       :: Float64
    thickness :: Array{Float64, 1}
end

end
