module GAAEMTools

#forward modelling using the C++ library
include("gaaem_forward.jl");
include("make_synthetic.jl");

#plotting and vis for rjmcmc outputs
include("rjmcmc_pmap_tools.jl");
include("plot_section.jl");

#both (postprocessing with forward modelling of ensemble members)
include("plot_ensemble_models_response.jl")

end