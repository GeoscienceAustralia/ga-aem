using NetCDF

#Data structure to hold everything in the pmap.
struct Pmap
    depth::Array{Float64,1}
    value::Array{Float64,1}
    layer::Array{UInt,1}
    lchist::Array{UInt,2}
    cphist::Array{UInt,1}
    nlhist::Array{UInt,1}
    observations::Array{Float64,1}
    errors::Array{Float64,1}
    ndata::UInt
    chain::Array{UInt,1}
    nchains::UInt
    cvs::Array{UInt,1}
    temperature::Array{Float64,2}
    misfit::Array{Float64,2}
    nlayers::Array{UInt,2}
    logppd::Array{Float64,2}
    mean_model::Array{Float64,1}
    mode_model::Array{Float64,1}
    p10_model::Array{Float64,1}
    p50_model::Array{Float64,1}
    p90_model::Array{Float64,1}
    ar_valuechange::Array{Float64,2}
    ar_move::Array{Float64,2}
    ar_birth::Array{Float64,2}
    ar_death::Array{Float64,2}
    swap_histogram::Array{UInt,2}

    ar_noisechange::Array{Float64,2}

    noise_bins::Array{Float64,2}
    noise_counts::Array{UInt,2}
    nnoises::UInt
end

function read_rjmcmc_pmap(ncfilename::String)::Pmap

    # ncinfo(ncfilename);

    depth = ncread(ncfilename,"depth");
    value = ncread(ncfilename,"value");
    layer = ncread(ncfilename,"layer");

    lchist = permutedims(ncread(ncfilename,"log10conductivity_histogram"));
    cphist = ncread(ncfilename,"interface_depth_histogram");
    nlhist = ncread(ncfilename,"nlayers_histogram");

    observations = ncread(ncfilename,"observations");
    errors       = ncread(ncfilename,"errors");
    ndata  = length(observations);

    chain = ncread(ncfilename,"chain");
    nchains = size(chain,1);
    cvs = ncread(ncfilename,"convergence_sample");
    temperature = permutedims(ncread(ncfilename,"temperature"));
    misfit = permutedims(ncread(ncfilename,"misfit"));
    nlayers = permutedims(ncread(ncfilename,"nlayers"));
    logppd = permutedims(ncread(ncfilename,"logppd"));

    mean_model = ncread(ncfilename,"mean_model");
    mode_model = ncread(ncfilename,"mode_model");
    p10_model  = ncread(ncfilename,"p10_model");
    p50_model  = ncread(ncfilename,"p50_model");
    p90_model  = ncread(ncfilename,"p90_model");

    ar_valuechange = permutedims(ncread(ncfilename,"ar_valuechange"));
    ar_move        = permutedims(ncread(ncfilename,"ar_move"));
    ar_birth       = permutedims(ncread(ncfilename,"ar_birth"));
    ar_death       = permutedims(ncread(ncfilename,"ar_death"));
    swap_histogram = permutedims(ncread(ncfilename,"swap_histogram"));

    local ar_noisechange,noise_bins,noise_counts,nnoises;
    try
        ar_noisechange = permutedims(ncread(ncfilename,"ar_noisechange"));
        noise_bins = ncread(ncfilename,"noise_bins");
        noise_counts = ncread(ncfilename,"noise_histogram");
        nnoises = size(noise_counts,1);
    catch
        nnoises = 0;
        ar_noisechange=[];
        noise_bins = [];
        noise_counts = [];
    end

    Pmap(depth,value,layer,lchist,cphist,nlhist,observations,errors,ndata,chain,nchains,cvs,
        temperature,misfit,nlayers,logppd,mean_model,mode_model,p10_model,p50_model,p90_model,
        ar_valuechange,ar_move,ar_birth,ar_death,swap_histogram,ar_noisechange,
        noise_bins,noise_counts,nnoises)
end
