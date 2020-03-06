using NetCDF

function read_rjmcmc_pmap(ncfilename)

    # ncinfo(ncfilename);

    P = Dict()

    P["depth"] = ncread(ncfilename,"depth");
    P["value"] = ncread(ncfilename,"value");
    P["layer"]  = ncread(ncfilename,"layer");

    P["lchist"] = permutedims(ncread(ncfilename,"log10conductivity_histogram"));
    P["cphist"] = ncread(ncfilename,"interface_depth_histogram");
    P["nlhist"] = ncread(ncfilename,"nlayers_histogram");

    P["observations"] = ncread(ncfilename,"observations");
    P["errors"]       = ncread(ncfilename,"errors");
    P["ndata"]  = length(P["observations"]);

    P["chain"] = ncread(ncfilename,"chain");
    P["nchains"] = size(P["chain"],1);
    P["cvs"] = ncread(ncfilename,"convergence_sample");
    P["temperature"] = permutedims(ncread(ncfilename,"temperature"));
    P["misfit"] = permutedims(ncread(ncfilename,"misfit"));
    P["nlayers"] = permutedims(ncread(ncfilename,"nlayers"));
    P["logppd"] = permutedims(ncread(ncfilename,"logppd"));

    P["mean_model"] = ncread(ncfilename,"mean_model");
    P["mode_model"] = ncread(ncfilename,"mode_model");
    P["p10_model"]  = ncread(ncfilename,"p10_model");
    P["p50_model"]  = ncread(ncfilename,"p50_model");
    P["p90_model"]  = ncread(ncfilename,"p90_model");

    P["ar_valuechange"] = permutedims(ncread(ncfilename,"ar_valuechange"));
    P["ar_move"]        = permutedims(ncread(ncfilename,"ar_move"));
    P["ar_birth"]       = permutedims(ncread(ncfilename,"ar_birth"));
    P["ar_death"]       = permutedims(ncread(ncfilename,"ar_death"));
    P["swap_histogram"] = permutedims(ncread(ncfilename,"swap_histogram"));

    P

end
