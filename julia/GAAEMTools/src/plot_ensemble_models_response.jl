using NCDatasets
using Statistics
using StatsBase
using PyPlot

function ct2cd(c::Array{<:Real,1},t::Array{<:Real,1},d::Array{<:Real,1})::Array{<:Real,1}
    interfaces = cumsum(t);
    layerinds = searchsortedfirst.(Ref(interfaces),d);
    c[layerinds]
end


function plot_swarm(ncfile::String, lmstm::String, hmstm::String)
    nc_pmap = Dataset(ncfile);
    nsamples = length(nc_pmap["ensemble_cond"]);

    #figure out number of low and high moment gates in the system
    em = init_GA_AEM(; stmfile_LM = lmstm, stmfile_HM = hmstm);
    ndata_LM = length(em.tLM);
    ndata_HM = length(em.tHM);

    #read measured geometry params from file
    dztxrx = nc_pmap.attrib["txrx_dz"];
    rRx = sqrt(nc_pmap.attrib["txrx_dx"]^2 + nc_pmap.attrib["txrx_dy"]^2);
    zTx = nc_pmap.attrib["tx_height"];

    @assert(nsamples >= 50);
    sample_idxs = sample(1:nsamples, 50, replace = false);

    cvec = map(x-> 10 .^ x, nc_pmap["ensemble_cond"][:][sample_idxs]);
    tvec = nc_pmap["ensemble_thick"][:][sample_idxs];
    d = nc_pmap["depth"][:];

    #transparency for plotting individual samples
    α = 6/length(cvec);

    fig, ax = subplots(1,2,figsize=(14,7));

    ax[1].set_xscale("log");
    ax[1].plot(ct2cd(cvec[1][:],tvec[1][:],d), d, "-k", alpha=α, label = "ensemble samples");
    for i=2:length(cvec)
        ax[1].plot(ct2cd(cvec[i][:],tvec[i][:],d), d, "-k", alpha=α);
    end

    #plot summary models over the thing
    mc = nc_pmap["p50_model"][:];
    lc = nc_pmap["p10_model"][:];
    uc = nc_pmap["p90_model"][:];
    ax[1].plot(10 .^ mc, d, "--b", label = "median model");
    ax[1].plot(10 .^ lc, d, ":b", label = "10th/90th percentile models");
    ax[1].plot(10 .^ uc, d, ":b");
    ax[1].invert_yaxis();
    ax[1].set_xlabel("Conductivity (S/m)");
    ax[1].set_ylabel("Depth (m)");
    ax[1].legend(loc="best");

    # if tx_height exists pop these in for the forward model
    tx_height_vec = nothing;
    mzTx = nothing;
    if haskey(nc_pmap, "ensemble_nuisances")
        tx_height_vec = nc_pmap["ensemble_nuisances"][:][sample_idxs];
        mzTx = median(nc_pmap["ensemble_nuisances"][:]);
    end


    #do the forward modelling and plot
    g = Geometry(ztxLM = zTx, ztxHM = zTx, rrx = -rRx,
        dzrxLM = dztxrx, dzrxHM = dztxrx);

    ax[2].set_xscale("log");
    ax[2].set_yscale("log");
    for i = 1:length(cvec)
        if !isnothing(tx_height_vec)
            g = Geometry(ztxLM = tx_height_vec[i],
                ztxHM = tx_height_vec[i], rrx = -rRx, dzrxLM = dztxrx, dzrxHM = dztxrx);
        end
        forward_data = forward_model(g, cvec[i][:], tvec[i][:],
            stmfile_LM = lmstm, stmfile_HM = hmstm);
        ax[2].plot(forward_data[1:ndata_LM,2],-forward_data[1:ndata_LM,1],"-k", alpha=α);
        ax[2].plot(forward_data[end-ndata_HM+1:end,2],-forward_data[end-ndata_HM+1:end,1],"-k", alpha=α);
    end

    if !isnothing(mzTx)
        g = Geometry(ztxLM = mzTx, ztxHM = mzTx, rrx = -rRx,
        dzrxLM = dztxrx, dzrxHM = dztxrx);
    end

    forward_data = forward_model(g, Float64(10) .^ mc, d[2:end] - d[1:end-1],
        stmfile_LM = lmstm, stmfile_HM = hmstm);
    ax[2].plot(forward_data[1:ndata_LM,2],-forward_data[1:ndata_LM,1],"--b");
    ax[2].plot(forward_data[end-ndata_HM+1:end,2],-forward_data[end-ndata_HM+1:end,1],"--b");


    ax[2].set_xlabel("Time (s)");
    ax[2].set_ylabel(raw"Signal (pV/A m$^4)$");

    obs = nc_pmap["observations"][:];
    total_noise = nc_pmap["errors"][:];

    if haskey(nc_pmap, "ensemble_noises")
        multnoises = nc_pmap["ensemble_noises"][:];
        multnoisemat = [repeat(nc_pmap["ensemble_noises"][:][1,:]',18,1);repeat(nc_pmap["ensemble_noises"][:][2,:]',23,1)];
        total_noise = median(sqrt.((multnoisemat.*obs).^2 .+ total_noise.^2),dims=2);
    end

    errorbar(forward_data[:,2], -obs[:], yerr=2*total_noise, ls="", marker=".", color="red");
end
