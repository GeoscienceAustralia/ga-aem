using NetCDF, PyPlot

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

    ar_nuisancechange::Array{Float64,2}

    nuisance_bins::Array{Float64,2}
    nuisance_counts::Array{UInt,2}
    nnuisances::UInt
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
    #read noise vars
    local ar_noisechange,noise_bins,noise_counts,nnoises;
    try
        ar_noisechange = permutedims(ncread(ncfilename,"ar_noisechange"));
        noise_bins = ncread(ncfilename,"noise_bins");
        noise_counts = ncread(ncfilename,"noise_histogram");
        nnoises = size(noise_counts,2);
    catch
        nnoises = 0;
        ar_noisechange=Array{Float64,2}(undef,0,0);
        noise_bins = Array{Float64,2}(undef,0,0);
        noise_counts = Array{UInt,2}(undef,0,0);
    end
    #read nuisance vars
    local ar_nuisancechange, nuisance_bins, nuisance_counts, nnuisances;
    try
        ar_nuisancechange = permutedims(ncread(ncfilename,"ar_nuisancechange"));
        nuisance_bins = ncread(ncfilename,"nuisance_bins");
        nuisance_counts = ncread(ncfilename,"nuisance_histogram");
        nnuisances = size(nuisance_counts,2);
    catch
        nnuisances = 0;
        ar_nuisancechange=Array{Float64,2}(undef,0,0);
        nuisance_bins = Array{Float64,2}(undef,0,0);
        nuisance_counts = Array{UInt,2}(undef,0,0);
    end

    Pmap(depth,value,layer,lchist,cphist,nlhist,observations,errors,ndata,chain,nchains,cvs,
        temperature,misfit,nlayers,logppd,mean_model,mode_model,p10_model,p50_model,p90_model,
        ar_valuechange,ar_move,ar_birth,ar_death,swap_histogram,ar_noisechange,
        noise_bins,noise_counts,nnoises, ar_nuisancechange, nuisance_bins,
        nuisance_counts, nnuisances)
end


function ct2cd(cond,thick,depth)
    interfaces = cumsum(thick)
    layerinds = searchsortedfirst.(Ref(interfaces),depth)
    cond[layerinds]
end

function view_rjmcmc_pmap(P::Pmap, TM=Dict())

    p = [
        0.1 0.50 0.775 0.980;
        0.1 0.50 0.550 0.755;
        0.1 0.50 0.325 0.530;
        0.1 0.50 0.100 0.305;

        0.57 0.99 0.82 0.98;
        0.57 0.94 0.10 0.76;
        0.95 0.99 0.10 0.76;
        ]

    if (P.nnoises > 0 && P.nnuisances > 0)
        p[6,:] = [0.57 0.94 0.42 0.76];
        p[7,:] = [0.95 0.99 0.42 0.76];
        p = [p; 0.6 0.99 0.10 0.19; 0.6 0.99 0.26 0.35];
    elseif (P.nnoises > 0 || P.nnuisances > 0)
        p[6,:] = [0.57 0.94 0.32 0.76];
        p[7,:] = [0.95 0.99 0.32 0.76];
        p = [p; 0.6 0.99 0.10 0.25];
    end


    fig = figure(figsize=(10,8))

    ax = Array{Any}(undef,9);

    for i=1:size(p,1)
        ap = [p[i,1];p[i,3];p[i,2]-p[i,1];p[i,4]-p[i,3]]
        ax[i] = PyPlot.axes(position=ap)
    end

    cmap = ColorMap("jet")
    #set zero values of the histogram to black
    #matlab lets you just do this with arrays
    #but PyPlot via julia has to make things hard
    cmap._i_under = 0
    cmap.set_under((0,0,0,1))

    sca(ax[6])
    pcolormesh(10 .^ P.value,P.depth,P.lchist,cmap=cmap)
    plot(10 .^ P.mean_model,P.depth,"-w");
    plot(10 .^ P.p10_model,P.depth,":m")
    plot(10 .^ P.p50_model,P.depth,"-m")
    plot(10 .^ P.p90_model,P.depth,":m")
    if length(TM) > 0
        tmcmap = ct2cd(TM["conductivity"],TM["thickness"],P.depth)
        plot(tmcmap,P.depth,"-g")
    end
    ax[6].set_xscale("log")
    ylabel("Depth (m)")
    xlabel("Conductivity (S/m)")
    gca().invert_yaxis()
    #gca().xaxis.tick_top()
    #gca().xaxis.set_label_position("top")

    sca(ax[7])
    pcolormesh([0 1],P.depth, reshape(P.cphist,(size(P.cphist)...,1)),cmap=cmap)
    gca().invert_yaxis()

    sca(ax[5])
    bar(P.layer,P.nlhist,color="r")
    xlabel("Number of layers")
    ylim([0;maximum(P.nlhist)])

    sca(ax[1]);
    semilogy(P.cvs,permutedims(P.misfit/P.ndata));
    ylabel("Normalised misfit");
    ylim([0.09;100]);
    xticks([])

    sca(ax[2]);
    semilogy(P.cvs,permutedims(P.temperature));
    ylim([0.9; maximum(P.temperature)*1.1]);
    ylabel("Temperature");
    xticks([])

    sca(ax[3]);
    plot(P.cvs,permutedims(P.nlayers))
    ylabel("#Layers");
    xticks([])

    sca(ax[4]);
    plot(P.cvs,P.ar_valuechange[1,:],"-m")
    plot(P.cvs,P.ar_move[1,:],"-b")
    plot(P.cvs,P.ar_birth[1,:],"-g")
    plot(P.cvs,P.ar_death[1,:],"-r")
    if P.nnoises > 0
        plot(P.cvs,P.ar_noisechange[1,:],"-k")
    end
    if P.nnuisances > 0
        plot(P.cvs,P.ar_nuisancechange[1,:],"-y")
    end


    xlabel("Sample#");
    ylabel("Accept rate");
    ylim([0;100]);
    legend(["change";"move";"birth";"death"],ncol=4)

    if (P.nnoises > 0)
        sca(ax[8]);
        for i=1:P.nnoises
            plot(P.noise_bins[:,i],P.noise_counts[:,i]);
        end
        xlabel("multiplicative noise magnitude");
        ylabel("counts");
    end
    if (P.nnuisances > 0)
        sca(ax[8+(P.nnoises > 0)]);
        for i=1:P.nnuisances
            plot(P.nuisance_bins[:,i],P.nuisance_counts[:,i]);
        end
        xlabel("nuisance magnitude");
        ylabel("counts");
    end

    gcf()
end
