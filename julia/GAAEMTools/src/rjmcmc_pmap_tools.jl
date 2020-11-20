using NCDatasets, PyPlot

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

	#metadata
	fiducial::Float64
	x::Float64
	y::Float64
	elevation::Float64
end

function read_rjmcmc_pmap(ncfilename::String)::Pmap

    nc_pmap = Dataset(ncfilename);

    depth = nc_pmap["depth"][:];
    value = nc_pmap["value"][:];
    layer = nc_pmap["layer"][:];

    lchist = permutedims(nc_pmap["log10conductivity_histogram"][:]);
    cphist = nc_pmap["interface_depth_histogram"][:];
    nlhist = nc_pmap["nlayers_histogram"][:];

    observations = nc_pmap["observations"][:];
    errors       = nc_pmap["errors"][:];
    ndata  = length(observations);

    chain = nc_pmap["chain"][:];
    nchains = size(chain,1);
    cvs = nc_pmap["convergence_sample"][:];
    temperature = permutedims(nc_pmap["temperature"][:]);
    misfit = permutedims(nc_pmap["misfit"][:]);
    nlayers = permutedims(nc_pmap["nlayers"][:]);
    logppd = permutedims(nc_pmap["logppd"][:]);

    mean_model = nc_pmap["mean_model"][:];
    mode_model = nc_pmap["mode_model"][:];
    p10_model  = nc_pmap["p10_model"][:];
    p50_model  = nc_pmap["p50_model"][:];
    p90_model  = nc_pmap["p90_model"][:];

    ar_valuechange = permutedims(nc_pmap["ar_valuechange"][:]);
    ar_move        = permutedims(nc_pmap["ar_move"][:]);
    ar_birth       = permutedims(nc_pmap["ar_birth"][:]);
    ar_death       = permutedims(nc_pmap["ar_death"][:]);
    swap_histogram = permutedims(nc_pmap["swap_histogram"][:]);
    #read noise vars
    local ar_noisechange,noise_bins,noise_counts,nnoises;
    try
        ar_noisechange = permutedims(nc_pmap["ar_noisechange"][:]);
        noise_bins = nc_pmap["noise_bins"][:];
        noise_counts = nc_pmap["noise_histogram"][:];
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
        ar_nuisancechange = permutedims(nc_pmap["ar_nuisancechange"][:]);
        nuisance_bins = nc_pmap["nuisance_bins"][:];
        nuisance_counts = nc_pmap["nuisance_histogram"][:];
        nnuisances = size(nuisance_counts,2);
    catch
        nnuisances = 0;
        ar_nuisancechange=Array{Float64,2}(undef,0,0);
        nuisance_bins = Array{Float64,2}(undef,0,0);
        nuisance_counts = Array{UInt,2}(undef,0,0);
    end

	fid = nc_pmap.attrib["fiducial"];
	x = nc_pmap.attrib["x"];
	y = nc_pmap.attrib["y"];
	elevation = nc_pmap.attrib["elevation"];

	close(nc_pmap)

    Pmap(depth,value,layer,lchist,cphist,nlhist,observations,errors,ndata,chain,nchains,cvs,
        temperature,misfit,nlayers,logppd,mean_model,mode_model,p10_model,p50_model,p90_model,
        ar_valuechange,ar_move,ar_birth,ar_death,swap_histogram,ar_noisechange,
        noise_bins,noise_counts,nnoises, ar_nuisancechange, nuisance_bins,
        nuisance_counts, nnuisances, fid, x, y, elevation)
end


function ct2cd(cond,thick,depth)
    interfaces = cumsum(thick)
    layerinds = searchsortedfirst.(Ref(interfaces),depth)
    cond[layerinds]
end

function d2t(depth::Array{Float64,1})::Array{Float64,1}
    depth[2:end] - depth[1:end-1]
end

function mean_ct(P::Pmap)::Tuple{Array{Float64,1},Array{Float64,1}}
    P.mean_model, d2t(P.depth)
end

function median_ct(P::Pmap)::Tuple{Array{Float64,1},Array{Float64,1}}
    P.p50_model, d2t(P.depth)
end

function view_rjmcmc_pmap(P::Pmap, TM=Dict())

    p = [
        0.066 0.33  0.775 0.980;
        0.066 0.284 0.550 0.755;
        0.066 0.284 0.325 0.530;
        0.066 0.284 0.100 0.305;

        0.376 0.65 0.82 0.98;
        #0.33  0.57 0.10 0.76;
	0.69 0.92 0.10 0.98 ;
        #0.581 0.65 0.10 0.76;
	0.93 0.99 0.10 0.98;
        ]

    if (P.nnoises > 0 && P.nnuisances > 0)
        #p[6,:] = [0.5 0.87 0.42 0.76];
        #p[7,:] = [0.88 0.99 0.42 0.76];
        p = [p; 0.33 0.65 0.10 0.41; 0.33 0.65 0.45 0.76];
    elseif (P.nnoises > 0 || P.nnuisances > 0)
        #p[6,:] = [0.5 0.87 0.32 0.76];
        #p[7,:] = [0.88 0.99 0.32 0.76];
        p = [p; 0.33 0.65 0.10 0.76];
    end


    fig = figure(figsize=(15,8))

    ax = Array{Any}(undef,9);

    for i=1:size(p,1)
        ap = [p[i,1];p[i,3];p[i,2]-p[i,1];p[i,4]-p[i,3]]
        ax[i] = PyPlot.axes(position=ap)
    end

    cmap = ColorMap("inferno")
    #set zero values of the histogram to black
    #matlab lets you just do this with arrays
    #but PyPlot via julia has to make things hard
    #cmap._i_under = 0
    #cmap.set_under((1,1,1,1))

    sca(ax[6])
    pcolormesh(10 .^ P.value,P.depth,P.lchist,cmap=cmap)
    #plot(10 .^ P.mean_model,P.depth,"-",color="grey", linewidth = 3, label = "mean");
    plot(10 .^ P.p10_model,P.depth,":w", linewidth = 3, label = "p10")
    #plot(10 .^ P.p50_model,P.depth,"-w", linewidth = 3, label = "median")
    plot(10 .^ P.p90_model,P.depth,":w", linewidth = 3, label = "p90")
    if length(TM) > 0
        tmcmap = ct2cd(TM["conductivity"],TM["thickness"],P.depth)
        plot(tmcmap,P.depth,"-g",linewidth = 3, label = "true model")
    end
    ax[6].set_xscale("log")
    ylabel("Depth (m)")
    xlabel("Conductivity (S/m)")
    gca().invert_yaxis()
    legend(loc="lower right", framealpha=0.3)
    #gca().xaxis.tick_top()
    #gca().xaxis.set_label_position("top")

    sca(ax[7])
    plot(P.cphist,P.depth)
    ylim([minimum(P.depth);maximum(P.depth)])
    gca().invert_yaxis()
    xticks([])
    yticks([])
    xlabel("Change histogram")

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
