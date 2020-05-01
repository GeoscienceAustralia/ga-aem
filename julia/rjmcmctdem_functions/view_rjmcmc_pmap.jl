using PyPlot

function ct2cd(cond,thick,depth)
    interfaces = cumsum(thickness)
    layerinds = searchsortedfirst.(Ref(interfaces),depth)
    cond[layerinds]
end

function view_rjmcmc_pmap(P, TM=Dict(), noise=true)

    p = [
        0.05 0.50 0.775 0.980;
        0.05 0.50 0.550 0.755;
        0.05 0.50 0.325 0.530;
        0.05 0.50 0.100 0.305;

        0.57 0.99 0.82 0.98;
        0.57 0.94 0.05 0.74;
        0.95 0.99 0.05 0.74;

        ]

    fig = figure(figsize=(10,8))

    ax = Array{Any}(undef,7)

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
    pcolormesh(10 .^ P["value"],P["depth"],P["lchist"],cmap=cmap)
    plot(10 .^ P["mean_model"],P["depth"],"-w");
    plot(10 .^ P["p10_model"],P["depth"],":m")
    plot(10 .^ P["p50_model"],P["depth"],"-m")
    plot(10 .^ P["p90_model"],P["depth"],":m")
    if length(TM) > 0
        tmcmap = ct2cd(TM["conductivity"],TM["thickness"],P["depth"])
        plot(tmcmap,P["depth"],"-g")
    end
    ax[6].set_xscale("log")
    ylabel("Depth (m)")
    xlabel("Conductivity (S/m)")
    gca().invert_yaxis()
    #gca().xaxis.tick_top()
    #gca().xaxis.set_label_position("top")

    sca(ax[7])
    pcolormesh([0 1],P["depth"], reshape(P["cphist"],(size(P["cphist"])...,1)),cmap=cmap)
    gca().invert_yaxis()

    sca(ax[5])
    bar(P["layer"],P["nlhist"],color="r")
    xlabel("Number of layers")
    ylim([0;maximum(P["nlhist"])])

    sca(ax[1]);
    semilogy(P["cvs"],permutedims(P["misfit"]/P["ndata"]));
    ylabel("Normalised misfit");
    ylim([0.09;100]);
    xticks([])

    sca(ax[2]);
    semilogy(P["cvs"],permutedims(P["temperature"]));
    ylim([0.9; maximum(P["temperature"])*1.1]);
    ylabel("Temperature");
    xticks([])

    sca(ax[3]);
    plot(P["cvs"],permutedims(P["nlayers"]))
    ylabel("#Layers");
    xticks([])

    sca(ax[4]);
    plot(P["cvs"],P["ar_valuechange"][1,:],"-m")
    plot(P["cvs"],P["ar_move"][1,:],"-b")
    plot(P["cvs"],P["ar_birth"][1,:],"-g")
    plot(P["cvs"],P["ar_death"][1,:],"-r")
    if noise
        plot(P["cvs"],P["ar_noisechange"][1,:],"-k")
    end


    xlabel("Sample#");
    ylabel("Accept rate");
    ylim([0;100]);
    legend(["change";"move";"birth";"death"],ncol=4)

    gcf()
end

""
