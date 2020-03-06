using PyPlot

function view_rjmcmc_pmap(P)

    p = [
        0.05 0.50 0.775 0.980;
        0.05 0.50 0.550 0.755;
        0.05 0.50 0.325 0.530;
        0.05 0.50 0.100 0.305;

        0.55 0.99 0.82 0.98;
        0.55 0.94 0.05 0.70;
        0.95 0.99 0.05 0.70;

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
    pcolor(P["value"],P["depth"],P["lchist"],cmap=cmap)
    plot(P["mean_model"],P["depth"],"-w");
    plot(P["p10_model"],P["depth"],":m")
    plot(P["p50_model"],P["depth"],"-m")
    plot(P["p90_model"],P["depth"],":m")
    ylabel("Depth (m)")
    xlabel("LOG10(Conductivity (S/m))")
    gca().invert_yaxis()
    #gca().xaxis.tick_top()
    #gca().xaxis.set_label_position("top")

    sca(ax[7])
    pcolor([0 1],P["depth"], reshape(P["cphist"],(size(P["cphist"])...,1)),cmap=cmap)
    gca().invert_yaxis()

    sca(ax[5])
    bar(P["layer"],P["nlhist"],color="r")
    xlabel("Number of layers")
    ylim([0;maximum(P["nlhist"])])

    sca(ax[1]);
    semilogy(P["cvs"],permutedims(P["misfit"]/P["ndata"]));
    ylabel("Normalised misfit");
    ylim([0.09;100]);

    sca(ax[2]);
    semilogy(P["cvs"],permutedims(P["temperature"]));
    ylim([0.9; maximum(P["temperature"])*1.1]);
    ylabel("Temperature");

    sca(ax[3]);
    plot(P["cvs"],permutedims(P["nlayers"]))
    ylabel("#Layers");

    sca(ax[4]);
    plot(P["cvs"],P["ar_valuechange"][1,:],"-m")
    plot(P["cvs"],P["ar_move"][1,:],"-b")
    plot(P["cvs"],P["ar_birth"][1,:],"-g")
    plot(P["cvs"],P["ar_death"][1,:],"-r")


    xlabel("Sample#");
    ylabel("Accept rate");
    ylim([0;100]);
    legend(["change";"move";"birth";"death"],ncol=4)

    gcf()
end

""