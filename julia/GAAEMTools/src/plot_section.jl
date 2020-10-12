using PyPlot

function sorted_pmaps(line_nc_dir::String)::Array{Pmap,1}
    """
    sorted_pmaps(line_nc_dir)

    read all Pmap netCDF files in the given directory and return an array of Pmap
    objects, sorted according to distance along the line. Directory is assumed to
    contain only one line of Pmap netCDFs.

    ## Parameters:
    line_nc_dir: string, location of directory containing pmaps.

    ## Returns:
    pmaps: sorted array of Pmap objects.

    """

    #make the list of netCDFs
    files = readdir(line_nc_dir);
    ncs = filter(x -> endswith(x,".nc"), files);

    ncs = joinpath.(line_nc_dir,ncs);

    if length(ncs) < 2
        error("need at least 2 soundings to plot a section.")
    end

    #read everything
    pmaps = map(ncfile -> read_rjmcmc_pmap(ncfile), ncs);

    #sorting by "distance along line" is a little confusing.
    #first we see whether the data is spread more along a north-south
    #or east-west axis. Then we choose the most northerly or most westerly
    #sounding (on the most spread-out coordinate), and sort according to
    #distance from this sounding

    north = map(pmap -> pmap.y, pmaps);
    east = map(pmap -> pmap.x, pmaps);

    minsdg = maximum(north) - minimum(north) > maximum(east) - minimum(east) ? argmax(north) : argmin(east);

    fcoords = [pmaps[minsdg].y, pmaps[minsdg].x];
    sortfunc = pmap -> (pmap.y - fcoords[1])^2 + (pmap.x - fcoords[2])^2;

    sort!(pmaps, by=sortfunc);

    return pmaps
end

function section_arrays(line_nc_dir::String)
    """
    section_arrays(line_nc_dir)

    process a directory of netCDFs into a tuple of arrays that can be used
    to plot an rjmcmc section.

    ## Parameters:
    line_nc_dir: string location of directory, which should contain netCDF files
    with the garjmcmc outputs from soundings on a single line.

    ## Returns:

    The first array is 1D, of size (number of soundings).

    line_dist: spatial coordinate corresponding to each sounding
    in the directory. Calculated
    according to distance from the last (maximum) fiducial in the line.

    The next set of arrays are of dimension (number of soundings x number of depth
    cells) and are sorted along the first axis according to line_dist in
    ascending order.

    plotting_elevation: map between depth cells and altitude, so that the
    section can accurately represent topography.

    median_models: the median conductivity for each depth cell in each sounding

    spread: difference between 10th and 90th percentile conductivity in each
    depth cell

    The last array is of dimension (2 x number of soundings x number of
    tx height histogram cells) and enables plotting of inverted tx height
    posterior histograms along the line.

    nuisances: tx_height histogram for each sounding. nuisances[1,:,:] are the
    bin centres and nuisances[2,:,:] are the counts.

    """
    #make the list of netCDFs
    files = readdir(line_nc_dir);
    ncs = filter(x -> endswith(x,".nc"), files);

    ncs = joinpath.(line_nc_dir,ncs);


    #read each nc and set data arrays for plotting
    nsoundings = length(ncs);

    if nsoundings < 2
        error("need at least 2 soundings to plot a section.")
    end

    eastings = zeros(nsoundings);
    northings = zeros(nsoundings);
    elevations = zeros(nsoundings);

    P0 = read_rjmcmc_pmap(ncs[1]);
    depth = P0.depth;
    ndepthcells = length(depth);

    median_models = zeros(nsoundings,ndepthcells);
    spread = zeros(nsoundings, ndepthcells);

    nbins = length(P0.nuisance_bins);
    nuisances = zeros(2,nsoundings, nbins);

    maxfid = P0.fiducial;
    maxsdg = 1;

    for (sounding, ncfile) in enumerate(ncs)
        println("Reading sounding "*string(sounding)*" of "*string(nsoundings));
        P = read_rjmcmc_pmap(ncfile);
        elevations[sounding] = P.elevation;
        eastings[sounding] = P.x;
        northings[sounding] = P.y;

        if P.fiducial > maxfid
            maxfid = P.fiducial;
            maxsdg = sounding;
        end

        median_models[sounding,:] = P.p50_model;
        spread[sounding,:] = P.p90_model - P.p10_model;

        nuisances[2,sounding,:] = P.nuisance_counts;
        nuisances[1,sounding,:] = P.nuisance_bins;
    end

    line_dist = sqrt.((eastings .- eastings[maxsdg]).^2 .+ (northings .- northings[maxsdg]).^2)

    #reorder so plotting proceeds smoothly
    idxs = sortperm(line_dist);
    line_dist = line_dist[idxs];
    elevations = elevations[idxs];
    median_models = median_models[idxs,:];
    spread = spread[idxs,:];
    nuisances = nuisances[:,idxs,:];

    #model cell elevations rel sea level
    plotting_elevation = transpose(depth) .- elevations;

    line_dist, plotting_elevation, median_models, spread, nuisances
end

function noise_hists()

end

function nuisance_section(x::Array{Float64,1}, yz::Array{Float64,3}; ax=nothing)
    """
    nuisance_section(x,yz)

    Plot a section of the inverted transmitter height using pcolormesh.

    Parameters:
    x: spatial plotting coordinate for each sounding, in sorted order.
    yz: yz[1,:,:] are the histogram bins for the inverted tx height.
        yz[2,:,:] are the histogram counts for the inverted tx height.

    ax: axes object for the plot. if none, new axes are created.

    Returns:

    fig: handle to the resultant plot. is nothing if axes are passed.

    """
    #compatible shapes for pcolormesh
    x = repeat(x,1,size(yz,3));

    fig=nothing
    if isnothing(ax)
        fig, ax = subplots(1,1,figsize=(30,5));
    end
    ax.pcolormesh(x,yz[1,:,:],yz[2,:,:],cmap="viridis");
    fig
end

function plot_section(x::Array{Float64,1}, y::Array{Float64,2},
                      z::Array{Float64,2}, dz::Array{Float64,2};
                      lb_dz::Real = 1.0,
                      ub_dz::Real = 3.0,
                      ax=nothing)
    """
    plot_section(x,y,z,dz)

    Plot a section of output, given arrays representing sounding positions,
    cell elevations, conductivity and uncertainty in the conductivity.

    Parameters:
    x, y, z, dz are arrays of shape (number of soundings x number of depth cells).
    x - lateral position of each cell.
    y - elevation of each cell (assumed relative to sea level)
    z - inverted conductivity
    dz - uncertainty in inverted conductivity

    Optional:
    lb_dz: lower bound for alpha blending of section according to posterior
           spread.
    ub_dz: upper bound for same.
    ax: axes object for the plot. If none, new axes are created.

    Returns:
    fig - handle to a PyPlot figure containing the resultant plot. Colours map
          to conductivity and alpha values to uncertainty. Is nothing if axes
          are passed.
    """
    fig=nothing
    #long figure
    if isnothing(ax)
        fig, ax = plt.subplots(1,1,figsize=(30,5));
    end
    alpha = (dz .- lb_dz)./(ub_dz - lb_dz);
    alpha = min.(1,max.(0,alpha));
    alphacmap = ColorMap("alphacmap",[(0.0,1.0,1.0),(1.0,1.0,1.0)],
                            [(0.0,1.0,1.0),(1.0,1.0,1.0)],
                            [(0.0,1.0,1.0),(1.0,1.0,1.0)],
                            [(0.0,0.0,0.0),(1.0,1.0,1.0)]);

    # repeat x coord along second dimension of
    # y coord so their shapes match
    x = repeat(x, 1, size(y,2));

    # second pcolormesh masks the first with white.
    # pcolormesh does not accept non-scalar alpha
    # channels or multivariate colour mapping,
    # so this is the easiest way to achieve the
    # masking effect
    ax.pcolormesh(x,y,z,cmap="jet",vmin=log10(0.002),vmax=log10(5.0),shading="gouraud");
    ax.pcolormesh(x,y,alpha,cmap=alphacmap);

    #ground line
    plot(x,y[:,1],"-k");
    ylim(minimum(y),minimum(maximum(y,dims=2)));
    ax.invert_yaxis();

    fig

end

function plot_section_with_nuisance(x,y,cond,spread,nuisance)
    """
    plot_section_with_nuisance(x,y,cond,spread,nuisance)

    Plot a section with the inverted transmitter height on top.
    Parameters:
    x: plotting spatial coordinate for each sounding
    y: plotting elevation
    cond: conductivity model
    spread: 80% credible interval size (10-90 %ile) fot each sounding
    nuisance: nuisance histogram per sounding

    Returns:
    fig: the resulting figure with nuisances plotted above the section
    """

    fig, ax = plt.subplots(2,1,figsize=(30,10));

    plot_section(x,y,cond,spread,ax=ax[2]);
    nuisance_section(x,nuisance,ax=ax[1]);

    ax[1].set_ylabel("Transmitter height");
    ax[2].set_ylabel("Depth (mAHD)");
    ax[2].set_xlabel("Distance along line (m)");

    fig

end
