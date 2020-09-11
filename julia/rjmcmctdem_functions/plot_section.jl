using PyPlot

function section_arrays(line_nc_dir::String)
    """
    section_arrays(line_nc_dir)

    process a directory of netCDFs into a tuple of arrays that can be used
    to plot an rjmcmc section.

    ## Parameters:
    line_nc_dir: string location of directory, which should contain netCDF files
    with the garjmcmc outputs from soundings on a single line.

    ## Returns:
    All arrays are of dimension (number of soundings x number of depth
    cells) and are sorted along the first axis according to the easting in
    ascending order

    plotting_dist: spatial coordinate corresponding to each sounding
    in the directory. Repeated for each depth cell for compatibility with
    matplotlib routines that need (x,y,z) values for all cells. calculated
    according to distance from the last (maximum) fiducial in the line.

    plotting_elevation: map between depth cells and altitude, so that the
    section can accurately represent topography.

    median_models: the median conductivity for each depth cell in each sounding

    spread: difference between 10th and 90th percentile conductivity in each
    depth cell
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

    maxfid = P0.fiducial;
    maxsdg = 1;

    for (sounding, ncfile) in enumerate(ncs)
        println("Reading sounding "*string(sounding)*" of "*string(nsoundings));
        P = read_rjmcmc_pmap(ncfile);
        elevations[sounding] = P.elevation;
        eastings[sounding] = P.x;
        northings[sounding] = P.y;

        if P.fiducial > maxsdg
            maxfid = P.fiducial;
            maxsdg = sounding;
        end

        median_models[sounding,:] = P.p50_model;
        spread[sounding,:] = P.p90_model - P.p10_model;
    end

    lineDist = sqrt.((eastings .- eastings[maxsdg]).^2 .+ (northings .- northings[maxsdg]).^2)

    #reorder so plotting proceeds smoothly
    idxs = sortperm(lineDist);
    elevations = elevations[idxs];
    eastings = eastings[idxs];
    northings = northings[idxs];
    median_models = median_models[idxs,:];
    spread = spread[idxs,:];

    #model cell elevations rel sea level
    plotting_elevation = transpose(depth) .- elevations;
    plotting_dist = repeat(lineDist[idxs],1,ndepthcells);

    plotting_dist, plotting_elevation, median_models, spread
end

function nuisance_sections(line_nc_dir::String)
    """
        nuisance_sections(line_nc_dir)

        given a directory with netCDF files for an AEM line,
        returns a tuple of arrays with all the data needed
        to plot a "section" of the inverted transmitter height
        for that line.

        Parameters:
        line_nc_dir: string location of directory containing pmap netcdfs
        for the line you want a section plotted for.

        Returns:
        All arrays of shape (number of soundings x number of tx_height bins)
        nuisance_counts: histogram counts for tx height, returned per sounding.
        nuisance_bins: histogram bins for tx height, per sounding.
        line_dist: distance along line per sounding (x axis for plotting
        results).

    """

end

function plot_section(x::Array{Float64,2}, y::Array{Float64,2},
                      z::Array{Float64,2}, dz::Array{Float64,2};
                      lb_dz::Real = 1.0,
                      ub_dz::Real = 3.0)
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

    Returns:
    fig - handle to a PyPlot figure containing the resultant plot. Colours map
          to conductivity and alpha values to uncertainty.
    """
    #long
    fig, ax = plt.subplots(1,1,figsize=(30,5));
    alpha = (dz .- lb_dz)./(ub_dz - lb_dz);
    alpha = min.(1,max.(0,alpha));
    alphacmap = ColorMap("alphacmap",[(0.0,1.0,1.0),(1.0,1.0,1.0)],
                            [(0.0,1.0,1.0),(1.0,1.0,1.0)],
                            [(0.0,1.0,1.0),(1.0,1.0,1.0)],
                            [(0.0,0.0,0.0),(1.0,1.0,1.0)]);

    # second pcolormesh masks the first with white.
    # pcolormesh does not accept non-scalar alpha
    # channels or multivariate colour mapping,
    # so this is the easiest way to achieve the
    # masking effect
    ax.pcolormesh(x,y,z,cmap="jet");
    ax.pcolormesh(x,y,alpha,cmap=alphacmap);

    #ground line
    plot(x,y[:,1],"-k");
    ylim(minimum(y),minimum(maximum(y,dims=2)));
    ax.invert_yaxis();

    fig

end
