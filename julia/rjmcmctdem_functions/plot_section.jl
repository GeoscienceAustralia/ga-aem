using PyPlot

function plot_section(line_nc_dir::String)
    """
    plot_section(line_nc_dir)

    process a directory of netCDFs into a tuple of arrays that can be used
    to plot an rjmcmc section.

    ## Parameters:
    line_nc_dir: string location of directory, which should contain netCDF files
    with the garjmcmc outputs from soundings on a single line.

    ## Returns:
    All arrays are of dimension (number of soundings x number of depth
    cells) and are sorted along the first axis according to the easting in
    ascending order

    plotting_easting: spatial easting coordinate corresponding to each sounding
    in the directory. Repeated for each depth cell for compatibility with
    matplotlib routines that need (x,y,z) values for all cells.

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

    for (sounding, ncfile) in enumerate(ncs)
        println("Reading sounding "*string(sounding)*" of "*string(nsoundings));
        P = read_rjmcmc_pmap(ncfile);
        elevations[sounding] = P.elevation;
        eastings[sounding] = P.x;
        northings[sounding] = P.y;

        median_models[sounding,:] = P.p50_model;
        spread[sounding,:] = P.p90_model - P.p10_model;
    end

    #reorder so plotting proceeds smoothly
    idxs = sortperm(eastings);
    elevations = elevations[idxs];
    eastings = eastings[idxs];
    northings = northings[idxs];
    median_models = median_models[idxs,:];

    #model cell elevations rel sea level
    plotting_elevation = transpose(depth) .- elevations;
    plotting_easting = repeat(eastings,1,ndepthcells);

    plotting_easting, plotting_elevation, median_models, spread
end

function nuisance_sections(line_nc_dir::String)

end
