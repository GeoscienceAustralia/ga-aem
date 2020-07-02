any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using UseGA_AEM

function forward_model(g::UseGA_AEM.Geometry,
    conductivity::Array{Float64,1},
    thickness::Array{Float64,1})::Array{Float64,1}

    # Initialize an AEM system (Tx/Rx electronics specs hardwired for now)
    em = UseGA_AEM.init_GA_AEM();
    # Initialize a forward operator with the given geometry
    f = UseGA_AEM.EMoperator(em, g);
    # Compute forward for the given system geometry and electronics
    f(ztxLM, ztxHM, conductivity, thickness);

    [f.em.SZLM; f.em.SZHM]
end

function write_noisy_forward(g::UseGA_AEM.Geometry,
    conductivity::Array{Float64,1},
    thickness::Array{Float64,1},
    mnoise::Float64,
    anoise::Array{Float64,1},
    fname::String)

    data = forward_model(g, conductivity, thickness);

    #sample the noise processes
    mnoise_r = 1. .+ mnoise .* randn(length(data));
    anoise_r = anoise .* randn(length(data));

    data = mnoise_r .* data .+ anoise_r;

    file = open(fname, "w");
    for d in data
        write(file, string(d, ' '));
    end

    #write the noise in as well
    for a in anoise
        write(file, string(a, ' '));
    end

    close(file)

end
