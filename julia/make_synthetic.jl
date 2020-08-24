any((@__DIR__) .== LOAD_PATH) || push!(LOAD_PATH, @__DIR__)
using UseGA_AEM

function forward_model(g::UseGA_AEM.Geometry,
    conductivity::Array{Float64,1},
    thickness::Array{Float64,1};
    stmfile_LM = "",
    stmfile_HM = "")::Array{Float64,2}

    # Initialize an AEM system (Tx/Rx electronics specs hardwired for now)
    stmdict = Dict()
    if (length(stmfile_LM)>0)
        stmdict[:stmfile_LM] = stmfile_LM
    end
    if (length(stmfile_HM)>0)
        stmdict[:stmfile_HM] = stmfile_HM
    end
    em = UseGA_AEM.init_GA_AEM(;stmdict...);
    # Initialize a forward operator with the given geometry
    f = UseGA_AEM.EMoperator(em, g);
    # Compute forward for the given system geometry and electronics
    f(g.geometryLM[1], g.geometryHM[1], conductivity, thickness);

    [f.em.SZLM f.em.tLM; f.em.SZHM f.em.tHM]
end

function write_noisy_forward(g::UseGA_AEM.Geometry,
    conductivity::Array{Float64,1},
    thickness::Array{Float64,1},
    mnoise::Float64,
    fname::String;
    anoise::Array{Float64,1}=[],
    stmfile_LM::String="",
    stmfile_HM::String="")

    data = forward_model(g, conductivity, thickness;
                        stmfile_LM=stmfile_LM, stmfile_HM=stmfile_HM);
    data = data[:,1];

    #sample the noise processes
    length(anoise) > 0 || (anoise = zeros(length(data)))
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
