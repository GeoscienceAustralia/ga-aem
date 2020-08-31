## includes
include("../rjmcmctdem_functions/rjmcmc_pmap_tools.jl")
## list all netcdf files in output directory
pmapdir = "output-nn/pmaps/"
files = readdir(pmapdir)
ncs = filter(x -> endswith(x,".nc"), files)

## centre of bin in which median occurs from a histogram
struct HistStats
    counts::Array{UInt,1}
    bins::Array{Real,1}
end

function histmedian(h::HistStats)::Float64
    nb = h.bins
    nc = cumsum(h.counts)
    med_count = last(nc)/2
    h2mask = nc .>= med_count
    pushfirst!(nc,0)
    pop!(nc)
    h1mask = nc .< med_count
    med_bin = findall(h1mask .& h2mask)[1]
    return nb[med_bin]
end

## additive noise sigmas
add_sigma_LM = [6.4087e+00, 4.0463e+00, 4.2526e+00, 3.0585e+00, 2.3754e+00,
          2.0971e+00, 1.7919e+00, 1.6001e+00, 1.3299e+00, 9.6707e-01,
          7.7027e-01, 7.4084e-01, 8.1148e-01, 6.4811e-01, 6.4059e-01,
          6.3400e-01, 5.9044e-01]
add_sigma_HM = [1.0652e-02, 1.0436e-02, 1.0185e-02, 1.0033e-02, 9.9405e-03,
          9.8320e-03, 9.6964e-03, 9.4945e-03, 9.2693e-03, 8.9339e-03,
          8.5681e-03, 8.2023e-03, 8.4102e-03, 8.9711e-03, 8.3685e-03,
          8.1345e-03, 7.8504e-03, 7.4964e-03, 7.0504e-03, 6.7447e-03,
          6.7160e-03, 5.8913e-03, 5.2938e-03, 4.7887e-03, 4.3874e-03]
## read the files and store sigma/d
nsoundings = length(ncs)
nw_lm = length(add_sigma_LM)
nw_hm = length(add_sigma_HM)
sig_d_lm = zeros(nsoundings, nw_lm)
sig_d_hm = zeros(nsoundings, nw_hm)

for idx in 1:nsoundings
    f = ncs[idx]
    println(f)
    P = read_rjmcmc_pmap(joinpath(pmapdir,f))
    d = P.observations
    nh_lm = HistStats(P.noise_counts[:,1], P.noise_bins[:,1])
    nh_hm = HistStats(P.noise_counts[:,2], P.noise_bins[:,2])
    noise_lm = sqrt.(histmedian(nh_lm)^2 .+ (add_sigma_LM./d[1:nw_lm]).^2)
    noise_hm = sqrt.(histmedian(nh_hm)^2 .+ (add_sigma_HM./d[nw_lm+1:end]).^2)
    sig_d_lm[idx,:] = noise_lm
    sig_d_hm[idx,:] = noise_hm
end

## plot
qranks = [0.05,0.25,0.5,0.75,0.95]
quants_lm = transpose(mapslices(x->quantile(x,qranks),sig_d_lm;dims=1))
quants_hm = transpose(mapslices(x->quantile(x,qranks),sig_d_hm;dims=1))

fig, ax = subplots(1,2, figsize=(20,10))
sca(ax[1])
yscale("log")
ylabel("relative uncertainty")
xlabel("LM gate")
xticks(1:nw_lm)
plot(quants_lm,"o--")
legend("Quantile ".*string.(qranks))
sca(ax[2])
yscale("log")
ylabel("relative uncertainty")
xlabel("HM gate")
xticks(1:nw_hm)
plot(quants_hm,"o--")
gcf()
