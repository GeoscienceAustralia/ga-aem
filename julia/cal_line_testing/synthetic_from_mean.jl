## includes for rjmcmc tools scripts
include("../rjmcmctdem_functions/rjmcmc_pmap_tools.jl")
include("../make_synthetic.jl")
## load Julia interace for GA-AEM forward modelling
any(abspath("../") .== LOAD_PATH) || push!(LOAD_PATH, abspath("../"))
using UseGA_AEM
## read pmap and make the mean model
P = read_rjmcmc_pmap("seq.00001501.200311.610684.900000.nc")
c, t = median_ct(P)
c = 10 .^ c
## construct the AEM system object with UseGA_AEM
em = UseGA_AEM.init_GA_AEM(;stmfile_LM = "Skytem-LM-Menindee-Yancowinna.stm",
                           stmfile_HM = "Skytem-HM-Menindee-Yancowinna.stm")
## read tx height from file
line = 1501
lines = readlines("em_with_noise_calibration_lines.dat")
data = split(lines[line])

## calculate tx height and display (9th column)
tx_height = parse(Float64,data[9])
## make a system geometry object to generate the synthetic with
g = UseGA_AEM.Geometry(ztxHM = tx_height, ztxLM = tx_height, rrx = -17.0, dzrxLM = 2.0, dzrxHM = 0.2)
## collect additive noise from file (columns 110-151)
anoise = parse.(Float64, data[110:151])
## do the deterministic forward model (to check whether the additive noise makes sense)
synth = forward_model(g,c,t;stmfile_LM = "Skytem-LM-Menindee-Yancowinna.stm",
                           stmfile_HM = "Skytem-HM-Menindee-Yancowinna.stm")
## write the noisy forward
write_noisy_forward(g,c,t,0.03,"meanmodel_synthetic.dat";anoise=anoise,
                    stmfile_LM = "Skytem-LM-Menindee-Yancowinna.stm",
                    stmfile_HM = "Skytem-HM-Menindee-Yancowinna.stm")

# data file columns: 1-17 data LM, 18-42 data HM, 43-59 anoise LM, 60-84 anoise HM,
# 85 tx_height
## our inversion also uses the transmitter height, so write this back too
f = open("meanmodel_synthetic.dat","a")
write(f, string(tx_height))
close(f)
##
signal = parse.(Float64,data[26:67])
signal[1]/anoise[18]
anoise[18]
