any(pwd() .== LOAD_PATH) || push!(LOAD_PATH, pwd())
using UseGA_AEM, PyPlot

# System geometry
ztxLM    =  35.0 # height of transmitter above ground in m for low moment transmission
ztxHM    =  35.0 # height of transmitter above ground in m for high moment transmission
rrx      = -17.0 # how far behind the loop centre is the receiver in m
dzrxLM   =   2.0 # how far above the transmitter is the receiver for low moment gates
dzrxHM   =   0.2 # how far above the transmitter is the receiver for high moment gates

# make a geometry struct
g = UseGA_AEM.Geometry(ztxLM  = ztxLM,
                       ztxHM  = ztxHM,
                       rrx    = rrx,
                       dzrxLM = dzrxLM,
                       dzrxHM = dzrxHM)

# Initialize an AEM system (Tx/Rx electronics specs hardwired for now)
em = UseGA_AEM.init_GA_AEM()

# Initialize a forward operator with the given geometry
f = UseGA_AEM.EMoperator(em, g)

# Earth model to compute response for using f
conductivity = [0.01, 0.05, 0.01] # in S/m
thickness    = [50.0, 50.0] # in m

# Compute forward for the given system geometry and electronics
f(ztxLM, ztxHM, conductivity, thickness)

# plot
fig, ax = plt.subplots(1,2, figsize=(10,5))
ax[1].step(log10.([conductivity; conductivity[end]]), cumsum([0;thickness;10]))
ax[2].loglog(f.em.tLM, abs.(f.em.SZLM), ".k")
ax[2].loglog(f.em.tHM, abs.(f.em.SZHM), ".k")
ax[1].invert_yaxis()

# make 2nd layer more resistive and compute forward
# conductivity[2] = 0.1
# f(ztxLM, ztxHM, conductivity, thickness)

# plot
# ax[1].step(log10.([conductivity; conductivity[end]]), cumsum([0;thickness;10]), "--")
# ax[2].loglog(f.em.tLM, abs.(f.em.SZLM), "--k")
# ax[2].loglog(f.em.tHM, abs.(f.em.SZHM), "--k")
# ax[1].set_xlabel("log10 σ")
# ax[2].set_xlabel("time s")
# ax[1].set_ylabel("depth m")
# ax[2].set_ylabel("dBz/dt V/(Am^4)")
# ax[2].grid()
# ax[1].set_title("earth model")
# ax[2].set_title("earth response")
fig.tight_layout()
gcf()
