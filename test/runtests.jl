using Test, SpatioTemporalSystems, DrWatson
# Test:

sts = STS{:bk}(N = 100, T = 1000, Δt = 0.1)
d = makesim(sts)

U = d["U"]

# using PyPlot
# figure()
# imshow(U[16])

# %% Henon map
sts = STS{:henon1d}(ic = 4, N = 200, T = 400, S = 1000)
d = makesim(sts)
u = d["U"]
#
# using PyPlot
# figure()
# U = hcat(u...)
# imshow(U)

# %% Kuramoto
L = 200  # physical size of system
N = 512  # number of gridpoints
Δt = 1/4
T = 100.0
x = L*(0:N-1)/N # position vector
u0 = cos.(x) + 0.1*sin.(x/8) + 0.01*cos.((2*pi/L)*x);

sts = STS{:ksiva}(N = N, ic = u0, T = T, Δt = Δt, p = (L = L,))
savename(sts)

sim = makesim(sts)
u = sim["u"]
# @unpack u, t = sim

# using PyPlot
# ax = gca()
# ax.clear()
# U = hcat(u...)
# ax.imshow(U, extent = (0, T, 0, L), cmap = "inferno")
# ax.set_aspect("auto")
# ax.set_ylabel("x")
# ax.set_xlabel("t")
# gcf().tight_layout()
