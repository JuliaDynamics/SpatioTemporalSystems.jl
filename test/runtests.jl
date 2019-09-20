using Test, SpatioTemporalSystems
# Test:

L = 200  # physical size of system
Q = 512  # number of gridpoints
dt = 1/4 # dt (accepted as keyword)
T = 1100.0
Tskip = 100.0
x = L*(0:Q-1)/Q # position vector
u0 = cos.(x) + 0.1*sin.(x/8) + 0.01*cos.((2*pi/L)*x);

t, u = kuramoto_sivashinsky(u0, L; tspan = (0.0, T), dt = dt, saveat = Tskip:dt:T)

# %%
using PyPlot
ax = gca()
ax.clear()
U = hcat(u...)
ax.imshow(U, extent = (Tskip, T, 0, L))
ax.set_aspect("auto")
ax.set_ylabel("x")
ax.set_xlabel("t")
gcf().tight_layout()
