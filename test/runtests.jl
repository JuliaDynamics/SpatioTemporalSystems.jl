using Test, SpatioTemporalSystems
# Test:

L = 200  # physical size of system
Q = 512  # number of gridpoints
dt = 1/4 # dt (accepted as keyword)
nsave = 1
Nt = 100000
Nskip = 10000
x = L*(0:Q-1)/Q
u0 = cos.(x) + 0.1*sin.(x/8) + 0.01*cos.((2*pi/L)*x);

t,u = ksintegrateDiffEq(u0, L; dt = dt Nt+Nskip, nsave)
t = t[Nskip÷nsave+1:end]
u = u[Nskip÷nsave+1:end]
# path = joinpath(datadir(), "sim","KS_T$(Nt)_L$(L)_Q$(Q).jld2")
# save(path, "t", t, "u", u)
