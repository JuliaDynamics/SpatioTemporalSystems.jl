using Test, SpatioTemporalSystems, DrWatson
# Test:

sts = STS{:bocf}(N = 100, T = 10, Δt = 0.1, S = 0, ic = "spiral")
d = makesim(sts)

U = d["U"]

using PyPlot
close("all")
fig, axs = subplots(1, 3, figsize = (10, 4))

for i in 1:3
    axs[i].imshow(U[i])
end
