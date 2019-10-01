using DrWatson

export SpatioTemporalSystem, STS, makesim

@with_kw struct SpatioTemporalSystem{Model}
    N                 # Intr or NTuple{Int}: spatial extent (in pixels!)
    T::Real = 1000    # total time to evolve for (in real time)
    Δt::Real = 1      # sampling time in real time units
    ic = nothing      # initial condition: string, seed, array...
    p  = NamedTuple() # parameters of system (expected as NamedTuple)
end
const STS = SpatioTemporalSystem

# Extend savename:
DrWatson.default_prefix(c::STS{M}) where {M} = string(M)
DrWatson.default_allowed(c::STS) = (Real, Tuple, NamedTuple, Dict, String)
DrWatson.default_expand(c::STS) = ["p"]
#TODO: does the expand also need "N" ?

"""
    makesim(sts::SpatioTemporalSystem) -> dict
Make a simulation of a spatiotemporal system. Return a **dictionary**
that contains the timevector `t`, corresponding evolved fields
and the system itself `sts`. Each system saves its fields with
different names.
"""
function makesim end



function
    ret = Dict{String, Any}("STS" => c)
    if c.model == "barkley"
        seed = c.ic == nothing ? rand(1:1000) : c.ic
        Nx, Ny = typeof(c.N) == Tuple ? c.N : (c.N, c.N)
        t, U, V = barkley(c.T, c.Δt, get(c.p, :periodic, true);
                       seed = seed, Nx = Nx, Ny = Ny, c.p...)
        ret["U"] = U; ret["V"] = V; ret["t"] = t
    elseif c.model == "henon1D"
        seed = c.ic == nothing ? rand(1:1000) : c.ic
        U, V = coupled_henon1D(c.N, c.T; seed = seed, c.p...)
        ret["U"] = U; ret["V"] = V
    elseif c.model == "ksiva"
        # What is the initial condition here?
        c.ic = nothing # keep nothing until we implement initial condition
    else
        error("$(c.model) is an unknown model")
    end
    return ret
end
