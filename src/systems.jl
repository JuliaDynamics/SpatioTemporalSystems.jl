using DrWatson

export SpatioTemporalSystem, STS, makesim

@with_kw struct SpatioTemporalSystem{Model}
    N                 # Intr or NTuple{Int}: spatial extent (in pixels!)
    T::Real = 100     # total time to record timeseries for (in real time)
    Δt::Real = 1      # sampling time in real time units
    S::Real = 10      # Skip time. saving data starts _after_ evolving for S
    ic = nothing      # initial condition: string, seed, array...
    p  = NamedTuple() # parameters of system (expected as NamedTuple)
end
const STS = SpatioTemporalSystem

# Extend savename:
DrWatson.allaccess(::STS) = (:N, :T, :Δt, :ic, :p)
DrWatson.default_prefix(c::STS{M}) where {M} = string(M)
DrWatson.default_allowed(c::STS) = (Real, Tuple, NamedTuple, Dict, String)
DrWatson.default_expand(c::STS) = ["p"]

"""
    makesim(sts::SpatioTemporalSystem) -> dict
Make a simulation of a spatiotemporal system. Return a **dictionary**
that contains the timevector `t`, corresponding evolved fields
and the system itself `sts`. Each system saves its fields with
different names.
"""
function makesim end
