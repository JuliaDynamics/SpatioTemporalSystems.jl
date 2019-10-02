module SpatioTemporalSystems

using Parameters
const ICERROR = "Incorrect type for `ic`."

include("systems.jl")
include("barkley.jl")
include("henon.jl")
include("kuramoto-sivashinsky.jl")
include("kdv.jl")

end
