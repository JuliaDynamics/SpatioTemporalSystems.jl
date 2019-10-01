# SpatioTemporalSystems.jl
This package provides a unified interface for simulating
the temporal evolution of spatiotemporal dynamical systems (both deterministic as well as stochastic).

The provided API interplays excellently with DrWatson, see below.

## API
The API of this package is provided by one struct `SpatioTemporalSystem` and one function `makesim(sts)` with `sts` an instance of `SpatioTemporalSystem`.

### `SpatioTemporalSystem`
This `struct` contains all information necessary to generate spatiotemporal timeseries for a system.

```julia
struct SpatioTemporalSystem{Model}
    N                 # Intr or NTuple{Int}: spatial extent (in pixels!)
    T::Real = 1000    # total time to record timeseries for (in real time)
    Δt::Real = 1      # sampling time in real time units
    S::Real = 10      # Skip time. saving data starts _after_ evolving for S
    ic = nothing      # initial condition: string, seed, array...
    p  = NamedTuple() # parameters of system (expected as NamedTuple)
end
const STS = SpatioTemporalSystem
```

Here `Model` is a `Symbol`, which specifies the "codeword" of the system.
The (currently) available systems are listed at the end.
Most fields of `SpatioTemporalSystem` are self-explaining. Notice that
`N` represents the spatial extend of the system in pixels (the array size), and not
the physical size of the system (the physical size, if needed, can be
given as part of the parameters).

All parameters of the system are contained in the field `p`, which is a dictionary that is expanded as keyword arguments to the low level functions that evolve the system.

The field `ic` needs some discussion. If `ic` is `nothing` (default), each system uses
its own, dedicated default initial condition (always the same). Else,
`ic` can be an initial condition, for example if the system has one field,
then `ic` will be the array `u0`, but if there are more fields then `ic` will
be a Tuple of arrays.

The important thing is that `ic` can also be an integer, which represents a "seed".
Then `ic` is used to seed a `MersenneTwister` generator, which will generate the
initial condition (either randomly, or by evolving a pre-defined initial condition
for a random amount of time).
`ic` could even be a `String`, if the underlying low-level functions for
a specific system provide such option.

`SpatioTemporalSystem` fully customizes DrWatson's `savename` (`S` is not printed in `savename`), and has been appropriately configured to play well with `produce_or_load` via the provided function `makesim`.
This ensures that you don't produce identical spatiotemporal timeseries more than once.

### `makesim`

Evolution of a spatiotemporal system is done with `makesim(sts)`,
which itself dispatches on other low-level functions based on the `Model` (a `Symbol` codeword):

```julia
makesim(sts::SpatioTemporalSystem) -> dict
```
Make a simulation of a spatiotemporal system. Return a **dictionary**
that contains the timevector `t`, corresponding evolved fields
and the system itself `sts`. Each system saves its fields with
different names.


## Available Systems
`makesim` dispatches on the `Model` that parameterizes the `SpatioTemporalSystem`. All systems have a long name stated in camel case, as well as a shorthand name in small letters.
Each system also has a low-level function that explicitly describes the parameters of the system and other information.

Currently available systems are:

* `:Barkley, :bk` : see the `bakley` function.
* `:KuramotoSivashinsky, :ksiva` : see the `kuramoto_sivashinsky` function.
* `:CoupledHenon1D, :henon1d`, see the `henon1d` function.
