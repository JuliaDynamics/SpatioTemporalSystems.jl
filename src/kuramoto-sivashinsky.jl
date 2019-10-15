using ApproxFun, OrdinaryDiffEq
using LinearAlgebra
using DiffEqOperators

function makesim(sts::Union{STS{:KuramotoSivashinsky}, STS{:ksiva}})
    @unpack ic, T, S, Δt, N = sts
    L = try
            L = sts.p.L
        catch e
            L = 22
        end
    # prepare initial condition
    if ic isa Union{Nothing, Integer}
        # here we integrate a pre-determined initial condition
        # for a random amount of time
        rng = MersenneTwister(isnothing(ic) ? 42 : ic)
        S += Δt*rand(rng, 10:1000)
        x = L*(0:N-1)/N
        u0 = @. cos(x) + 0.1*sin(x/8) + 0.01*cos((2π/L)*x)
    elseif ic isa AbstractVector{<:Real}
        u0 = ic
        @assert length(u0) == N
    else
        error(ICERROR)
    end
    saveat = S:Δt:(T+S)
    t, u = kuramoto_sivashinsky(u0, L, saveat)
    t = 0:Δt:T
    return @strdict u t sts
end


"""
    kuramoto_sivashinsky(u, L, saveat)

Solve the Kuramoto Sivashinsky system given by the PDE
```math
u_t = -u*u_x - u_{xx} - u_{xxxx}
```
with periodic boundary conditions and initial condition `u`.

`L` is the physical spatial length (with `length(u)` points sampled
equidistantly in it).

`saveat` dictates at which timepoints to save the system (which is solved
using `OrdinaryDiffEq`).
"""
function kuramoto_sivashinsky(u, L, saveat)
    n = length(u)                  # number of gridpoints
    F = Fourier(0..L)
    T = ApproxFun.plan_transform(F, n)
    Ti = ApproxFun.plan_itransform(F, n)

    #Linear Part
    D  = (Derivative(F) → F)[1:n,1:n]
    D2 = Derivative(F,2)[1:n,1:n]
    D4 = Derivative(F,4)[1:n,1:n]
    A = DiffEqArrayOperator(Diagonal(-D2-D4))

    #Nonlinear Part
    ks = KSFunctor(zeros(n))

    params = (D, T, Ti)
    prob = SplitODEProblem(A, ks, T*u, (0.0, saveat[end]), params)
    sol = solve(prob, ETDRK4(), dt=Float64(saveat.step), saveat=saveat,
                save_start = false, save_end = false)
    # @assert sol.t == saveat
    return sol.t, map(u -> Ti*u, sol.u)
end

struct KSFunctor{T}; tmp::Vector{T}; end
function (ks::KSFunctor)(du,u,p,t)
    D, T, Ti = p
    mul!(du,Ti,u)
    @. du = -1/2*du^2
    mul!(ks.tmp, T, du)
    mul!(du, D, ks.tmp)
end
