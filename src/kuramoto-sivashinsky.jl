export kuramoto_sivashinsky

"""
```
kuramoto_sivashinsky(u0::AbstractVector{X}, L::Real;
solver = ETDRK4(), kwargs...)
```
"""
function kuramoto_sivashinsky(u0::AbstractVector{X}, L::Real;
    solver = ETDRK4(), tspan = (zero(T), one(T)), kwargs...) where {X}
    n = length(u0)                  # number of gridpoints
    S = Fourier(0..L)
    T = ApproxFun.plan_transform(S, n)
    Ti = ApproxFun.plan_itransform(S, n)

    #Linear Part
    D  = (Derivative(S) → S)[1:n,1:n]
    D2 = Derivative(S,2)[1:n,1:n]
    D4 = Derivative(S,4)[1:n,1:n]
    A = DiffEqArrayOperator(Diagonal(-D2-D4))
    #Nonlinear Part
    ksnonlinear = KSFunctor(zeros(X, n))

    params = (D, T, Ti)
    prob = SplitODEProblem(A, ksnonlinear, T*u0, tspan, params)
    sol = solve(prob, solver; kwargs...)
    return sol.t, map(u -> Ti*u, sol.u), prob
end

struct KSFunctor{T}; tmp::Vector{T}; end
function (ks::KSFunctor)(du,u,p,t)
    D, T, Ti = p
    mul!(du,Ti,u)
    @. du = -1/2*du^2
    mul!(ks.tmp, T, du)
    mul!(du, D, ks.tmp)
end
