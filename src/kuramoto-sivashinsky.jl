export KuramotoSivashinsky

"""
ksintegrateDiffEq: integrate kuramoto-sivashinsky equation (Julia)
       u_t = -u*u_x - u_xx - u_xxxx, domain x in [0,Lx], periodic BCs
 inputs
          u = initial condition (vector of u(x) values on uniform gridpoints))
         Lx = domain length
         dt = time step
         Nt = number of integration timesteps
         nsave = save every n-th timestep
 outputs
          u = final state, vector of u(x, Nt*dt) at uniform x gridpoints
This an implementation using ApproxFun and OrdinaryDiffEq.
"""
function KuramotoSivashinsky(u0::AbstractVector{X}, L;
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
    return sol.t, Ti .* sol.u, prob
end

struct KSFunctor{T}; tmp::Vector{T}; end
function (ks::KSFunctor)(du,u,p,t)
    D, T, Ti = p
    mul!(du,Ti,u)
    @. du = -1/2*du^2
    mul!(ks.tmp, T, du)
    mul!(du, D, ks.tmp)
end
