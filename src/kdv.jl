using ApproxFun, OrdinaryDiffEq
using LinearAlgebra
using DiffEqOperators

# The problem is that kdv needs infinite amount of time to solve...

function makesim(sts::Union{STS{:KortewegDeVries}, STS{:kdv}})
    @unpack ic, T, S, Δt, N = sts
    c = try
            c = sts.p.c
        catch e
            c = 22
        end
    # prepare initial condition
    if ic isa Union{Nothing, Integer}
        # here we integrate a pre-determined initial condition
        # for a random amount of time
        rng = MersenneTwister(isnothing(ic) ? 42 : ic)
        S += Δt*rand(rng, 10:1000)
        x = range(0, 2π, length = N)
        u0 = cos.(x)
    elseif ic isa AbstractVector{<:Real}
        u0 = ic
        @assert length(u0) == N
    else
        error(ICERROR)
    end
    println("i.c. is resolved")

    saveat = S:Δt:(T+S)
    t, u = korteweg_de_vries(u0, saveat)
    t = 0:Δt:T
    return @dict u t sts
end

function korteweg_de_vries(u, saveat)

    n = length(u)
    F = Fourier()
    T = ApproxFun.plan_transform(F, n)
    Ti = ApproxFun.plan_itransform(F, n)

    #Linear Part
    D  = (Derivative(F) → F)[1:n,1:n]
    D3 = (Derivative(F,3) → F)[1:n,1:n]
    A = DiffEqArrayOperator(-D3)

    #Nonlinear Part
    function kdv(dû,û,p,t)
        D,T,Ti,u,tmp = p
        mul!(u,D,û)
        mul!(tmp,Ti,u)
        mul!(u,Ti,û)
        @. tmp=u*tmp
        mul!(u,T,tmp)
        @.dû = 6*u
    end
    û = T*u
    tmp = similar(û)
    p = (D,T,Ti,similar(tmp),tmp)
    prob = SplitODEProblem(A, kdv, û, (0.0, saveat[end]), p)
    println("start")
    sol = solve(prob, RadauIIA5(autodiff=false), dt=Float64(saveat.step), saveat=saveat,
                save_start = false, save_end = false, reltol=1e-14,abstol=1e-14)
    # @assert sol.t == saveat
    return sol.t, map(u -> Ti*u, sol.u)
end
