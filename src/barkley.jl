using Random

function makesim(sts::Union{STS{:Barkley}, STS{:bk}})
    @unpack N, T, Δt, S, ic = sts
    Nx, Ny = typeof(N) <: Tuple ? (N[1], N[2]) : (N, N)
    FType = get(sts.p, :FType, Float32)
    # Assertions
    @assert Δt/BK_INTEG_DT == round(Int, Δt/BK_INTEG_DT) "Δt must be multiple of the integration step (0.01)"
    every = round(Int, Δt/BK_INTEG_DT)
    @assert Nx ≥ 40
    @assert Ny ≥ 40
    @assert Nx%10==0
    @assert Ny%10==0
    # prepare initial condition...
    if ic isa Union{Nothing, Integer}
        rng = MersenneTwister(isnothing(ic) ? 42 : ic)
        αα, ββ = 10, 10
        init = rand(rng,αα,ββ)
        u = zeros(FType,Nx,Ny)
        v = zeros(FType,Nx,Ny)
        v .= u .= repeat(init, inner=(Nx÷αα,Ny÷ββ))
        v .= v .> 0.2
    elseif ic isa Tuple
        u, v = ic
    else
        error(ICERROR)
    end
    t, U, V = barkley(u, v, S, T, Δt; sts.p...)
    return @strdict t sts U V
end

"""
    barkley(U0, V0, S, T, Δt; kwargs...)
Simulate the nonlinear Barkley model with initial conditions `U0, V0` for
the fields `U, V`.

## Keywords
```
periodic = true, FType = Float32,
skip=100, a=0.75, b=0.06, ε=0.08, D=1/50, h=0.1
```
The simulation uses `periodic` boundary conditions by default, or absorbing otherwise.
`skip` is an integer, multiple of `Δt`, which notes how many
`Δt` frames to skip before starting to save data. `FType` is the
number type of the output fields.
`a,b,ε,D` are formal parameters of the system while `h` is the spatial sampling
(the evolution algorithm is the one from the Scholarpedia entry:
http://www.scholarpedia.org/article/Barkley_model )
(notice that the integration time is fixed to 0.01, irrespectively of `Δt`).

`a` can also be a `Matrix` (of size `Nx×Ny`), to e.g. represent inhomogeneity.
"""
function barkley(u::Matrix, v::Matrix, S, T, Δt; periodic=true,
                 skip=100, a=0.75, b=0.06, ε=0.08, D=1/50, h=0.1,
                 )

    FType = eltype(u)
    bk_integ_dt = FType(BK_INTEG_DT)
    total = Int(T÷Δt) # how many steps to save
    Nx, Ny = size(u)
    @assert size(u) == size(v)
    u, v = copy(u), copy(v)
    U = Matrix{FType}[]
    V = Matrix{FType}[]
    sizehint!(U, total); sizehint!(V, total)

    Σ = zeros(FType, Nx, Ny, 2)
    r, s = 1, 2
    every = Δt÷bk_integ_dt

    function F(u, uth)
        if u < uth
            u/(1-(bk_integ_dt/ε)*(1-u)*(u-uth))
        else
            (u + (bk_integ_dt/ε)*u*(u-uth))/(1+(bk_integ_dt/ε)*u*(u-uth))
        end
    end

    function periodic_step!(u,v,Σ,a,b,Nx,Ny,s,r,D,bk_integ_dt,h)
        @inbounds for i=1:Nx, j=1:Ny
            uth = (v[i,j] + b)/_get_a(a, i, j)
            v[i,j] = v[i,j] + bk_integ_dt*(u[i,j]^3 - v[i,j])
            u[i,j] = F(u[i,j], uth) + D*bk_integ_dt/h^2 *Σ[i,j,r]
            Σ[i,j,s] -= 4u[i,j]
            Σ[  mod(i-1-1,Nx)+1,j,s] += u[i,j]
            Σ[  mod(i+1-1,Nx)+1,j,s] += u[i,j]
            Σ[i,mod(j-1-1,Ny)+1,  s] += u[i,j]
            Σ[i,mod(j+1-1,Ny)+1,  s] += u[i,j]
            Σ[i,j,r] = 0
        end
    end
    function constant_step!(u,v,Σ,a,b,Nx,Ny,s,r,D,bk_integ_dt,h)
        @inbounds for i=1:Nx, j=1:Ny
            uth = (v[i,j] + b)/_get_a(a, i, j)
            v[i,j] = v[i,j] + bk_integ_dt*(u[i,j]^3 - v[i,j])
            u[i,j] = F(u[i,j], uth) + D*bk_integ_dt/h^2 *Σ[i,j,r]
            Σ[i,j,s] -= 4*u[i,j]
            i > 1  && (Σ[i-1,j,s] += u[i,j])
            i < Nx && (Σ[i+1,j,s] += u[i,j])
            j > 1  && (Σ[i,j-1,s] += u[i,j])
            j < Ny && (Σ[i,j+1,s] += u[i,j])
            Σ[i,j,r] = 0
        end
    end

    # evolve for skip time
    for j in 1:S÷bk_integ_dt
        periodic ? periodic_step!(u,v,Σ,a,b,Nx,Ny,s,r,D,bk_integ_dt,h) :
                   constant_step!(u,v,Σ,a,b,Nx,Ny,s,r,D,bk_integ_dt,h)
        r,s = s,r
    end
    push!(U, copy(u)); push!(V, copy(v)) # actual initial condition
    for j in 1:total
        for i in 1:every
            periodic ? periodic_step!(u,v,Σ,a,b,Nx,Ny,s,r,D,bk_integ_dt,h) :
                       constant_step!(u,v,Σ,a,b,Nx,Ny,s,r,D,bk_integ_dt,h)
            r,s = s,r
        end
        push!(U, copy(u)); push!(V, copy(v))
    end
    return 0:Δt:T, U, V
end

_get_a(a::Number, i, j) = a
_get_a(a::AbstractMatrix, i, j) = @inbounds a[i,j]
const BK_INTEG_DT = 0.01
