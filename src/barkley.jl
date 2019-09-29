using Random

barkley(T::Number, args...; kwargs...) =
barkley(Float32, T, args...; kwargs...)

function barkley(::Type{Ftype}, T, Δt, periodic=true;
                tskip=100, Nx = 50, Ny = Nx,
                a=0.75, b=0.06, ε=0.08, D=1/50, h=0.1,
                seed = rand(1:1000), kwargs...
                ) where {Ftype<:AbstractFloat}


    @assert Δt/bk_integ_dt == round(Int, Δt/bk_integ_dt) "Δt must be multiple of the integration step (0.01)"
    every = round(Int, Δt/bk_integ_dt)
    @assert Nx ≥ 40
    @assert Ny ≥ 40
    @assert Nx%10==0
    @assert Ny%10==0
    U = Vector{Array{Ftype,2}}()
    V = Vector{Array{Ftype,2}}()
    u = zeros(Ftype,Nx,Ny)
    v = zeros(Ftype,Nx,Ny)

    Random.seed!(seed)
    init = rand(10,10)
    Random.seed!()
    αα, ββ = size(init)

    v .= u .= repeat(init, inner=(Nx÷αα,Ny÷ββ))
    v .= v .> 0.2
    Σ = zeros(Ftype, Nx, Ny, 2)
    r,s = 1,2

    function F(u, uth)
        if u < uth
            u/(1-(bk_integ_dt/ε)*(1-u)*(u-uth))
        else
            (u + (bk_integ_dt/ε)*u*(u-uth))/(1+(bk_integ_dt/ε)*u*(u-uth))
        end
    end

struct Barkley{T}
    periodic::Bool
    a::T
    b::T
    D::T
    h::T
    Nx::Int
    Ny::Int
    s::Int8
    r::Int8
    Σ::Array{T, 3}
end
 # this is so fucking complicated, what the hell is Σ

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

    for m=1:(T*every+tskip*every)
        periodic ? periodic_step!(u,v,Σ,a,b,Nx,Ny,s,r,D,bk_integ_dt,h) :
                   constant_step!(u,v,Σ,a,b,Nx,Ny,s,r,D,bk_integ_dt,h)
        r,s = s,r
        if m > tskip*every && (m-tskip*every) % every == 0
            push!(U,copy(u))
            push!(V,copy(v))
        end
    end
    return U,V
end

_get_a(a::Number, i, j) = a
_get_a(a::AbstractMatrix, i, j) = @inbounds a[i,j]
