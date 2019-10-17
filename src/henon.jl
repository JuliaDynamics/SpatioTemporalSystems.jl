using Random

function makesim(sts::Union{STS{:CoupleHenon1D}, STS{:henon1d}})
    @unpack N, T, S, ic = sts
    @assert T isa Integer
    @assert S isa Integer
    # prepare initial condition...
    if ic isa Union{Nothing, Integer}
        rng = MersenneTwister(isnothing(ic) ? 27 : ic)
        u = 0.1rand(rng, N)
        v = 0.1rand(rng, N)
    elseif ic isa Tuple
        u, v = ic
    else
        error(ICERROR)
    end
    t, U, V = coupled_henon_1D(u, v, S, T; sts.p...)
    return @strdict t sts U V
end

"""
    coupled_henon_1D(u, v, S, T; a = 1.4, b = 0.3)
Coupled henon maps in one dimension.
"""
function coupled_henon_1D(u, v, S, T; a = 1.4, b = 0.3)
    @assert length(u) == length(v)
    N = length(u)
    u = copy(u); v = copy(v)

    function henon(U,V)
        Un = copy(U)
        Vn = copy(V)
        Un[1] = Un[N] = 0.5
        Vn[1] = Vn[N] = 0
        for m=2:N-1
            @inbounds Un[m] = 1 - a*(.5*U[m]+ .25*U[m-1] + .25*U[m+1])^2 + b*V[m]
            @inbounds Vn[m] = U[m]
        end
        return Un, Vn
    end

    for n in 1:S
        u, v = henon(u, v)
    end
    U = Vector{Vector{Float32}}(undef,T)
    V = Vector{Vector{Float32}}(undef,T)
    U[1] = u; V[1] = v
    for n = 2:T
        U[n], V[n] = henon(U[n-1],V[n-1])
    end
    return 1:T, U,V
end
