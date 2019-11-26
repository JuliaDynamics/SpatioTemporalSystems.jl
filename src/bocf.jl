include("fast_convolutions.jl")

function makesim(sts::Union{STS{:BuenoOrovioCherryFenton}, STS{:bocf}})
    @unpack N, T, Δt, S, ic = sts
    Ttr = S
    Nx, Ny = typeof(N) <: Tuple ? (N[1], N[2]) : (N, N)

    # TODO: Change the model to use Float32 numbers instead

    u, v, w, s = [zeros(Nx, Ny) for i in 1:4]

    pset = get(sts.p, :pset, "tnpp")

    # here we get the initial condition
    if ic isa Union{Integer, Nothing} # random seed
        initialize_random!(u, ic isa Integer ? ic : 42)
    elseif ic == "spiral_virtheart"
        initialize_spiral_virtheart!(u, v)
    elseif ic == "spiral"
        initialize_spiral!(u, v, w)
    else error("unknown i.c.")
    end

    t, U, V, W, S = bocf(u, v, s, w, Ttr, T, Δt; sts.p...)
    @strdict t sts U V W S
end

################################################
# PDE Stepping
################################################
"""
    bocf(u, v, s, w, Ttr, T, Δt; kwargs...) → t, U, V, W, S
Simulate the Bueno-Orovio, Cherry & Fenton model for cardiac dynamics
with initial conditions `u, v, s, w`. Simulate for total time `T`
with sampling time `Δt`, while transiently integrating for `Ttr`.

A note on the initial condition, `sts.ic`.
The initial condition is by default `ic = 42`.
If number, random initial condition with `ic` seed is generated.
Otherwise `ic` can be `"spiral_virtheart"` or `"spiral"`, which initializes
a specific pattern.

The keywords of the BOCF model are (given to `sts.p` as a `NamedTuple`)
* `D = 0.2` the diffusion constant
* `pset = "tnpp"` the parameter set. Can be `"tnpp", "thomas", "epi", "virtheart", "pb"`.
  See the source code for the actual parameters (too many to list).
* `dt = Δt` : integration timestep (NOT the sampling time `Δt`).
"""
function bocf(u, v, s, w, Ttr, T, Δt;
        pset = "tnpp", D = 0.2, dt = Δt)
        # TODO: what is a good default dt?

    @assert Δt ≥ dt
    simpars = bocf_parameters(pset)
    pushfirst!(simpars, D) # add diffusion constant

    Nx, Ny = size(u)
    olds = [copy(u), copy(v), copy(w), copy(s)]
    derivs = [copy(u), copy(v), copy(w), copy(s)]
    aux = [zeros(Nx, Ny) for i in 1:7]
    currents = [zeros(Nx, Ny) for i in 1:3]

    TU = typeof(u)
    U, V, W, S = TU[], TU[], TU[], TU[]

    for j in 1:Ttr÷dt
        bocf_step!(u, v, w, s, olds, aux, currents, derivs, simpars, dt)
    end
    push!(U, copy(w)); push!(V, copy(v)); push!(W, copy(w)); push!(S, copy(s))
    for j in 1:Int(T÷Δt) # how many data to save
        for i in 1:Δt÷dt # how often to sample
            bocf_step!(u, v, w, s, olds, aux, currents, derivs, simpars, dt)
        end
        push!(U, copy(w)); push!(V, copy(v)); push!(W, copy(w)); push!(S, copy(s))
    end
    return 0:Δt:T, U, V, W, S

end

function bocf_step!(u, v, w, s, olds, aux, currents, derivs, simpars, dt)

    D,u_o,u_u,theta_v,theta_w,theta_v_minus,theta_o,tau_v1_minus,tau_v2_minus,
    tau_v_plus,tau_w1_minus,tau_w2_minus,k_w_minus,u_w_minus,tau_w_plus,
    tau_fi,tau_o1,tau_o2,tau_so1,tau_so2,k_so,u_so,tau_s1,tau_s2,k_s,
    u_s,tau_si,tau_winfinity,w_inf_star = simpars

    # update old fields (only necessary for boundary conditions)
    # TODO: this can be optimized to only update the values
    # used in the _set_boundaries! function
    for (new, old) in zip((u, v, w, s), olds)
        old .= new
    end

    _update_constants!(aux, u, v, w, s, simpars)
    tau_v_minus, tau_w_minus, tau_so, tau_s, tau_o, v_infinity, w_infinity = aux
    J_fi, J_so, J_si = currents

    # Compute currents
    _current_fi!(J_fi, u, v, theta_v, u_u, tau_fi)
    _current_so!(J_so, u, u_o, theta_w, tau_o, tau_so)
    _current_si!(J_si, u, theta_w, w, s, tau_si)

    dudt, dvdt, dwdt, dsdt = derivs
    lap = laplace(u)
    @. dudt = D * lap - (J_fi + J_so + J_si)
    @. dvdt = (1.0 - H(u - theta_v)) * (v_infinity - v) / tau_v_minus - H(u - theta_v) * v / tau_v_plus
    @. dwdt = (1.0 - H(u - theta_w)) * (w_infinity - w) / tau_w_minus - H(u - theta_w) * w / tau_w_plus
    @. dsdt = ((1.0 + tanh(k_s * (u - u_s))) / 2.0 - s) / tau_s

    @. u += dt * dudt
    @. v += dt * dvdt
    @. w += dt * dwdt
    @. s += dt * dsdt

    _set_boundaries!((u, v, w, s), olds)
    return u, v, w, s
end


function _set_boundaries!(fields, old_fields)
    # TODO: We have to check if this is correct, because I took this from
    # Python code, which interchanges columns with rows
    # TODO: What is happening in this function seems really weird for me.
    # Is it really correct? The end values before the step go the the near
    # end values after the step? That's so weird.
    for (field, old_field) in zip(fields, old_fields)
        field[:, 1] .= @view old_field[:, 2]
        field[:, end] .= @view old_field[:, end-1]
        field[1, :] .= @view old_field[2, :]
        field[end, :] .= @view old_field[end-1, :]
    end
end

# TODO: `laplace` creates arrays. It could be optmized by replacing the
# convolution with an explicit loop, because
# laplace = -4*copy(field)
# laplace += roll(field, +1, axis=0)
# laplace += roll(field, -1, axis=0)
# laplace += roll(field, +1, axis=1)
# laplace += roll(field, -1, axis=1)
# notice however that our fastconv is using an optimized FFT algorithm.
# I doubt that normal loops will be faster.
# Laplacian matrix
const Λ = [1  4  1;
           4 -20 4;
           1  4  1]./6
laplace(u) = fastconv(u, Λ)[2:end-1, 2:end-1]

H(y) = Float64(y > 0)

# TODO: All currents can be hugely optimized by normal loop over
# the operations, and if statements on the result of H(), so that
# all other operations are not done

function _current_fi!(cfi, u, v, theta_v, u_u, tau_fi)
    @. cfi = -v * H(u - theta_v) * (u - theta_v) * (u_u - u) / tau_fi
end

function _current_so!(cso, u, u_o, theta_w, tau_o, tau_so)
    @. cso = (u - u_o) * (1 - H(u - theta_w)) / tau_o + H(u - theta_w) / tau_so
end

function _current_si!(csi, u, theta_w, w, s, tau_si)
    @. csi = -H(u - theta_w) * w * s / tau_si
end

function _update_constants!(aux, u, v, w, s, simpars)
    tau_v_minus, tau_w_minus, tau_so, tau_s, tau_o, v_infinity, w_infinity = aux

    D,u_o,u_u,theta_v,theta_w,theta_v_minus,theta_o,tau_v1_minus,tau_v2_minus,
    tau_v_plus,tau_w1_minus,tau_w2_minus,k_w_minus,u_w_minus,tau_w_plus,
    tau_fi,tau_o1,tau_o2,tau_so1,tau_so2,k_so,u_so,tau_s1,tau_s2,k_s,
    u_s,tau_si,tau_winfinity,w_inf_star = simpars

    @. tau_v_minus = (1.0 - H(u - theta_v_minus)) * tau_v1_minus + H(u - theta_v_minus) * tau_v2_minus
    @. tau_w_minus = tau_w1_minus + (tau_w2_minus - tau_w1_minus) * (1 + tanh(k_w_minus * (u - u_w_minus))) / 2.0
    @. tau_so      = tau_so1      + (tau_so2      - tau_so1)      * (1 + tanh(k_so      * (u - u_so)))      / 2.0
    @. tau_s = (1.0 - H(u - theta_w)) * tau_s1 + H(u - theta_w) * tau_s2
    @. tau_o = (1.0 - H(u - theta_o)) * tau_o1 + H(u - theta_o) * tau_o2

    @. v_infinity = Float64(u .< theta_v_minus)
    @. w_infinity = (1.0 - H(u - theta_o)) * (1.0 - u / tau_winfinity) + H(u - theta_o) * w_inf_star
    return
end

################################################
# Initial conditions
################################################
function initialize_random!(u, ic)
    Nx, Ny = size(u)
    rng = MersenneTwister(ic)
    n = 1 # ceil(Int, 1/deltaX), what is deltaX ?
    u .= rand(rng, Nx, Ny)
end

function initialize_spiral_virtheart!(u, v)
    Nx, Ny = size(u)
    for i in 1:Ny
        for j in (Nx÷2):Nx
            t  = (i - Ny÷2)*0.1
            t2 = (i - Ny÷2 + 20)*0.05
            u[j, i] = 1.5*exp(-t*t)
            v[j, i] = 1.0 - 0.9*exp(-t2*t2)
        end
    end
end

function initialize_spiral!(u, v, w)
    Nx, Ny = size(u)
    for x in 1:(Nx÷4)
        for y in 1:Ny
            u[x, y] = 0.5
        end
    end
    for x in (3Nx÷4):Nx
        for y in 1:Ny
            u[x, y] = 0.5
        end
    end
    for x in 1:(Nx÷2)
        for y in 1:(2Ny÷5)
            v[x, y] = 0.5
        end
    end
    for x in (Nx÷2):Nx
        for y in (3Ny÷5):Ny
            v[x, y] = 0.5
        end
    end
    fill!(w, 1.0)
end

################################################
# Parameter sets
################################################
function bocf_parameters(ic)
    if ic == "tnpp"
        u_o = 0.0
        u_u = 1.58
        theta_v = 0.3
        theta_w = 0.015
        theta_v_minus = 0.015
        theta_o = 0.006
        tau_v1_minus = 60
        tau_v2_minus = 1150
        tau_v_plus = 1.4506
        tau_w1_minus = 70
        tau_w2_minus = 20
        k_w_minus = 65
        u_w_minus = 0.03
        tau_w_plus = 280
        tau_fi = 0.11
        tau_o1 = 6
        tau_o2 = 6
        tau_so1 = 43
        tau_so2 = 0.2
        k_so = 2
        u_so = 0.65
        tau_s1 = 2.7342
        tau_s2 = 3
        k_s = 2.0994
        u_s = 0.9087
        tau_si = 2.8723
        tau_winfinity = 0.07
        w_inf_star = 0.94
    elseif ic == "virtheart"
        u_o = 0
        u_u = 1.58
        theta_v = 0.3
        theta_w = 0.015
        theta_v_minus = 0.015
        theta_o = 0.006
        tau_v1_minus = 60
        tau_v2_minus = 60
        tau_v_plus = 1.4506
        tau_w1_minus = 170
        tau_w2_minus = 120
        k_w_minus = 65
        u_w_minus = 0.03
        tau_w_plus = 280
        tau_fi = 0.2
        tau_o1 = 6
        tau_o2 = 6
        tau_so1 = 43
        tau_so2 = 0.2
        k_so = 2
        u_so = 0.65
        tau_s1 = 2.7342
        tau_s2 = 3
        k_s = 2.0994
        u_s = 0.9087
        tau_si = 3.8723
        tau_winfinity = 0.07
        w_inf_star = 0.94
    elseif ic == "thomas"
        u_o = 0#13.03
        tau_o1 = 33.25
        theta_w = 800
        tau_so2 = 0.85
        tau_v1_minus = 0.45
        tau_s1 = 0.04
        tau_w1_minus = 0.45
        u_s = 0.04
        u_w_minus = 0.45
        w_inf_star = 0.04
        tau_fi = 12.5
        theta_v = 1250
        tau_so1 = 0.13
        theta_o = 0.45
        u_so = 0.04
        tau_v_plus = 0.45
        k_s = 0.04
        k_w_minus = 0.45
        tau_winfinity = 0.04
        u_u = 19.6
        tau_o2 = 29
        theta_v_minus = 40
        k_so = 0.04
        tau_v2_minus = 0.45
        tau_s2 = 0.04
        tau_w2_minus = 0.45
        tau_si = 0.04
        tau_w_plus = 0.45
    elseif ic == "epi"
        u_o = 0
        u_u = 1.55
        theta_v = 0.3
        theta_w = 0.13
        theta_v_minus = 0.006
        theta_o = 0.006
        tau_v1_minus = 60
        tau_v2_minus = 1150
        tau_v_plus = 1.4506
        tau_w1_minus = 60
        tau_w2_minus = 15
        k_w_minus = 65
        u_w_minus = 0.03
        tau_w_plus = 200
        tau_fi = 0.11
        tau_o1 = 400
        tau_o2 = 6
        tau_so1 = 30.0181
        tau_so2 = 0.9957
        k_so = 2.0458
        u_so = 0.65
        tau_s1 = 2.7342
        tau_s2 = 16
        k_s = 2.0994
        u_s = 0.9087
        tau_si = 1.8875
        tau_winfinity = 0.07
        w_inf_star = 0.94
    elseif ic == "pb"
        u_o = 0
        u_u = 1.45
        theta_v = 0.35
        theta_w = 0.13
        theta_v_minus = 0.175
        theta_o = 0.006
        tau_v1_minus = 10
        tau_v2_minus = 1150
        tau_v_plus = 1.4506
        tau_w1_minus = 140
        tau_w2_minus = 6.25
        k_w_minus = 65
        u_w_minus = 0.015
        tau_w_plus = 326
        tau_fi = 0.105
        tau_o1 = 400
        tau_o2 = 6
        tau_so1 = 30.0181
        tau_so2 = 0.9957
        k_so = 2.0458
        u_so = 0.65
        tau_s1 = 2.7342
        tau_s2 = 16
        k_s = 2.0994
        u_s = 0.9087
        tau_si = 1.8875
        tau_winfinity = 0.175
        w_inf_star = 0.9
    else
        error("Unknown parameter set")
    end

    simpars =
    [u_o,u_u,theta_v,theta_w,theta_v_minus,theta_o,tau_v1_minus,tau_v2_minus,
    tau_v_plus,tau_w1_minus,tau_w2_minus,k_w_minus,u_w_minus,tau_w_plus,
    tau_fi,tau_o1,tau_o2,tau_so1,tau_so2,k_so,u_so,tau_s1,tau_s2,k_s,
    u_s,tau_si,tau_winfinity,w_inf_star]
 end
