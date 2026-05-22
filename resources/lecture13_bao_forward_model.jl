using LinearAlgebra

"""
Pedagogical Julia script for the Lecture 13 BAO forward model.
It implements the logic
    theory multipoles on fine grid -> stack -> window convolution ->
    broadband terms -> chi^2
with a compact toy BAO model.
"""

struct BAOParams
    alpha_perp::Float64
    alpha_par::Float64
    B::Float64
    beta::Float64
    Sigma_par::Float64
    Sigma_perp::Float64
    Sigma_s::Float64
end

struct BroadbandParams
    a01::Float64
    a02::Float64
    a03::Float64
    a04::Float64
    a05::Float64
    a21::Float64
    a22::Float64
    a23::Float64
    a24::Float64
    a25::Float64
end

legendre0(mu) = 1.0
legendre2(mu) = 0.5 * (3.0 * mu^2 - 1.0)
legendre4(mu) = (35.0 * mu^4 - 30.0 * mu^2 + 3.0) / 8.0

function trapz(y::AbstractVector, x::AbstractVector)
    s = 0.0
    @inbounds for i in 1:(length(x) - 1)
        s += 0.5 * (y[i] + y[i + 1]) * (x[i + 1] - x[i])
    end
    return s
end

function kprime_muprime(k, mu, p::BAOParams)
    alpha = p.alpha_perp^(2 / 3) * p.alpha_par^(1 / 3)
    one_plus_eps = (p.alpha_par / p.alpha_perp)^(1 / 3)
    fac = sqrt(1.0 + mu^2 * (one_plus_eps^(-6) - 1.0))
    kp = k * one_plus_eps / alpha * fac
    mup = mu / one_plus_eps^3 / fac
    return kp, mup
end

# Toy no-wiggle spectrum and wiggles, just for demonstration.
Pnw_lin(k) = k * exp(-k / 0.25)
Owig(k) = 0.06 * sin(110.0 * k) * exp(-(k / 0.35)^2)

function Pg(k, mu, p::BAOParams)
    kp, mup = kprime_muprime(k, mu, p)
    fog = 1.0 / (1.0 + kp^2 * mup^2 * p.Sigma_s^2 / 2.0)
    Pnw = p.B^2 * (1.0 + p.beta * mup^2)^2 * Pnw_lin(kp) * fog
    damping = exp(-0.5 * kp^2 * (mup^2 * p.Sigma_par^2 + (1.0 - mup^2) * p.Sigma_perp^2))
    return Pnw * (1.0 + Owig(kp) * damping)
end

function multipoles(kfine::AbstractVector, p::BAOParams; nmu::Int = 401)
    mus = collect(range(-1.0, 1.0; length = nmu))
    P0 = similar(kfine, Float64)
    P2 = similar(kfine, Float64)
    P4 = similar(kfine, Float64)

    @inbounds for (ik, k) in pairs(kfine)
        vals = [Pg(k, mu, p) for mu in mus]
        P0[ik] = 0.5 * trapz([vals[i] * legendre0(mus[i]) for i in eachindex(mus)], mus)
        P2[ik] = 2.5 * trapz([vals[i] * legendre2(mus[i]) for i in eachindex(mus)], mus)
        P4[ik] = 4.5 * trapz([vals[i] * legendre4(mus[i]) for i in eachindex(mus)], mus)
    end

    return P0, P2, P4
end

function broadband_vector(kdata::AbstractVector, b::BroadbandParams)
    bb0 = b.a01 ./ kdata.^3 .+ b.a02 ./ kdata.^2 .+ b.a03 ./ kdata .+ b.a04 .+ b.a05 .* kdata
    bb2 = b.a21 ./ kdata.^3 .+ b.a22 ./ kdata.^2 .+ b.a23 ./ kdata .+ b.a24 .+ b.a25 .* kdata
    return vcat(bb0, bb2)
end

function evaluate_chi2(data::AbstractVector,
                       Cinv::AbstractMatrix,
                       W::AbstractMatrix,
                       kfine::AbstractVector,
                       kdata::AbstractVector,
                       p::BAOParams,
                       b::BroadbandParams)
    P0, P2, P4 = multipoles(kfine, p)
    t = vcat(P0, P2, P4)
    model = W * t + broadband_vector(kdata, b)
    r = data - model
    return dot(r, Cinv * r), model, t
end

function toy_window_matrix(Nd::Int, Nt::Int)
    W = zeros(2Nd, 3Nt)
    fine_per_data = max(1, cld(Nt, Nd))

    # Main monopole and quadrupole response.
    for i in 1:Nd
        j0 = min(Nt, max(1, (i - 1) * fine_per_data + 1))
        for dj in 0:2
            j = min(Nt, j0 + dj)
            W[i, j] += 0.25
            W[Nd + i, Nt + j] += 0.25
        end
    end

    # A small amount of multipole leakage from P4 into observed P0 and P2.
    for i in 1:Nd
        j = min(Nt, max(1, Int(round(i * Nt / Nd))))
        W[i, 2Nt + j] += 0.03
        W[Nd + i, 2Nt + j] += 0.06
    end

    return W
end

function main()
    Nd = 24
    Nt = 120
    kdata = collect(range(0.03, 0.27; length = Nd))
    kfine = collect(range(0.001, 0.40; length = Nt))

    p = BAOParams(1.01, 1.03, 1.6, 0.4, 6.0, 3.0, 3.0)
    b = BroadbandParams(0.0, 0.0, 0.0, 0.0, 0.0,
                        0.0, 0.0, 0.0, 0.0, 0.0)

    W = toy_window_matrix(Nd, Nt)
    data = zeros(2Nd)
    Cinv = Matrix{Float64}(I, 2Nd, 2Nd)

    chi2, model, theory = evaluate_chi2(data, Cinv, W, kfine, kdata, p, b)

    println("Toy Lecture 13 BAO forward model")
    println("length(data)   = ", length(data))
    println("size(W)        = ", size(W))
    println("length(theory) = ", length(theory))
    println("chi2           = ", chi2)
    println("first five model entries = ", model[1:5])
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
