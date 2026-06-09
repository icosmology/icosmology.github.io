# lecture18_desi_bao_hmc_example.jl
# ---------------------------------------------------------------
# Pedagogical Julia example for Lecture 18.
#
# Goal:
#   1. use public DESI DR2 BAO mean/covariance values,
#   2. fit the flat-LCDM parameter pair (Omega_m, H0*rd),
#   3. demonstrate automatic differentiation with ForwardDiff,
#   4. implement a compact Hamiltonian Monte Carlo sampler,
#   5. reproduce a DESI-style contour figure in the (Omega_m, H0*rd) plane.
#
# Notes:
#   - This is a teaching script, not a production inference pipeline.
#   - The combined posterior is sampled with HMC.
#   - Individual-tracer contours are plotted from direct chi^2 grids because
#     that is the clearest way to reproduce the DESI-style figure.
#   - Public data values are taken from the DESI DR2 BAO files distributed in
#     the Cobaya BAO repository.
#
# Suggested packages:
#   using Pkg
#   Pkg.add(["ForwardDiff", "Plots"])
#
# ---------------------------------------------------------------

using LinearAlgebra
using Random
using Statistics
using Printf
using ForwardDiff
using Plots

# ---------------------------------------------------------------
# Public DESI DR2 BAO data (combined mean vector and covariance)
# ---------------------------------------------------------------

const c_light = 299792.458  # km / s

const z_data = [
    0.295,
    0.510, 0.510,
    0.706, 0.706,
    0.934, 0.934,
    1.321, 1.321,
    1.484, 1.484,
    2.33,  2.33,
]

# kind[i] is one of :DV, :DM, :DH
const kind_data = [
    :DV,
    :DM, :DH,
    :DM, :DH,
    :DM, :DH,
    :DM, :DH,
    :DM, :DH,
    :DH, :DM,
]

const data_vector = [
    7.94167639,
    13.58758434, 21.86294686,
    17.35069094, 19.45534918,
    21.57563956, 17.64149464,
    27.60085612, 14.17602155,
    30.51190063, 12.81699964,
    8.631545674846294, 38.988973961958784,
]

const covariance = [
    0.00578998687 0 0 0 0 0 0 0 0 0 0 0 0;
    0 0.0283473742 -0.0326062007 0 0 0 0 0 0 0 0 0 0;
    0 -0.0326062007 0.18392804 0 0 0 0 0 0 0 0 0 0;
    0 0 0 0.0323752442 -0.0237445646 0 0 0 0 0 0 0 0;
    0 0 0 -0.0237445646 0.111469198 0 0 0 0 0 0 0 0;
    0 0 0 0 0 0.0261732816 -0.0112938006 0 0 0 0 0 0;
    0 0 0 0 0 -0.0112938006 0.0404183878 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0.105336516 -0.0290308418 0 0 0 0;
    0 0 0 0 0 0 0 -0.0290308418 0.0504233092 0 0 0 0;
    0 0 0 0 0 0 0 0 0 0.583020277 -0.195215562 0 0;
    0 0 0 0 0 0 0 0 0 -0.195215562 0.268336193 0 0;
    0 0 0 0 0 0 0 0 0 0 0 0.0102136194 -0.0231395216;
    0 0 0 0 0 0 0 0 0 0 0 -0.0231395216 0.282685779
]

const Cinv = inv(covariance)

# Tracer blocks used for the DESI-style figure.
# Indices are 1-based, as in Julia.
const blocks = [
    ("BGS",        [1]),
    ("LRG1",       [2, 3]),
    ("LRG2",       [4, 5]),
    ("LRG3+ELG1",  [6, 7]),
    ("ELG2",       [8, 9]),
    ("QSO",        [10, 11]),
    ("Lyα",        [12, 13]),
    ("All DESI DR2 BAO", collect(1:length(data_vector))),
]

const block_colors = Dict(
    "BGS" => :yellowgreen,
    "LRG1" => :orange,
    "LRG2" => :orangered,
    "LRG3+ELG1" => :slateblue,
    "ELG2" => :steelblue,
    "QSO" => :seagreen,
    "Lyα" => :purple,
    "All DESI DR2 BAO" => :black,
)

# ---------------------------------------------------------------
# Flat-LCDM BAO theory
# ---------------------------------------------------------------

sigmoid(x) = inv(one(x) + exp(-x))

E_of_z(z, Omega_m) = sqrt(Omega_m * (1 + z)^3 + (one(Omega_m) - Omega_m))

"""
    chiint(z, Omega_m; n=800)

Numerical approximation to ∫_0^z dz' / E(z') using Simpson's rule.
This implementation is friendly to ForwardDiff dual numbers because it
uses ordinary arithmetic only.
"""
function chiint(z, Omega_m; n = 800)
    @assert iseven(n) "Simpson rule needs an even number of subintervals."
    h = z / n
    s = inv(E_of_z(zero(z), Omega_m)) + inv(E_of_z(z, Omega_m))
    for i in 1:(n - 1)
        zi = i * h
        s += (isodd(i) ? 4 : 2) / E_of_z(zi, Omega_m)
    end
    return h * s / 3
end

function bao_model_vector(Omega_m, H0rd)
    A = c_light / H0rd
    T = typeof(A + Omega_m)
    μ = Vector{T}(undef, length(data_vector))
    for i in eachindex(data_vector)
        z = z_data[i]
        Ez = E_of_z(z, Omega_m)
        χ = chiint(z, Omega_m)
        if kind_data[i] == :DM
            μ[i] = A * χ
        elseif kind_data[i] == :DH
            μ[i] = A / Ez
        elseif kind_data[i] == :DV
            μ[i] = A * (z * χ^2 / Ez)^(one(T) / 3)
        else
            error("Unknown BAO observable type")
        end
    end
    return μ
end

function chi2(indices, Omega_m, H0rd)
    μ = bao_model_vector(Omega_m, H0rd)[indices]
    d = data_vector[indices]
    C = covariance[indices, indices]
    Δ = d .- μ
    return dot(Δ, inv(C) * Δ)
end

# ---------------------------------------------------------------
# Posterior parameterization for HMC
# ---------------------------------------------------------------
# We sample unconstrained u = (u1, u2) in R^2 and map to physical space via
# logistic transforms. This avoids hard boundaries in HMC.
#
#   Omega_m in [0.05, 0.60]
#   H0rd    in [9000, 11500] km/s
#
# A uniform prior in physical space corresponds to adding the log Jacobian.
# ---------------------------------------------------------------

function unpack(u)
    s1, s2 = sigmoid(u[1]), sigmoid(u[2])
    Omega_m = 0.05 + 0.55 * s1
    H0rd = 9000.0 + 2500.0 * s2
    logJ = log(0.55 * s1 * (1 - s1)) + log(2500.0 * s2 * (1 - s2))
    return Omega_m, H0rd, logJ
end

function logposterior_u(u)
    Omega_m, H0rd, logJ = unpack(u)
    μ = bao_model_vector(Omega_m, H0rd)
    Δ = data_vector .- μ
    return -0.5 * dot(Δ, Cinv * Δ) + logJ
end

U(u) = -logposterior_u(u)
gradU(u) = ForwardDiff.gradient(U, u)

# ---------------------------------------------------------------
# A compact pedagogical HMC implementation
# ---------------------------------------------------------------

function leapfrog(q, p, ϵ, L)
    qnew = copy(q)
    pnew = copy(p)

    pnew .-= 0.5 * ϵ * gradU(qnew)
    for i in 1:L
        qnew .+= ϵ * pnew
        if i != L
            pnew .-= ϵ * gradU(qnew)
        end
    end
    pnew .-= 0.5 * ϵ * gradU(qnew)
    return qnew, -pnew
end

function hmc(q0; nsamples = 6000, nwarmup = 2000, ϵ = 0.03, L = 25, rng = MersenneTwister(42))
    q = copy(q0)
    chain = Matrix{Float64}(undef, nsamples, length(q0))
    accepted = 0

    for n in 1:(nwarmup + nsamples)
        p0 = randn(rng, length(q0))
        current_H = U(q) + 0.5 * dot(p0, p0)

        qprop, pprop = leapfrog(q, p0, ϵ, L)
        prop_H = U(qprop) + 0.5 * dot(pprop, pprop)
        α = min(1.0, exp(current_H - prop_H))

        if rand(rng) < α
            q = qprop
            if n > nwarmup
                accepted += 1
            end
        end

        if n > nwarmup
            chain[n - nwarmup, :] .= q
        end
    end

    return chain, accepted / nsamples
end

function to_physical(chain_u)
    N = size(chain_u, 1)
    Ω = Vector{Float64}(undef, N)
    Hrd = Vector{Float64}(undef, N)
    for i in 1:N
        Ω[i], Hrd[i], _ = unpack(view(chain_u, i, :))
    end
    return Ω, Hrd
end

# ---------------------------------------------------------------
# Grid calculations for DESI-style contour plotting
# ---------------------------------------------------------------

function chi2_grid(indices, Ωgrid, Hgrid)
    χ2 = Matrix{Float64}(undef, length(Hgrid), length(Ωgrid))
    C = covariance[indices, indices]
    Cinv_block = inv(C)
    d = data_vector[indices]

    for (j, H0rd) in enumerate(Hgrid)
        for (i, Ωm) in enumerate(Ωgrid)
            μ = bao_model_vector(Ωm, H0rd)[indices]
            Δ = d .- μ
            χ2[j, i] = dot(Δ, Cinv_block * Δ)
        end
    end
    χ2 .-= minimum(χ2)
    return χ2
end

const Δχ2_levels = [2.30, 6.18]  # 68% and 95% for 2 parameters

function make_desi_style_plot(; Ωgrid = range(0.16, 0.45, length = 220),
                                Hgrid = range(9400.0, 10950.0, length = 220),
                                outfile = "lecture18_desi_bao_hmc_demo.png")
    plt = plot(
        xlabel = "Ω_m",
        ylabel = "H₀ r_d [km s⁻¹]",
        title = "DESI DR2 BAO-only constraints in flat ΛCDM",
        legend = :bottomleft,
        size = (1100, 800),
        linewidth = 2.0,
        guidefontsize = 18,
        tickfontsize = 13,
        titlefontsize = 22,
        legendfontsize = 12,
        xlims = (first(Ωgrid), last(Ωgrid)),
        ylims = (first(Hgrid), last(Hgrid)),
    )

    for (name, inds) in blocks
        χ2 = chi2_grid(inds, Ωgrid, Hgrid)
        ls = name == "All DESI DR2 BAO" ? :solid : :dash
        lw = name == "All DESI DR2 BAO" ? 3.0 : 1.8
        contour!(
            plt,
            Ωgrid,
            Hgrid,
            χ2,
            levels = Δχ2_levels,
            color = block_colors[name],
            linestyle = ls,
            linewidth = lw,
            label = name,
        )
    end

    savefig(plt, outfile)
    return plt
end

# ---------------------------------------------------------------
# Simple posterior summary
# ---------------------------------------------------------------

function summarize_samples(Ω, Hrd)
    μΩ = mean(Ω)
    μH = mean(Hrd)
    σΩ = std(Ω)
    σH = std(Hrd)
    corr = cor(Ω, Hrd)

    println("Posterior summary from HMC chain")
    println("--------------------------------")
    @printf("Omega_m  = %.6f ± %.6f\n", μΩ, σΩ)
    @printf("H0 * rd  = %.2f ± %.2f km/s\n", μH, σH)
    @printf("corr     = %.4f\n", corr)
end

# ---------------------------------------------------------------
# Main script
# ---------------------------------------------------------------

function main()
    println("Running pedagogical HMC fit for DESI DR2 BAO-only flat-LCDM example...")

    # A reasonable starting point near the combined contour.
    q0 = [0.0, 0.0]

    chain_u, accept_rate = hmc(q0; nsamples = 7000, nwarmup = 3000, ϵ = 0.028, L = 28)
    Ω, Hrd = to_physical(chain_u)

    println()
    @printf("Acceptance rate = %.3f\n\n", accept_rate)
    summarize_samples(Ω, Hrd)

    println()
    println("Making DESI-style contour plot from direct chi^2 grids...")
    make_desi_style_plot(outfile = "lecture18_desi_bao_hmc_demo.png")
    println("Saved figure to lecture18_desi_bao_hmc_demo.png")

    println()
    println("A few practical comments:")
    println("  * The combined posterior is sampled with HMC using ForwardDiff gradients.")
    println("  * The contour figure is drawn from chi^2 grids for clarity and robustness.")
    println("  * In a production analysis, one would normally rely on Cobaya/CosmoMC or")
    println("    another mature inference framework rather than a hand-written sampler.")
end

main()
