# Pedagogical Julia implementation of weighted DD, DR, and RR pair counts
# for anisotropic two-point correlation-function measurements.
#
# This version is intentionally simple and uses O(N^2) / O(N_D N_R) loops.
# For real surveys one would replace the brute-force loops with a kd-tree,
# cell-pair search, or another accelerated neighbor-finding method.

using LinearAlgebra

struct Catalog
    x::Matrix{Float64}   # N x 3 Cartesian coordinates
    w::Vector{Float64}   # object weights
end

function bin_index(value::Float64, edges::AbstractVector{<:Real})
    i = searchsortedlast(edges, value)
    return (1 <= i < length(edges)) ? Int(i) : 0
end

function auto_norm(weights::AbstractVector{<:Real})
    s1 = sum(weights)
    s2 = sum(abs2, weights)
    return 0.5 * (s1^2 - s2)
end

cross_norm(w1::AbstractVector{<:Real}, w2::AbstractVector{<:Real}) = sum(w1) * sum(w2)

function paircounts_DD(cat::Catalog,
                       s_edges::AbstractVector{<:Real},
                       mu_edges::AbstractVector{<:Real})
    ns = length(s_edges) - 1
    nm = length(mu_edges) - 1
    hist = zeros(Float64, ns, nm)

    x = cat.x
    w = cat.w
    N = size(x, 1)
    normsum = auto_norm(w)

    @inbounds for i in 1:(N - 1)
        x1 = x[i, 1]
        y1 = x[i, 2]
        z1 = x[i, 3]
        wi = w[i]

        for j in (i + 1):N
            dx = x[j, 1] - x1
            dy = x[j, 2] - y1
            dz = x[j, 3] - z1
            s = sqrt(dx^2 + dy^2 + dz^2)
            s == 0.0 && continue

            nx = x[j, 1] + x1
            ny = x[j, 2] + y1
            nz = x[j, 3] + z1
            nmag = sqrt(nx^2 + ny^2 + nz^2)
            nmag == 0.0 && continue

            mu = abs((dx * nx + dy * ny + dz * nz) / (s * nmag))

            is = bin_index(s, s_edges)
            im = bin_index(mu, mu_edges)
            if is > 0 && im > 0
                hist[is, im] += wi * w[j]
            end
        end
    end

    return hist / normsum
end

function paircounts_DR(data::Catalog, randoms::Catalog,
                       s_edges::AbstractVector{<:Real},
                       mu_edges::AbstractVector{<:Real})
    ns = length(s_edges) - 1
    nm = length(mu_edges) - 1
    hist = zeros(Float64, ns, nm)

    xd = data.x
    wd = data.w
    xr = randoms.x
    wr = randoms.w

    Nd = size(xd, 1)
    Nr = size(xr, 1)
    normsum = cross_norm(wd, wr)

    @inbounds for i in 1:Nd
        x1 = xd[i, 1]
        y1 = xd[i, 2]
        z1 = xd[i, 3]
        wi = wd[i]

        for j in 1:Nr
            dx = xr[j, 1] - x1
            dy = xr[j, 2] - y1
            dz = xr[j, 3] - z1
            s = sqrt(dx^2 + dy^2 + dz^2)
            s == 0.0 && continue

            nx = xr[j, 1] + x1
            ny = xr[j, 2] + y1
            nz = xr[j, 3] + z1
            nmag = sqrt(nx^2 + ny^2 + nz^2)
            nmag == 0.0 && continue

            mu = abs((dx * nx + dy * ny + dz * nz) / (s * nmag))

            is = bin_index(s, s_edges)
            im = bin_index(mu, mu_edges)
            if is > 0 && im > 0
                hist[is, im] += wi * wr[j]
            end
        end
    end

    return hist / normsum
end

function paircounts_RR(randoms::Catalog,
                       s_edges::AbstractVector{<:Real},
                       mu_edges::AbstractVector{<:Real})
    return paircounts_DD(randoms, s_edges, mu_edges)
end

landy_szalay(DD, DR, RR) = (DD .- 2.0 .* DR .+ RR) ./ RR

# Example workflow:
#   data = Catalog(x_data, w_data)
#   rand = Catalog(x_rand, w_rand)
#   DD = paircounts_DD(data, s_edges, mu_edges)
#   DR = paircounts_DR(data, rand, s_edges, mu_edges)
#   RR = paircounts_RR(rand, s_edges, mu_edges)
#   xi = landy_szalay(DD, DR, RR)
