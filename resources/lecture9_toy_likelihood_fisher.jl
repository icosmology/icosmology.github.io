using LinearAlgebra
using Plots

chi2(x, y) = x^2 + ((y - 0.45 * x^2) / 0.6)^2
nll(x, y) = 0.5 * chi2(x, y)

x0 = 0.0
y0 = 0.0
h = 1.0e-3

# Numerical Hessian of chi^2 / 2 at the best-fit point.
F11 = (nll(x0 + h, y0) - 2nll(x0, y0) + nll(x0 - h, y0)) / h^2
F22 = (nll(x0, y0 + h) - 2nll(x0, y0) + nll(x0, y0 - h)) / h^2
F12 = (nll(x0 + h, y0 + h) - nll(x0 + h, y0 - h)
     - nll(x0 - h, y0 + h) + nll(x0 - h, y0 - h)) / (4h^2)

F = [F11 F12; F12 F22]
C = inv(F)

function ellipse_points(C, deltachi2; n = 400)
    eig = eigen(Symmetric(C))
    order = sortperm(eig.values; rev = true)
    vals = eig.values[order]
    vecs = eig.vectors[:, order]

    t = range(0, 2pi; length = n)
    circle = [cos.(t)'; sin.(t)']
    axes = Diagonal(sqrt.(vals .* deltachi2))
    pts = vecs * axes * circle
    return pts[1, :], pts[2, :]
end

xs = range(-2.0, 2.0; length = 401)
ys = range(-1.2, 2.8; length = 401)
chi = [chi2(x, y) for y in ys, x in xs]

x68, y68 = ellipse_points(C, 2.30)
x95, y95 = ellipse_points(C, 6.18)

plt = contour(xs, ys, chi;
    levels = [2.30, 6.18],
    aspect_ratio = :equal,
    xlabel = "x",
    ylabel = "y",
    linewidth = 2,
    label = false)
plot!(plt, x68, y68; linewidth = 2, linestyle = :dash, label = "Fisher 68%")
plot!(plt, x95, y95; linewidth = 2, linestyle = :dashdot, label = "Fisher 95%")
scatter!(plt, [x0], [y0]; markersize = 4, label = "best fit")

savefig(plt, "lecture9_toy_likelihood_fisher.pdf")
