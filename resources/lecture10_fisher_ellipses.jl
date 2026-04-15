using LinearAlgebra
using Plots

F = [8.0 -3.0;
     -3.0 5.0]

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

x68, y68 = ellipse_points(C, 2.30)
x95, y95 = ellipse_points(C, 6.18)

plt = plot(x68, y68;
    aspect_ratio = :equal,
    xlabel = "parameter 1",
    ylabel = "parameter 2",
    linewidth = 2,
    label = "68% ellipse")
plot!(plt, x95, y95; linewidth = 2, linestyle = :dash, label = "95% ellipse")
scatter!(plt, [0.0], [0.0]; markersize = 4, label = "fiducial point")

savefig(plt, "lecture10_fisher_ellipses.pdf")
