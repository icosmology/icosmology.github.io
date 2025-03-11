using Plots
using LaTeXStrings

Om0=0.3
Ok0=0.01
OR0=0.0001
OL0=1-Om0-Ok0-OR0
amin=1e-10
amax=1

i = range(1, 10000, length=10000)
a = amin*(amax/amin).^((i .-1)/9999.0)
z = 1 ./a .-1

Om= Om0*a.^(-3)
OL= OL0
OK= Ok0*a.^(-2)
OR= OR0*a.^(-4)

E=Om .+ OL .+ OK .+ OR

Oma = Om ./E
OLa = OL ./E
OKa = OK ./E
ORa = OR ./E

plot!(z, Oma, xscale=:log10, label=L"\Omega_{\rm M}(z)", framestyle=:box, 
legend=:bottomright, grid=:false, minorticks=:true)
plot!(z, OLa, xscale=:log10, label=L"\Omega_{\Lambda}(z)")
plot!(z, OKa, xscale=:log10, label=L"\Omega_{\rm K}(z)")
plot!(z, ORa, xscale=:log10, label=L"\Omega_{\rm R}(z)")

xlims!(0.01, 1e9)
ylims!(0, 1.1)
xlabel!(L"redshift $z$")
ylabel!("Relative fraction")

