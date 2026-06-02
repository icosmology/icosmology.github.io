# lecture16_camb_workflow.jl
#
# A compact Julia example for the workflow discussed in Lecture 16.
# This uses PythonCall.jl to access the Python CAMB package, because CAMB's
# public frontend is Python/Fortran rather than native Julia.
#
# Requirements:
#   - Julia package PythonCall.jl
#   - Python package camb installed in the Python environment used by PythonCall
#
# Example setup in Julia:
#   using Pkg
#   Pkg.add("PythonCall")
#
# Then make sure CAMB is installed in the corresponding Python environment.

using PythonCall
using DelimitedFiles

camb = pyimport("camb")

# 1. Build the main CAMB parameter object
pars = camb.CAMBparams()

# 2. Set the background cosmology
pars.set_cosmology(H0=67.4, ombh2=0.0224, omch2=0.120, mnu=0.06, omk=0.0)

# 3. Set primordial power-spectrum parameters
pars.InitPower.set_params(As=2.1e-9, ns=0.965)

# 4. Set the CMB calculation range
pars.set_for_lmax(2500, lens_potential_accuracy=1)

# 5. Request transfer functions / matter power
pars.WantTransfer = true
pars.set_matter_power(redshifts=[0.0], kmax=2.0)

# 6. Run CAMB
results = camb.get_results(pars)

# 7. Get CMB power spectra
cls_py = results.get_cmb_power_spectra(pars, CMB_unit="muK")
total_cls = pyconvert(Array{Float64,2}, cls_py["total"])

# CAMB returns columns [TT, EE, BB, TE] in the "total" block.
ells = collect(0:size(total_cls, 1)-1)
tt = total_cls[:, 1]

# 8. Get the linear matter power spectrum at z = 0
kh_py, z_py, pk_py = results.get_matter_power_spectrum(
    minkh=1e-4,
    maxkh=1.0,
    npoints=200
)

kh = pyconvert(Vector{Float64}, kh_py)
zvals = pyconvert(Vector{Float64}, z_py)
pk = pyconvert(Array{Float64,2}, pk_py)

# With redshifts=[0.0], the first row corresponds to z = 0.
pk_z0 = vec(pk[1, :])

# 9. Save simple text outputs that students can inspect or plot later
writedlm("camb_tt_cls.txt", hcat(ells, tt))
writedlm("camb_pk_z0.txt", hcat(kh, pk_z0))

println("Saved TT spectrum to camb_tt_cls.txt")
println("Saved linear matter power spectrum at z = 0 to camb_pk_z0.txt")
println("Returned redshift grid = ", zvals)
