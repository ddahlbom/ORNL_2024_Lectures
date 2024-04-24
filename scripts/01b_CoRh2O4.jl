# This script shows the basic steps necessary to calculate a finite temperature
# powder averaged spectrum for the material CoRh₂O₄. By manipulating the
# temperature, it is possible to reproduce figures similar to the experimental
# results reported in L. Ge et al, "Spin order and dynamics in the
# diamond-lattice Heisenberg antiferromagnets...", PRB 98 (2018):
# https://doi.org/10.1103/PhysRevB.96.064413

using DrWatson
@quickactivate "ORNL_2024_LECTURES"
using Sunny, GLMakie

# Define a diamond-lattice
a = 8.5031 # (Å)
latvecs = lattice_vectors(a, a, a, 90, 90, 90)
cryst = Crystal(latvecs, [[0,0,0]], 227, setting="1")
view_crystal(cryst)

# Make a system and add NN AFM exchange.
latsize = (1, 1, 1)
S = 3/2
J = 7.5413*meV_per_K # (~ 0.65 meV)
sys = System(cryst, latsize, [SpinInfo(1; S, g=2)], :dipole; seed=0)
set_exchange!(sys, J, Bond(1, 3, [0,0,0]))

# For computational efficiency, reshape into the primitive cell.
shape = [0 1 1;
         1 0 1;
         1 1 0] / 2
sys_prim = reshape_supercell(sys, shape)
randomize_spins!(sys_prim)
minimize_energy!(sys_prim)
plot_spins(sys_prim)

# Now make a larger supercell.
sys_large = repeat_periodically(sys_prim, (14, 14, 4))
minimize_energy!(sys_large)
plot_spins(sys_large)


# Set up numerical integrator to generate finite temperature dynamics.
kT = 30.0 # K
integrator = Langevin(; damping=0.1, kT=kT*Sunny.meV_per_K)
suggest_timestep(sys, integrator; tol=1e-2)
integrator.dt = dt = 0.025


# Set up an object to contain spin-spin correlation data and begin accumulate samples.
# This may take about a minute.
nω = 200
ωmax = 6.0
sc = dynamical_correlations(sys_large; dt, nω, ωmax)
@time for _ in 1:10
    for _ in 1:1000
        step!(sys_large, integrator)
    end
    add_sample!(sc, sys_large)
end

# Now tell Sunny how to calculate the intensities at any wave vector. 
formfactors = [FormFactor("Co2")]
formula = intensity_formula(sc, :perp; formfactors)

# Perform a powder average.
radii = 0.01:0.02:3 # (1/Å)
output = zeros(Float64, length(radii), length(available_energies(sc)))
for (i, radius) in enumerate(radii)
    n = 300
    qs = reciprocal_space_shell(cryst, radius, n)
    is = intensities_interpolated(sc, qs, formula)
    output[i, :] = sum(is, dims=1) / size(is, 1)
end

# Plot the results
fig = Figure()
ax = Axis(fig[1,1]; xlabel="Q (Å⁻¹)", ylabel="ω (meV)")
heatmap!(ax, radii, energies, output, colormap=:gnuplot2, colorrange=(0.0, 0.1))
fig


# It is a useful exercise to try this at several different temperatures.