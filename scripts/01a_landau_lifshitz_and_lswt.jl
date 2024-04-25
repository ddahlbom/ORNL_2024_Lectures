# This script contains the simple examples shown on the first day of lectures.
# It demonstrates how to use Landau-Lifshitz (LL) dynamics in Sunny. It then
# gives a demonstration of the fact that, near T=0, LL reproduces linear SWT.

using DrWatson
@quickactivate "ORNL_2024_LECTURES"
using Sunny, LinearAlgebra, GLMakie

# Set up a "crystal" to support of spin system. Here we just want a linear spin
# chain, so we artificially set up a crystal that breaks symmetries along the b
# and c axes.
latvecs = lattice_vectors(1, 1.1, 2, 90, 90, 90)
positions = [[0, 0, 0]]
crystal = Crystal(latvecs, positions)
view_crystal(crystal)
print_symmetry_table(crystal, 1.1)


# Make a spin system based on this lattice.
dims = (10, 1, 1)
sys = System(crystal, dims, [SpinInfo(1; S=1, g=1)], :dipole_large_S; units=Sunny.Units.theory)

# Add ferromagnetic exchange.
J = -1.0
set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))

# Add an external field
h = 1.0
B = (0, 0, h)
set_external_field!(sys, B)

# Look at the system
plot_spins(sys)

# Go to random initial state
randomize_spins!(sys)
plot_spins(sys)

# Create a numerical integrator for the LL dynamics. We add damping so that the dynamics
# finds the ground state.
integrator = Langevin(; damping=0.1, kT=0.0)
suggest_timestep(sys, integrator; tol=1e-2)
integrator.dt = 0.05

# Now we'll create an animation, observing the system relax into the ground state.
fig = plot_spins(sys; colorfn=i->sys.dipoles[i][3], colorrange=(-1, 1))
for _ in 1:1000
    for _ in 1:5
        step!(sys, integrator)
    end
    notify(fig)
    sleep(1/60)
end

# For convenience, we define a function to create a spin chain.
function spin_chain(; J=1.0, D=0.0, h=0.0, dims=(10,1,1))
    latvecs = lattice_vectors(1, 1.1, 2, 90, 90, 90)
    positions = [[0, 0, 0]]
    crystal = Crystal(latvecs, positions)
    sys = System(crystal, dims, [SpinInfo(1; S=1, g=1)], :dipole_large_S; units=Sunny.Units.theory)
    set_exchange!(sys, J, Bond(1, 1, [1, 0, 0]))
    B = (0, 0, h)
    set_external_field!(sys, B)
    S = spin_matrices(Inf)
    set_onsite_coupling!(sys, D*S[3]^2, 1)

    return sys
end

# Let's make a system and find the ground state using the `minimize_energy!`
# function instead of the LL dynamics.
sys = spin_chain(; J=-1.0, D=-0.1, h=0.1, dims=(50, 1, 1))
minimize_energy!(sys)
fig = plot_spins(sys; colorfn=i->sys.dipoles[i][3], colorrange=(-1, 1))

# Let's thermalize slightly, i.e., perturb from the fully polarized ground
# state. To do this, we set the temperature field of the integrator and then run
# the dynamics.
kT = 0.05
integrator.kT = 0.05
for _ in 1:1000
    step!(sys, integrator)
    notify(fig)
    sleep(1/120)
end

# Then let's run a dissipationaless trajectory. This is the sort of trajectory
# that will be generated when calculating spin-spin correlations.
dt = 0.05
midpoint = ImplicitMidpoint(dt)
for _ in 1:1000
    step!(sys, midpoint)
    notify(fig)
    sleep(1/120)
end


# Use Sunny's tools to calculate the correlations of such trajectories.
nω = 200        # Number of energy bins
ωmax = 6.0      # Maximum resolved energy
sc = dynamical_correlations(sys; nω, ωmax, dt)  # Object to calculate and store correlations
add_sample!(sc, sys)

# Let's look at S(q,ω). We'll define a path in reciprocal space and tell Sunny how to calculate intensities.
points = [
    (-1, 0, 0),
    (-0.5, 0, 0),
    (0, 0, 0),
    (0.5, 0, 0),
    (1, 0, 0),
]
qs, xticks = reciprocal_space_path(sys.crystal, points, 100)
formula = intensity_formula(sc, :trace; kT)
is_classical = intensities_interpolated(sc, qs, formula)
heatmap(1:size(is_classical, 1), available_energies(sc), is_classical; colorrange=(0.0, 0.1), axis=(xticks=xticks,))

# Note that this is not very "clean." We only have one sample. We'll add 100 more.
for _ in 1:100
    for _ in 1:1000
        step!(sys, integrator)
    end
    add_sample!(sc, sys)
end

# Now plot the results.
is_classical = intensities_interpolated(sc, qs, formula)
heatmap(1:size(is_classical, 1), available_energies(sc), is_classical; colorrange=(0.0, 0.1), axis=(xticks=xticks,))


# Now let's compare to spin wave theory, which we have argued will produce essentially
# identical results.
sys_swt = spin_chain(; J=-1.0, D=-0.1, h=0.1, dims=(1, 1, 1))  # Magnetic unit cell contains a single spin.
swt = SpinWaveTheory(sys_swt)
formula_swt = intensity_formula(swt, :trace; kernel=lorentzian(0.05))
energies = range(0.0, ωmax, nω)
ΔE = energies[2]-energies[1]
is_swt = intensities_broadened(swt, qs, energies, formula_swt)

begin
    fig = Figure()
    ax_classical = Axis(fig[1,1]; xticks, title="Classical")
    ax_swt = Axis(fig[1,2]; xticks, title="LSWT")
    colorrange = (0.0, 0.1)

    heatmap!(ax_classical, 1:size(is_classical, 1), available_energies(sc), is_classical; colorrange)
    heatmap!(ax_swt, 1:size(is_swt, 1), energies, is_swt .* ΔE; colorrange)

    fig
end


# An advantage of the LL approach is that we can immediately proceed to examine
# finite temperature results, which take into account all the nonlinearities of
# the LL equations.

# We'll create a new correlations object and add 100 samples.
sc_hiT = dynamical_correlations(sys; nω, ωmax, dt)
integrator.kT = 1.0  # Set the temperature to 1 J.
for _ in 1:100
    for _ in 1:1000
        step!(sys, integrator)
    end
    add_sample!(sc_hiT, sys)
end

is_hiT = intensities_interpolated(sc_hiT, qs, formula)
heatmap(1:size(is_classical, 1), available_energies(sc), is_hiT; colorrange=(0.0, 0.25), axis=(xticks=xticks, title="High T"))

# Note that the broadening effect observed here is not artificial but an
# intrinsic product of the full nonlinear dynamics.