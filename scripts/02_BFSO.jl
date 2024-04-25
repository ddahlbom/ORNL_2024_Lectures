# This script demonstrates how to reproduce the magnetization curve reported for
# Ba₂FeSi₂O₇ in M. Lee et al., "Field-induced level crossings...", PRB 107
# (2023): https://doi.org/10.1103/PhysRevB.107.144427

using DrWatson
@quickactivate "ORNL_2024_LECTURES"

using Sunny, GLMakie, LinearAlgebra

# Function that creates a model of Ba₂FeSi₂O₇
function BFSO(dims; J=1.028, B = [0.0 0.0 0.0], g = 1.93)
    mpK = Sunny.meV_per_K

    # Set up crystal
    a = b =  8.3194; c =  5.3336
    α = β = γ = 90 
    latvecs = lattice_vectors(a, b, c, α, β, γ)
    types = ["Fe"]
    basisvecs = [[0.0, 0, 0]]
    cryst = Crystal(latvecs, basisvecs, 113; types)

    # Set up system
    sys = System(cryst, dims, [SpinInfo(1; S=2, g)], :SUN)

    # Nearest neighbor in-plane interactions
    J₁ = J * mpK
    set_exchange!(sys, J₁, Bond(1, 2, [0, 0, 0]))

    # Next-nearest neighbor in-plane interactions
    J₂ = 0.1*J₁
    set_exchange!(sys, J₂, Bond(1, 1, [1, 0, 0]))

    # Nearest out-of-plane interaction
    J′₁ = 0.1*J₁
    set_exchange!(sys, J′₁ , Bond(1, 1, [0, 0, 1]))

    # Single-ion anisotropy
    A   = 1.16 * mpK   
    C   = -1.74 * mpK  
    D   = 28.65 * mpK 
    Sˣ, Sʸ, Sᶻ = spin_matrices(2) 
    Λ = D*(Sᶻ)^2 + A*((Sˣ)^4 + (Sʸ)^4) + C*(Sᶻ)^4
    set_onsite_coupling!(sys, Λ, 1)

    # External field
    set_external_field!(sys, B) 

    return sys
end

# Initialize the system and perform a quick optimization to determine an
# approximate ground state. This should be staggered ordering in the plane with
# reduced dipole moment.
sys = BFSO((5, 5, 5); B=[0, 0, 0.0])
plot_spins(sys) 
randomize_spins!(sys)
minimize_energy!(sys; maxiters=10_000)
plot_spins(sys)

# Function to calculation magnetization
function magnetization(sys)
    M_avg = sum([magnetic_moment(sys, site) for site in Sunny.eachsite(sys)]) / prod(size(sys.dipoles)) / sys.units.μB
    return norm(M_avg)
end

# Function to calculate that average dipole moment on each site.
function avg_dipole(sys)
    sum(norm.(sys.dipoles)) / prod(size(sys.dipoles))
end



# Iterate through a range of applied field values (in Tesla). Minimize the
# minimize the energy and calculate the resulting magnetization and average
# dipole magnitude.
Bs = range(0.0, 55.0, 50)  # Applied field values.
Ms = Float64[]             # Vector to record magnetization at each field value.
dips = Float64[]           # Vector to record average dipole magnitude.
for B in Bs
    set_external_field!(sys, [0,0,B])
    minimize_energy!(sys; maxiters=1_000)

    M = magnetization(sys)
    push!(Ms, M)

    dip = avg_dipole(sys)
    push!(dips, dip)
end

# Plot the results.
fig = Figure(; resolution = (1200,500))
ax = Axis(fig[1,1]; ylabel = "M (μ_B/Fe²⁺)", xlabel = "μ₀B (T)", title="Magnetization")
scatter!(ax, Bs, Ms)
ax = Axis(fig[1,2]; ylabel = "Mean |S|", xlabel = "μ₀B (T)", title="Average magnitude of dipole moments")
scatter!(ax, Bs, dips)
fig
