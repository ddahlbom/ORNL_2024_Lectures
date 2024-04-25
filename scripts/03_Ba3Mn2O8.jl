# This script shows how to reproduce the dispersion and intensities for the
# material Ba₃Mn₂O₈, as reported in M. Stone et al., "Singlet triplet dispersion
# reveals additional frustration...," PRL 100 (2008):
# https://doi.org/10.1103/PhysRevLett.100.237201
#
# Note that this script uses unreleased Sunny features. To execute this script
# (as of April 24, 2024), it is necessary to to work on the `entangled-units`
# development branch of Sunny.jl. A proper release with documentation and a
# tutorial will be made available in the coming months.


# Load the project environment. This will ensure you are using the correct Sunny
# branch.
using DrWatson
@quickactivate "ORNL_2024_LECTURES"

# Import the necessary packages and helper functions. See the `src` directory
# for `dispersion_relation.jl` and `model.jl`.
using Sunny, LinearAlgebra, GLMakie
includet(srcdir("dispersion_relation.jl"))
includet(srcdir("model.jl"))

# Prepare the system and get symmetry information.
sys = system_for_contraction()
view_crystal(sys.crystal)
print_symmetry_table(sys.crystal, 10.0)

# Reshape the system into the primitive cell 
(a1, a2, a3) = eachcol(sys.crystal.prim_latvecs)
shape = (sys.crystal.latvecs \ [-a1 a2 (-a3 + a1)])
sys_reshaped = reshape_supercell(sys, shape)
plot_spins(sys_reshaped)

# Make an entangled system
units = [(1,2)]
sys_entangled = EntangledSystem(sys_reshaped, units)

# Set the state to S=1 singlet (alternatively, use `randomize_spins!` followed by `minimize_energy!`)
gs_ref = Sunny.CVec{9}(
    0.0, 0.0, √3/3, 0.0, -√3/3, 0.0, √3/3, 0, 0
)
for site in Sunny.eachsite(sys_entangled)
    set_coherent!(sys_entangled, gs_ref, site)
end

# Confirm that we're in a singlet state (no expected dipole moment).
plot_spins(sys_entangled)

# Make an EntangledSpinWaveTheory
swt = EntangledSpinWaveTheory(sys_entangled)

# Calculate dispersion and intensities
begin
    # Define the broadening kernel
    FWHM = 0.295
    σ = FWHM/2.355
    kernel = ω -> 1/√(2π*σ^2)*exp(-ω^2/(2σ^2))

    # Define a path in reciprocal space and determine energies to examine.
    points_ref = [
        [0.175, 0.175, 1.5],
        [0.85, 0.85, 1.5],
        [0.85, 0.85, 3],
        [0.0, 0.0, 3],
        [0.0, 0.0, 8],
    ]
    path, xticks = reciprocal_space_path(sys.crystal, points_ref, 100)
    energies = range(0.0, 4.0, 400)

    # Tell Sunny how to calculate the intenties (use polarization factor).
    formfactors = [FormFactor("Mn5")] # Form factors for entangled units will be supported soon. Currently ignored.
    formula = intensity_formula(swt, :perp; kernel, formfactors)
    formula_bands = intensity_formula(swt, :perp; kernel=delta_function_kernel, formfactors)

    # Calculate the intensities and dispersion relations.
    is = intensities_broadened(swt, path, energies, formula)
    disp, is_bands = intensities_bands(swt, path, formula_bands)
end

# Plot results and compare to analytical solution.
begin
    fig = Figure(; size=(800,1200))
    ax1 = Axis(fig[1,1]; 
        xticks, 
        ylabel = "ω (meV)",
        xlabelsize = 24,
        ylabelsize = 24,
        xticklabelsize = 20,
        yticklabelsize = 20,
        title = "Ba₃Mn₂O₈",
        titlesize=26
    )
    hidexdecorations!(ax1; grid=false, ticks=false)

    # Plot dispersion curves.
    ylims!(ax1, 0.5, 3.5)
    for i in axes(disp, 2)
        if i == 1
            lines!(ax1, 1:size(disp, 1), disp[:,i]; color=:blue, label = "Sunny")
        else
            lines!(ax1, 1:size(disp, 1), disp[:,i]; color=:blue)
        end
    end

    # Calculate and plot analytical dispersion.
    d0 = disp0.(Js()..., -0.032, path)
    dm = dispm.(Js()..., -0.032, path)
    interval = 8
    scatter!(ax1, 1:interval:size(disp, 1), d0[1:interval:end]; marker=:xcross, color=:red, markersize=12.0, alpha=0.6, label="Analytical")
    scatter!(ax1, 1:interval:size(disp, 1), dm[1:interval:end]; marker=:xcross, color=:red, markersize=12.0, alpha=0.6)
    xlims!(ax1, 1, size(disp, 1))
    axislegend(ax1)


    # Plot intensities.
    ax2 = Axis(fig[2,1]; 
        xticklabelrotation=π/8, 
        ylabel = "ω (meV)",
        ylabelsize = 24,
        yticklabelsize = 20,
        titlesize=26
    )
    hidexdecorations!(ax2)
    ylims!(ax2, 0.5, 3.5)
    heatmap!(ax2, 1:size(is, 1), energies, is; colorrange=(0.0, 11.0), colormap=:gnuplot2)
    vpoints = xticks[1][2:end-1]
    for p in vpoints
        lines!(ax2, [p, p], [0.5, 3.5]; color=:white, alpha=0.5)
    end

    # Plot integrated intensities.
    ax3 = Axis(fig[3,1];
        xticks,
        xticklabelrotation=π/8, 
        xlabel = "[H, H, L] (RLU)", 
        ylabel = "Integrated Intensity (arb.)",
        xlabelsize = 24,
        ylabelsize = 24,
        xticklabelsize = 20,
        yticklabelsize = 20,
    )
    ylims!(ax3, 0.0, 75)
    xlims!(ax3, 1, size(disp, 1))
    is_integrated = sum(is, dims=2)
    is_integrated *= 65/maximum(is_integrated) # Arbitrary intensity scaling to match publication
    lines!(ax3, 1:size(is, 1), is_integrated[:])

    fig
end