# Model of Ba₃Mn₂O₈ using parameters reported in M. Stone et al., "Singlet
# triplet dispersion reveals additional frustration...," PRL 100 (2008):
# https://doi.org/10.1103/PhysRevLett.100.237201

function system_for_contraction(; dims=(4, 4, 4), g=1.0,
    D = -0.032, 
    J0 = 1.642, 
    J1 = 0.118, 
    J2 = 0.256, 
    J3 = 0.142, 
    J4 = 0.037
)

    # Crystal and system
    crystal = Crystal(datadir("Ba3Mn2O8_OCD_2008132.cif"); symprec=0.01)
    crystal_Mn = subcrystal(crystal, "Mn1")
    sys = System(crystal_Mn, dims, [SpinInfo(1; S=1, g)], :SUN)

    # Exchange
    S = spin_matrices(1)
    S1, S2 = to_product_space(S, S)
    set_pair_coupling!(sys, J0 * (S1' * S2), Bond(1, 2, [0, 0, 0]))
    set_pair_coupling!(sys, J1 * (S1' * S2), Bond(3, 2, [0, 0, 0]))
    set_pair_coupling!(sys, J2 * (S1' * S2), Bond(1, 1, [1, 0, 0]))
    set_pair_coupling!(sys, J3 * (S1' * S2), Bond(1, 2, [1, 0, 0]))
    set_pair_coupling!(sys, J4 * (S1' * S2), Bond(3, 2, [1, 0, 0]))

    # Single-ion anisotropy
    set_onsite_coupling!(sys, D*S[3]^2, 1)

    return sys
end