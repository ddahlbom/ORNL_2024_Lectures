# Analytical expressions for the dispersion relation of Ba₃Mn₂O₈, as reported in
# M. Stone et al., "Singlet triplet dispersion reveals additional
# frustration...," PRL 100 (2008):
# https://doi.org/10.1103/PhysRevLett.100.237201

Js() = (1.642, 0.118, 0.256, 0.142, 0.037)

function J_of_Q(_, J1, J2, J3, J4, q)
    h, k, l = q
    ω1(h, k, l) = cos((2π/3)*(-h + k + l)) + cos((2π/3)*(-h - 2k + l)) + cos((2π/3)*(2h + k + l))
    ω2(h, k, l) = cos(2π*k) + cos(2π*(h + k)) + cos(2π*h)
    ω4(h, k, l) = cos((2π/3)*(2h - 2k + l)) + cos((2π/3)*(2h + 4k + l)) + cos((2π/3)*(-4h - 2k + l))
    -J1*ω1(h, k, l) + 2(J2 - J3)*ω2(h, k, l) - J4*ω4(h, k, l)
end

function disp0(J0, J1, J2, J3, J4, D, q)
    Δ0 = J0 + 2D/3
    sqrt(Δ0^2 + (8/3)*Δ0*J_of_Q(J0, J1, J2, J3, J4, q) )
end

function dispp(J0, J1, J2, J3, J4, D, q)
    Δp = J0 - D/3
    sqrt(Δp^2 + (8/3)*Δp*J_of_Q(J0, J1, J2, J3, J4, q) )
end

function dispm(J0, J1, J2, J3, J4, D, q)
    Δm = J0 - D/3
    sqrt(Δm^2 + (8/3)*Δm*J_of_Q(J0, J1, J2, J3, J4, q) )
end