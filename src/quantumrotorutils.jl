
"""
    angdiff_to_ang(dconf)

Takes model `dconf` with BC `StAngleDifferenceDiscretization` and returns
`aconf` with `StadardDiscretization`, i.e. an absolute angle version with
`phi[1]=0.0`, such that `action(dconf)` equals `action(aconf)`.
"""
function angdiff_to_ang(dconf)
    aconf = QuantumRotor(I = dconf.params.I, iT = dconf.params.iT, BC =
                         dconf.params.BC, disc = StandardDiscretization)
    aconf.phi[1] = 0.0
    for i in 2:aconf.params.iT
        aconf.phi[i] = aconf.phi[i-1] + dconf.phi[i-1]
    end
    return aconf
end
