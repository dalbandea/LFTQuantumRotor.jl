import LFTSampling: action

action_discretization_factor(::Type{AbstractDiscretization}) = 0.5
action_discretization_factor(::Type{StAngleDifferenceDiscretization}) = 1.0
action_discretization_factor(::Type{CPAngleDifferenceDiscretization}) = 0.5

theta_term(qrws::QuantumRotor) = theta_term(qrws, qrws.params.disc)
function theta_term(qrws::QuantumRotor, ::Type{D}) where D <: AbstractDiscretization
    if qrws.params.theta == 0.0
        return 0.0
    else
        return top_charge(qrws)*qrws.params.theta
    end
end


function action(qrws::QuantumRotor)
    S = zero(qrws.PRC)

    for t in 1:qrws.params.iT-1
        S += action_t(qrws, t)
    end
    S += boundary_action(qrws)

    return qrws.params.I * S * action_discretization_factor(qrws.params.disc) + theta_term(qrws)
end

action_t(qrws::QuantumRotor, t::Int64) = action_t(qrws, qrws.params.disc, t)

function action_t(qrws::QuantumRotor, disc::Type{CPAngleDifferenceDiscretization}, t::Int64)
    ds = Mod(qrws.phi[t], 2pi)^2
    return ds
end

function action_t(qrws::QuantumRotor, disc::Type{StAngleDifferenceDiscretization}, t::Int64)
    ds = 1-cos(qrws.phi[t])
    return ds
end


boundary_action(qrws::QuantumRotor) = boundary_action(qrws, qrws.params.disc, qrws.params.BC)
boundary_action(qrws::QuantumRotor, disc::Type{CPAngleDifferenceDiscretization}, BC::Type{OpenBC}) = zero(qrws.PRC)
boundary_action(qrws::QuantumRotor, disc::Type{StAngleDifferenceDiscretization}, BC::Type{OpenBC}) = zero(qrws.PRC)
function boundary_action(qrws::QuantumRotor, disc::Type{CPAngleDifferenceDiscretization}, BC::Type{PeriodicBC})
    qt = zero(qrws.PRC)
    for t in 1:qrws.params.iT-1
        qt -= qrws.phi[t]
    end
    return Mod(qt, 2pi)^2
end

function boundary_action(qrws::QuantumRotor, disc::Type{StAngleDifferenceDiscretization}, BC::Type{PeriodicBC})
    qt = zero(qrws.PRC)
    for t in 1:qrws.params.iT-1
        qt -= qrws.phi[t]
    end
    return 1 - cos(qt)
end


# function action_t(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, t::Int64)
#     tu = right(qrws, t)
#     ds = Mod(qrws.phi[tu] - qrws.phi[t], 2pi)^2
#     return ds * qrws.params.I  / 2
# end


# """
#     daction_t

# Compute action of points coupled to the angle `t`.
# """
# daction_t(qrws::QuantumRotor, t::Int64) = daction_t(qrws, qrws.params.disc, t)
# function daction_t(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, t::Int64)
#     td = left(qrws, t)
#     ds = action_t(qrws, td) + action_t(qrws, t)
#     return ds
# end

# right(qrws::QuantumRotor, t) = right(qrws, qrws.params.BC, t)
# function right(qrws::QuantumRotor, BC::Type{B}, t::Int64) where B <: AbstractBoundaryCondition 
#     error("Function right not implemented for boundary condition $B")
# end
# right(qrws::QuantumRotor, BC::Type{PeriodicBC}, t::Int64) = mod1(t+1,
#                                                                  qrws.params.iT)

# left(qrws::QuantumRotor, t) = left(qrws, qrws.params.BC, t)
# function left(qrws::QuantumRotor, BC::Type{B}, t::Int64) where B <: AbstractBoundaryCondition 
#     error("Function left not implemented for boundary condition $B")
# end
# left(qrws::QuantumRotor, BC::Type{PeriodicBC}, t::Int64) = mod1(t-1,
                                                                 # qrws.params.iT)

# boundary_action(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, BC::Type{PeriodicBC}) = Mod(qrws.phi[1] - qrws.phi[end], 2pi)^2
# boundary_action(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, BC::Type{AntiperiodicBC}) = Mod(-qrws.phi[1] - qrws.phi[end], 2pi)^2
# boundary_action(qrws::QuantumRotor, disc::Type{ClassicalPerfectDiscretization}, BC::Type{OpenBC}) = zero(qrws.PRC)



