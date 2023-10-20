import LFTSampling: sample!, generate_momenta!, Hamiltonian, update_momenta!, update_fields!

function winding_step!(qrws::QuantumRotor)
    ws_cp = deepcopy(qrws)

    sini = action(qrws)

    r = rand()

    if r > 0.5
        winding!(qrws)
    else
        antiwinding!(qrws)
    end

    sfin = action(qrws)

    ds = sfin - sini

    LFTSampling.metropolis_accept_reject!(qrws, ws_cp, ds)

    return nothing
end

###########################
# Abstract Discretization #
###########################

function sample!(qrws::QuantumRotor, samplerws::AbstractHMC; do_winding = false)
    hmc!(qrws, samplerws)
    do_winding && winding_step!(qrws)
    return nothing
end

generate_momenta!(qrws::QuantumRotor, hmcws::QuantumRotorHMC) = generate_momenta!(qrws, hmcws, qrws.params.disc)
function generate_momenta!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{D}) where D <: AbstractDiscretization
    for i in 1:length(hmcws.mom)
        hmcws.mom[i] = randn(qrws.PRC) * hmcws.params.width
    end
    return nothing
end

function Hamiltonian(qrws::QuantumRotor, hmcws::QuantumRotorHMC)
    H = mapreduce(x -> x^2, +, hmcws.mom)/2.0/hmcws.params.width^2 + action(qrws)
    return H
end

function force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC)
    return force!(qrws, hmcws, qrws.params.disc, qrws.params.BC)
end

function theta_force(qrws::QuantumRotor, BC::Type{OpenBC}) 
    if qrws.params.theta == 0.0
        return 0.0
    else
        return one(qrws.PRC)*qrws.params.theta
    end
end

function update_momenta!(qrws::QuantumRotor, epsilon, hmcws::QuantumRotorHMC)

    # Load phi force
    force!(qrws, hmcws) 

    # Update phi momenta
    hmcws.mom .= hmcws.mom .+ epsilon .* hmcws.frc

    return nothing
end

update_fields!(qrws::QuantumRotor, epsilon, hmcws::QuantumRotorHMC) = update_fields!(qrws, epsilon, hmcws, qrws.aux)
function update_fields!(qrws::QuantumRotor, epsilon, hmcws::QuantumRotorHMC, aux::AbstractAuxFields)
    # Update phi field
    qrws.phi .= qrws.phi .+ epsilon .* hmcws.mom ./ hmcws.params.width^2
    return nothing
end

function flip_momenta_sign!(hmcws::QuantumRotorHMC)
    hmcws.mom .= .- hmcws.mom
    return nothing
end




###################################
# Angle Difference Discretization #
###################################

function generate_momenta!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceDiscretization})
    for i in 1:length(hmcws.mom)-1
        hmcws.mom[i] = randn(qrws.PRC) * hmcws.params.width
    end
    hmcws.mom[qrws.params.iT] = zero(qrws.PRC)
    return nothing
end

function force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{D}, BC::Type{B}) where {D <: AbstractAngleDifferenceDiscretization, B <: AbstractBoundaryCondition}

    for t in 1:qrws.params.iT-1
        hmcws.frc[t] = force_t(qrws, t, disc)
    end

    if B == PeriodicBC && D == StAngleDifferenceDiscretization
        sumphi = @views sum(qrws.phi[1:end-1])
        hmcws.frc .-= qrws.params.I * sin(sumphi) 
    end

    boundary_force!(qrws, hmcws, disc, BC)
    
    return nothing
end

function force_t(qrws::QuantumRotor, t::Int64, disc::Type{StAngleDifferenceDiscretization})
    return -qrws.params.I * sin(qrws.phi[t]) - theta_force(qrws, qrws.params.BC)
end

function force_t(qrws::QuantumRotor, t::Int64, disc::Type{CPAngleDifferenceDiscretization})
    return -qrws.params.I * Mod(qrws.phi[t],2pi) - theta_force(qrws, qrws.params.BC)
end

function boundary_force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{D}, BC::Type{OpenBC}) where D <: AbstractAngleDifferenceDiscretization
    hmcws.frc[qrws.params.iT] = zero(qrws.PRC)
    return nothing
end

function boundary_force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceDiscretization}, BC::Type{PeriodicBC})
    hmcws.frc[qrws.params.iT] = zero(qrws.PRC)
    return nothing
end

