import LFTSampling: sample!, generate_momenta!, Hamiltonian, update_momenta!, update_fields!

function sample!(qrws::QuantumRotor, samplerws::AbstractHMC; do_winding = false)
    hmc!(qrws, samplerws)
    do_winding && winding_step!(qrws)
    return nothing
end

generate_momenta!(qrws::QuantumRotor, hmcws::QuantumRotorHMC) = generate_momenta!(qrws, hmcws, qrws.params.disc)
function generate_momenta!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{AbstractDiscretization})
    for i in 1:length(hmcws.mom)
        hmcws.mom[i] = randn() * hmcws.params.width
    end
    return nothing
end

function generate_momenta!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceDiscretization})
    for i in 1:length(hmcws.mom)-1
        hmcws.mom[i] = randn() * hmcws.params.width
    end
    hmcws.mom[qrws.params.iT] = zero(qrws.PRC)
    return nothing
end

function generate_momenta!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceTrivMapDiscretization})
    for i in 1:length(hmcws.mom)-1
        hmcws.mom[i] = randn() * hmcws.params.width
    end
    hmcws.mom[qrws.params.iT] = zero(qrws.PRC)
    return nothing
end

function Hamiltonian(qrws::QuantumRotor, hmcws::QuantumRotorHMC)
    H = mapreduce(x -> x^2, +, hmcws.mom)/2.0/hmcws.params.width^2 + action(qrws)
    return H
end

function force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC)
    return force!(qrws, hmcws, qrws.params.disc, qrws.params.BC)
end

function force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceDiscretization}, BC::Type{OpenBC})

    for t in 1:qrws.params.iT-1
        hmcws.frc[t] = -qrws.params.I * sin(qrws.phi[t])
    end

    boundary_force!(qrws, hmcws, disc, BC)
    
    return nothing
end

function boundary_force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceDiscretization}, BC::Type{OpenBC})
    hmcws.frc[qrws.params.iT] = zero(qrws.PRC)
    return nothing
end

function force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceDiscretization}, BC::Type{PeriodicBC})

    sumphi = @views sum(qrws.phi[1:end-1])
    for t in 1:qrws.params.iT-1
        hmcws.frc[t] = -qrws.params.I * (sin(qrws.phi[t]) + sin(sumphi))
    end

    boundary_force!(qrws, hmcws, disc, BC)
    
    return nothing
end

function boundary_force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceDiscretization}, BC::Type{PeriodicBC})
    hmcws.frc[qrws.params.iT] = zero(qrws.PRC)
    return nothing
end

function force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceTrivMapDiscretization}, BC::Type{PeriodicBC})

    sumphi = @views sum(qrws.phi[1:end-1])
    for t in 1:qrws.params.iT-1
        hmcws.frc[t] = -qrws.params.I / sqrt(10*qrws.params.I) * (sin(qrws.phi[t]/sqrt(10*qrws.params.I)) + sin(sumphi/sqrt(10*qrws.params.I)))
    end

    boundary_force!(qrws, hmcws, disc, BC)
    
    return nothing
end

function boundary_force!(qrws::QuantumRotor, hmcws::QuantumRotorHMC, disc::Type{StAngleDifferenceTrivMapDiscretization}, BC::Type{PeriodicBC})
    hmcws.frc[qrws.params.iT] = zero(qrws.PRC)
    return nothing
end


function update_momenta!(qrws::QuantumRotor, epsilon, hmcws::QuantumRotorHMC)

    # Load phi force
    force!(qrws, hmcws) 

    # Update phi momenta
    hmcws.mom .= hmcws.mom .+ epsilon .* hmcws.frc

    return nothing
end


function update_fields!(qrws::QuantumRotor, epsilon, hmcws::QuantumRotorHMC)
    # Update phi field
    qrws.phi .= qrws.phi .+ epsilon .* hmcws.mom ./ hmcws.params.width^2
    return nothing
end

function flip_momenta_sign!(hmcws::QuantumRotorHMC)
    hmcws.mom .= .- hmcws.mom
    return nothing
end
