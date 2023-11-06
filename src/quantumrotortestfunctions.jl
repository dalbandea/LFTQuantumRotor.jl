
import LFTSampling: analytic_force, infinitesimal_transformation, get_field, flip_momenta_sign!
function analytic_force(qrws::QuantumRotor, hmcws::AbstractHMC)
    LFTQuantumRotor.force!(qrws, hmcws)
    frc = copy(hmcws.frc)
    return frc
end
function infinitesimal_transformation(field_elem, epsilon, qrws::QuantumRotor)
    return field_elem + epsilon
end
function get_field(qrws::QuantumRotor) 
    return qrws.phi
end
function flip_momenta_sign!(hmcws::LFTQuantumRotor.QuantumRotorHMC) 
    hmcws.mom .*= -1
end


