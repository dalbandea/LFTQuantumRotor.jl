import LFTSampling: sampler, copy!

struct QuantumRotorWorkspace{T1, T2, N, P <: LFTParm, AUX <: AbstractAuxFields} <: QuantumRotor
    PRC::Type{T1}
    phi::Array{T2, N}
    params::P
    aux::AUX
    function QuantumRotorWorkspace(::Type{T1}, ::Type{T2}=T1; aux::T3 = FallbackAuxField(), params::QuantumRotorParm) where {T1,T2,T3<:AbstractAuxFields}
        phi = Array{T2, 1}(undef, params.iT)
        return new{T1,T2, 1, typeof(params), typeof(aux)}(T1, phi, params, aux)
    end
end

function (::Type{QuantumRotor})(::Type{T1} = Float64, ::Type{T2} = T1; aux::T3 = FallbackAuxField(), kwargs...) where {T1,T2,T3<:AbstractAuxFields}
    return QuantumRotorWorkspace(T1, T2, aux = aux, params=QuantumRotorParm(;kwargs...))
end

struct QuantumRotorHMC{A <: AbstractArray} <: AbstractHMC
    params::HMCParams
    frc::A
    mom::A
end

function QuantumRotorHMC(qrws::QuantumRotor, hmcp::HMCParams)
    frc = similar(qrws.phi)
    mom = similar(qrws.phi)
    return QuantumRotorHMC{typeof(frc)}(hmcp, frc, mom)
end

sampler(lftws::QuantumRotor, hmcp::HMCParams) = QuantumRotorHMC(lftws, hmcp)

function copy!(qrws_dst::QuantumRotor, qrws_src::QuantumRotor)
    qrws_dst.phi .= qrws_src.phi
    return nothing
end

function randomize!(qrws::QuantumRotor)
    qrws.phi .= 2pi*Random.rand(qrws.PRC, size(qrws.phi)...) .- pi
    return nothing
end

function coldstart!(qrws::QuantumRotor)
    qrws.phi .= one(qrws.PRC)
end

#custom mod
#changes mod domain from (0 to z) to (-z/2 to z/2)
function Mod(x, z)

    a = mod(x,z)
    if a <= z/2
        return a
    else
        return a - z
    end

end

winding!(qrws::QuantumRotor) = winding!(qrws, qrws.params.disc)
winding!(qrws::QuantumRotor, L::Integer) = winding!(qrws, qrws.params.disc, L)
antiwinding!(qrws::QuantumRotor) = antiwinding!(qrws, qrws.params.disc)

function winding!(qrws::QuantumRotor, disc::Type{D}) where D <: AbstractDiscretization
    for i in 1:qrws.params.iT
        qrws.phi[i] = qrws.phi[i] + (i-1) * 2pi/qrws.params.iT
    end
end

function antiwinding!(qrws::QuantumRotor, disc::Type{D}) where D <: AbstractDiscretization
    for i in 1:qrws.params.iT
        qrws.phi[i] = qrws.phi[i] - (i-1) * 2pi/qrws.params.iT
    end
end

function winding!(qrws::QuantumRotor, disc::Type{D}) where D <: AbstractAngleDifferenceDiscretization
    for i in 1:qrws.params.iT
        qrws.phi[i] = qrws.phi[i] + 2pi/qrws.params.iT
    end
end

function winding!(qrws::QuantumRotor, disc::Type{D}, L::Integer) where D <: AbstractAngleDifferenceDiscretization
    for i in 1:L
        qrws.phi[i] = qrws.phi[i] + 2pi/L
    end
end

function antiwinding!(qrws::QuantumRotor, disc::Type{D}) where D <: AbstractAngleDifferenceDiscretization
    for i in 1:qrws.params.iT
        qrws.phi[i] = qrws.phi[i] - 2pi/qrws.params.iT
    end
end


function unwind!(qrws::QuantumRotor)
    i = 0.0
    while LFTs.top_charge(qrws) |> round != 0.0
        Q = LFTs.top_charge(qrws)
        if Q < 0
            LFTs.winding!(qrws)
        elseif Q > 0
            LFTs.antiwinding!(qrws)
        end
        if i == 10000
            println("Stopped unwinding after $i iterations")
            break
        end
        i+=1
    end
    return nothing
end
