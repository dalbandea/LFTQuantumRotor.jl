import LFTSampling: sampler, copy!

struct QuantumRotorWorkspace{T, N, P <: LFTParm} <: QuantumRotor
    PRC::Type{T}
    phi::Array{T, N}
    params::P
    function QuantumRotorWorkspace(::Type{T}, lp::QuantumRotorParm) where {T}
        phi = Array{T, 1}(undef, lp.iT)
        return new{T, 1, typeof(lp)}(T, phi, lp)
    end
end

function (::Type{QuantumRotor})(::Type{T} = Float64; kwargs...) where {T}
    return QuantumRotorWorkspace(T, QuantumRotorParm(;kwargs...))
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

function winding!(qrws::QuantumRotor, disc::Type{AbstractDiscretization})
    for i in 1:qrws.params.iT
        qrws.phi[i] = qrws.phi[i] + (i-1) * 2pi/qrws.params.iT
    end
end

function antiwinding!(qrws::QuantumRotor, disc::Type{AbstractDiscretization})
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
