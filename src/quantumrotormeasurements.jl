
function top_charge(qrws::QuantumRotor)
    Q = 0.0

    for i in 1:qrws.params.iT-1
        if qrws.params.theta == 0.0
            Q += topcharge_t(qrws, i)
        else
            Q += cinf_topcharge_t(qrws,i)
        end
    end

    if qrws.params.theta == 0.0
        Q += boundary_topcharge(qrws)
    else
        Q += boundary_cinf_topcharge(qrws)
    end

    return Q/2pi
end


topcharge_t(qrws::QuantumRotor, t::Int64) = topcharge_t(qrws, qrws.params.disc, t)
function topcharge_t(qrws::QuantumRotor, disc::Type{D}, t::Int64) where D <: AbstractAngleDifferenceDiscretization
    return Mod(qrws.phi[t], 2pi)
end

cinf_topcharge_t(qrws::QuantumRotor, t::Int64) = cinf_topcharge_t(qrws, t, qrws.params.disc)
function cinf_topcharge_t(qrws::QuantumRotor, t::Int64, disc::Type{D}) where D <: AbstractAngleDifferenceDiscretization
    return -im*log(exp(im*qrws.phi[t]))
end

boundary_topcharge(qrws::QuantumRotor) = boundary_topcharge(qrws, qrws.params.disc, qrws.params.BC)
boundary_topcharge(qrws::QuantumRotor, disc::Type{D}, BC::Type{OpenBC}) where D <: AbstractAngleDifferenceDiscretization = zero(qrws.PRC)
function boundary_topcharge(qrws::QuantumRotor, disc::Type{D}, BC::Type{PeriodicBC}) where D <: AbstractAngleDifferenceDiscretization
    qt = zero(qrws.PRC)
    for t in 1:qrws.params.iT-1
        qt -= qrws.phi[t]
    end
    return Mod(qt, 2pi)
end

boundary_cinf_topcharge(qrws::QuantumRotor) = cinf_boundary_topcharge(qrws, qrws.params.disc, qrws.params.BC)
cinf_boundary_topcharge(qrws::QuantumRotor, disc::Type{D}, BC::Type{OpenBC}) where D <: AbstractAngleDifferenceDiscretization = zero(qrws.PRC)
function cinf_boundary_topcharge(qrws::QuantumRotor, disc::Type{D}, BC::Type{PeriodicBC}) where D <: AbstractAngleDifferenceDiscretization
    qt = zero(qrws.PRC)
    for t in 1:qrws.params.iT-1
        qt -= qrws.phi[t]
    end
    return -im*log(exp(im*qt))
end

diff_top_charge(qrws::LFTQuantumRotor.QuantumRotor) = diff_top_charge(qrws, qrws.params.disc, qrws.params.BC)
function diff_top_charge(qrws::LFTQuantumRotor.QuantumRotor, disc::Type{D}, ::Type{BC}) where {D <: LFTQuantumRotor.AbstractAngleDifferenceDiscretization, BC <: LFTQuantumRotor.AbstractBoundaryCondition}
    Q = zero(eltype(qrws.phi))
    qt = zero(eltype(qrws.phi))
    for i in 1:qrws.params.iT-1
        qt -= qrws.phi[i]
        Q += sin(qrws.phi[i])
    end

    if BC == LFTQuantumRotor.PeriodicBC
        Q += sin(qt) 
    elseif BC == LFTQuantumRotor.OpenBC
    else
        error("This should not happen")
    end
    return Q/2pi
end

function diff_top_charge(qrws::LFTQuantumRotor.QuantumRotor, disc::Type{StandardDiscretization}, ::Type{BC}) where BC <: AbstractBoundaryCondition
    Q = zero(eltype(qrws.phi))
    for i in 1:qrws.params.iT-1
        Q += sin(qrws.phi[i+1]-qrws.phi[i])
    end

    if BC == LFTQuantumRotor.PeriodicBC
        Q += sin(qrws.phi[1]-qrws.phi[qrws.params.iT]) 
    end

    return Q/2pi
end


# abstract type Susceptibility <: AbstractObservable end

# struct QRMasterObs end

# mutable struct QRTopologicalCharge <: AbstractScalar
#     name::String
#     ID::String
#     filepath::String
#     result::Float64
#     history::Vector{Float64}
#     function QRTopologicalCharge(; wdir::String = "./results/trash/", 
#                               name::String = "Topological charge", 
#                               ID::String = "topcharge", 
#                               mesdir::String = "measurements/", 
#                               extension::String = ".txt")
#         filepath = joinpath(wdir, mesdir, ID*extension)
#         mkpath(joinpath(wdir,mesdir))
#         result = zero(Float64)
#         history = Vector{Float64}()
#         return new(name, ID, filepath, result, history)
#     end
# end
# export QRTopologicalCharge


# mutable struct QRMFTopologicalCharge <: AbstractCorrelator
#     name::String
#     ID::String
#     filepath::String
#     result::Vector{Float64}
#     history::Vector{Vector{Float64}}
#     R::Int64
#     function QRMFTopologicalCharge(R::Int64; wdir::String = "./results/trash/", 
#                               name::String = "Topological charge", 
#                               ID::String = "MFtopcharge", 
#                               mesdir::String = "measurements/", 
#                               extension::String = ".txt")
#         filepath = joinpath(wdir, mesdir, ID*extension)
#         mkpath(joinpath(wdir,mesdir))
#         result = Vector{Float64}()
#         history = Vector{Vector{Float64}}()
#         return new(name, ID, filepath, result, history, R)
#     end
# end
# export QRMFTopologicalCharge

# function top_charge(qrws::QuantumRotor)
#     Q = 0.0

#     for i in 1:qrws.params.iT-1
#         Q += Mod(qrws.phi[i+1]-qrws.phi[i], 2pi)
#     end

#     Q += Mod(qrws.phi[1]-qrws.phi[end], 2pi)

#     return Q/2pi
# end

# function (obs::QRTopologicalCharge)(qrws::QuantumRotor)
#     obs.result = top_charge(qrws)
#     return nothing
# end


# function top_charge_density(qrws::QuantumRotor, i::Int64)
#     iu = mod1(i+1, qrws.params.iT)
#     Q  = Mod(qrws.phi[iu]-qrws.phi[i], 2pi)

#     return Q/2pi
# end

# function suscep_R(qrws::QuantumRotor, R::Int64, t0::Int64)
#     q0 = top_charge_density(qrws, t0)
#     res = zero(Float64)
#     for t in t0-R:t0+R
#         qt = top_charge_density(qrws, mod1(t, qrws.params.iT))
#         res += qt * q0
#     end
#     return res
# end

# function (obs::QRMFTopologicalCharge)(qrws::QuantumRotor)
#     obs.result = similar(qrws.phi)
#     obs.result .= zero(qrws.PRC)
#     for t0 in 1:qrws.params.iT
#         obs.result[t0] = suscep_R(qrws, obs.R, t0)
#     end
# end

# # function (obs::QRMFTopologicalCharge)(qrws::QuantumRotor)
# #     obs.result = top_charge(qrws)
# #     return nothing
# # end
