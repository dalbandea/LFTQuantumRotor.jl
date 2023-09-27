
function top_charge(qrws::QuantumRotor)
    Q = 0.0

    for i in 1:qrws.params.iT-1
        Q += topcharge_t(qrws, i)
    end

    Q += boundary_topcharge(qrws)

    return Q/2pi
end


topcharge_t(qrws::QuantumRotor, t::Int64) = topcharge_t(qrws, qrws.params.disc, t)
function topcharge_t(qrws::QuantumRotor, disc::Type{D}, t::Int64) where D <: AbstractAngleDifferenceDiscretization
    return Mod(qrws.phi[t], 2pi)
end

boundary_topcharge(qrws::QuantumRotor) = boundary_topcharge(qrws, qrws.params.disc, qrws.params.BC)
function boundary_topcharge(qrws::QuantumRotor, disc::Type{D}, BC::Type{PeriodicBC}) where D <: AbstractAngleDifferenceDiscretization
    qt = zero(qrws.PRC)
    for t in 1:qrws.params.iT-1
        qt -= qrws.phi[t]
    end
    return Mod(qt, 2pi)
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
