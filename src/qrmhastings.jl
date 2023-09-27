Base.@kwdef mutable struct UnboundedQRPBCDistribution <: AbstractDistribution
    I::Float64
    T::Float64
    phi1dist = Distributions.Normal(0.0, 1/sqrt(I*T))
    phiidist = Distributions.Normal(0.0, 1/sqrt(I))
    winding_n = zero(Int64)
    winding_n_backup = zero(Int64)
    S = orthonormalizing_base(T)
end
export UnboundedQRPBCDistribution


# Interface functions

function sample!(qrws::QuantumRotor, dist::AbstractDistribution)
    qrws.phi .= [sample(dist) for _ in 1:qrws.params.iT]
    # sample(qrws.phi, samplerws.params.dist, qrws)
    return nothing
end

function log_prob(qrws::QuantumRotor, dist::AbstractDistribution)
    return @views sum(log_prob(dist, qrws.phi)[1:qrws.params.iT-1])
end

# function log_prob(qrws::QuantumRotor, samplerws::AbstractMetropolisHastings) 
#     log_prob(qrws, samplerws.params.dist)
#     return nothing
# end


## QuantumRotorOBCDistribution

function sample(dist::QuantumRotorOBCDistribution)
    return random_OBC_CP(dist.I)
end


function log_prob(dist::QuantumRotorOBCDistribution, x)
    logpdf = logprob_OBC_CP(dist.I)
    return logpdf(x)
end

function log_prob(dist::QuantumRotorOBCDistribution, x::Vector)
    logpdf = logprob_OBC_CP(dist.I)
    return logpdf.(x)
end

## UnboundedQRPBCDistribution

function sample!(qrws::QuantumRotor, dist::UnboundedQRPBCDistribution)
    psi = zeros(qrws.params.iT-1)
    psi[1] = rand(dist.phi1dist)
    for i in 2:qrws.params.iT-1
        psi[i] = rand(dist.phiidist)
    end
    qrws.phi[1:end-1] .= dist.S * psi
    dist.winding_n = Random.rand(0:0)
    do_windings!(qrws, dist)
    return nothing
end

function log_prob(qrws::QuantumRotor, dist::UnboundedQRPBCDistribution) 
    undo_windings!(qrws, dist)
    logp = log_prob(qrws.phi, dist::UnboundedQRPBCDistribution) 
    do_windings!(qrws, dist)
    return logp
end

# function log_prob(qrws::QuantumRotor, dist::UnboundedQRPBCDistribution) 
#     psi = transpose(dist.S) * qrws.phi[1:qrws.params.iT-1]
#     res = Distributions.logpdf(dist.phi1dist, psi[1])
#     # res = -1/2*psi[1]^2*qrws.params.I*qrws.params.iT
#     for i in 2:qrws.params.iT-1
#         res += Distributions.logpdf(dist.phiidist, psi[i])
#         # res += -1/2*psi[i]^2*qrws.params.I
#     end
#     return res
# end

function log_prob(phi, dist::UnboundedQRPBCDistribution) 
    T = size(phi,1)
    psi = transpose(dist.S) * phi[1:T-1]
    res = Distributions.logpdf(dist.phi1dist, psi[1])
    # res = -1/2*psi[1]^2*qrws.params.I*qrws.params.iT
    for i in 2:T-1
        res += Distributions.logpdf(dist.phiidist, psi[i])
        # res += -1/2*psi[i]^2*qrws.params.I
    end
    return res
end

function reset_sampler!(dist::UnboundedQRPBCDistribution, accepted)
    if accepted == 0
        dist.winding_n = dist.winding_n_backup
    end
end

function do_windings!(qrws::QuantumRotor, dist::UnboundedQRPBCDistribution)
    if dist.winding_n > 0
        for i in 1:dist.winding_n
            winding!(qrws)
        end
    elseif dist.winding_n < 0
        for i in 1:(-dist.winding_n)
            antiwinding!(qrws)
        end
    end
    return nothing
end

function undo_windings!(qrws::QuantumRotor, dist::UnboundedQRPBCDistribution)
    # dist.winding_n = round(top_charge(qrws))
    dist.winding_n *= -1
    do_windings!(qrws,dist)
    dist.winding_n *= -1
    return nothing
end


# Low-level functions

## QuantumRotorOBCDistribution

function random_OBC_CP(I)
    xi = rand()
    return sqrt(2/I) * erfinv((2*xi-1) * erf(sqrt(I/2)*pi))
end

function PDF_OBC_CP(I)
    Z = sqrt(2pi/I) * erf(sqrt(I/2)*pi)
    return x -> exp(-I/2 * x^2) / Z
end

function logprob_OBC_CP(I)
    Z = sqrt(2pi/I) * erf(sqrt(I/2)*pi)
    return x -> -I/2 * x^2 #- log(Z)
end

function PDF_configuration_OBC_CP(phi, I, T)
    factorized_PDF = PDF_OBC_CP(I)
    res = 1.0
    for i in 1:T
        res = res * factorized_PDF(phi[i])
    end
    return res
end

function logprob_configuration_OBC_CP(phi, I, T)
    factorized_logprob = logprob_OBC_CP(I)
    res = 0.0
    for i in 1:T
        res += factorized_logprob(phi[i])
    end
    return res
end

## UnboundedQRPBCDistribution

function orthonormalizing_base(T)
    S = zeros(T-1,T-1)
    S[1,:] .= -1.0
    S[:,1] .= 1.0
    for i in 2:T-1
        S[i,i] = 1.0
    end
    return copy(LinearAlgebra.qr(S).Q)
end
