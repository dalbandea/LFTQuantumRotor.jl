
abstract type AbstractBoundaryCondition end
abstract type PeriodicBC <: AbstractBoundaryCondition end
abstract type AntiperiodicBC <: AbstractBoundaryCondition end
abstract type OpenBC <: AbstractBoundaryCondition end

abstract type AbstractDiscretization end
abstract type ClassicalPerfectDiscretization <: AbstractDiscretization end

abstract type AbstractAngleDifferenceDiscretization <: AbstractDiscretization end
abstract type CPAngleDifferenceDiscretization <: AbstractAngleDifferenceDiscretization end
abstract type StAngleDifferenceDiscretization <: AbstractAngleDifferenceDiscretization end

abstract type AbstractAuxFields end
struct FallbackAuxField <: AbstractAuxFields end

struct QuantumRotorParm{B <: AbstractBoundaryCondition, D <: AbstractDiscretization, T} <: LFTParm
    iT::Int64
    I::Float64
    BC::Type{B}
    disc::Type{D}
    theta::T
end

function QuantumRotorParm(; iT, I, BC::Type{B} = PeriodicBC, disc::Type{D} = ClassicalPerfectDiscretization, theta = 0.0) where {B <: AbstractBoundaryCondition, D <: AbstractDiscretization}
    return QuantumRotorParm{BC, D, typeof(theta)}(iT, I, BC, disc, theta)
end


