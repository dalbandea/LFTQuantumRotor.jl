
abstract type AbstractBoundaryCondition end
abstract type PeriodicBC <: AbstractBoundaryCondition end
abstract type AntiperiodicBC <: AbstractBoundaryCondition end
abstract type OpenBC <: AbstractBoundaryCondition end

abstract type AbstractDiscretization end
abstract type ClassicalPerfectDiscretization <: AbstractDiscretization end

abstract type AbstractAngleDifferenceDiscretization <: AbstractDiscretization end
abstract type CPAngleDifferenceDiscretization <: AbstractAngleDifferenceDiscretization end
abstract type StAngleDifferenceDiscretization <: AbstractAngleDifferenceDiscretization end
abstract type StAngleDifferenceTrivMapDiscretization <: AbstractAngleDifferenceDiscretization end

struct QuantumRotorParm{B <: AbstractBoundaryCondition, D <: AbstractDiscretization} <: LFTParm
    iT::Int64
    I::Float64
    BC::Type{B}
    disc::Type{D}
end

struct QuantumRotorThetaParm{B <: AbstractBoundaryCondition, D <: AbstractDiscretization, TT} <: LFTParm
    iT::Int64
    I::Float64
    theta::TT
    BC::Type{B}
    disc::Type{D}
end

function QuantumRotorParm(; iT, I, BC::Type{B} = PeriodicBC, disc::Type{D} = ClassicalPerfectDiscretization) where {B <: AbstractBoundaryCondition, D <: AbstractDiscretization}
    return QuantumRotorParm{BC, D}(iT, I, BC, disc)
end

function QuantumRotorThetaParm(; iT, I, BC::Type{B} = PeriodicBC, disc::Type{D} = ClassicalPerfectDiscretization) where {B <: AbstractBoundaryCondition, D <: AbstractDiscretization}
    return QuantumRotorParm{BC, D}(iT, I, BC, disc)
end

