using Test
using LFTSampling
using LFTQuantumRotor

@testset verbose = true "Quantum Rotor tests" begin
    @testset verbose = true "Model tests" begin
        include("modeltests.jl")
    end

    @testset verbose = true "HMC tests" begin
        include("hmctests.jl")
    end

    @testset verbose = true "I/O tests" begin
        include("iotests.jl")
    end
end
