using Test
using LFTSampling
using LFTQuantumRotor

@testset verbose = true "HMC tests" begin
    include("hmctests.jl")
end
