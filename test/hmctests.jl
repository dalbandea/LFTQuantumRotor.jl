I = 5.0
iT = 100

for BC in [PeriodicBC, OpenBC]
    model = QuantumRotor(
                         Float64, 
                         I = I, 
                         iT = iT, 
                         BC = BC, 
                         disc = StAngleDifferenceDiscretization,
                        )

    randomize!(model)

    smplr = HMC(integrator = Leapfrog(1.0, 20))
    samplerws = LFTSampling.sampler(model, smplr)
    LFTSampling.generate_momenta!(model, samplerws)


    @testset "$(model.params.BC) HMC reversibility" begin
        model_bckp = deepcopy(model)
        LFTSampling.reversibility!(model, samplerws)
        dphi = model.phi .- model_bckp.phi
        @test isapprox(zero(model.PRC), mapreduce(x -> abs2(x), +, dphi), atol = 1e-15)
    end

    @testset verbose = true "$(model.params.BC) HMC force" begin
        dF = LFTSampling.force_test(model, samplerws, 1e-6)
        @test isapprox(zero(model.PRC), dF, atol = 1e-5)
    end

end


for BC in [PeriodicBC, OpenBC]
    model = QuantumRotor(
                         Float64, 
                         I = I, 
                         iT = iT, 
                         BC = BC, 
                         disc = StAngleDifferenceDiscretization,
                         theta = 1.0
                        )

    randomize!(model)

    smplr = HMC(integrator = Leapfrog(1.0, 20))
    samplerws = LFTSampling.sampler(model, smplr)
    LFTSampling.generate_momenta!(model, samplerws)


    @testset "θ=1.0 $(model.params.BC) HMC reversibility" begin
        model_bckp = deepcopy(model)
        LFTSampling.reversibility!(model, samplerws)
        dphi = model.phi .- model_bckp.phi
        @test isapprox(zero(model.PRC), mapreduce(x -> abs2(x), +, dphi), atol = 1e-15)
    end

    @testset verbose = true "θ=1.0 $(model.params.BC) HMC force" begin
        dF = LFTSampling.force_test(model, samplerws, 1e-6)
        @test isapprox(zero(model.PRC), dF, atol = 1e-5)
    end

end


for BC in [PeriodicBC, OpenBC]
    model = QuantumRotor(
                         Float64, 
                         I = I, 
                         iT = iT, 
                         BC = BC, 
                         disc = StandardDiscretization,
                        )

    randomize!(model)

    smplr = HMC(integrator = Leapfrog(1.0, 20))
    samplerws = LFTSampling.sampler(model, smplr)
    LFTSampling.generate_momenta!(model, samplerws)


    @testset "StDisc $(model.params.BC) HMC reversibility" begin
        model_bckp = deepcopy(model)
        LFTSampling.reversibility!(model, samplerws)
        dphi = model.phi .- model_bckp.phi
        @test isapprox(zero(model.PRC), mapreduce(x -> abs2(x), +, dphi), atol = 1e-15)
    end

    @testset verbose = true "StDisc $(model.params.BC) HMC force" begin
        dF = LFTSampling.force_test(model, samplerws, 1e-6)
        @test isapprox(zero(model.PRC), dF, atol = 1e-5)
    end

end


for BC in [PeriodicBC, OpenBC]
    model = QuantumRotor(
                         Float64, 
                         I = I, 
                         iT = iT, 
                         BC = BC, 
                         disc = StandardDiscretization,
                         theta = 1.0
                        )

    randomize!(model)

    smplr = HMC(integrator = Leapfrog(1.0, 20))
    samplerws = LFTSampling.sampler(model, smplr)
    LFTSampling.generate_momenta!(model, samplerws)

    @testset "StDisc $(model.params.BC) θ=1.0 HMC reversibility" begin
        model_bckp = deepcopy(model)
        LFTSampling.reversibility!(model, samplerws)
        dphi = model.phi .- model_bckp.phi
        @test isapprox(zero(model.PRC), mapreduce(x -> abs2(x), +, dphi), atol = 1e-15)
    end

    @testset verbose = true "StDisc $(model.params.BC) θ=1.0 HMC force" begin
        dF = LFTSampling.force_test(model, samplerws, 1e-6)
        @test isapprox(zero(model.PRC), dF, atol = 1e-5)
    end

end
