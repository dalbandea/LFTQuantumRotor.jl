
I = 5.0
iT = 10

modeldiff = QuantumRotor(I = I, iT = iT, BC = PeriodicBC, disc = StAngleDifferenceDiscretization)

randomize!(modeldiff)

model = LFTQuantumRotor.angdiff_to_ang(modeldiff)

@testset verbose = true "Angle Diff/Abs action" begin
    @test isapprox(action(model), action(modeldiff), atol = 1e-12)
end

@testset verbose = true "Angle Diff/Abs diff_top_charge" begin
    @test isapprox(diff_top_charge(model), diff_top_charge(modeldiff), atol = 1e-12)
end

