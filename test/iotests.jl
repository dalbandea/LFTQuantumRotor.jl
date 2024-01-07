
I = 5.0
iT = 10

model = QuantumRotor(I = I, iT = iT, BC = PeriodicBC, disc = StandardDiscretization)

fname = "qriotest.bdio"
isfile(fname) && error("File already exists!")

ens = [deepcopy(model) for i in 1:20]

for i in 1:length(ens)
    randomize!(ens[i])
    save_cnfg(fname, ens[i])
end

rens = LFTSampling.read_ensemble(fname, QuantumRotor)

for i in 1:length(ens)
    @test ens[i].params == rens[i].params
    @test ens[i].phi == rens[i].phi
end

rm(fname, force=true)


fname = "qriotest.bdio"
isfile(fname) && error("File already exists!")

for i in 1:length(ens)
    randomize!(ens[i])
    randomize!(rens[i])
end

save_ensemble(fname, ens)

rens = LFTSampling.read_ensemble(fname, QuantumRotor)

for i in 1:length(ens)
    @test ens[i].params == rens[i].params
    @test ens[i].phi == rens[i].phi
end

rm(fname, force=true)
