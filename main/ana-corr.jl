# Quantum Rotor
using Revise
import Pkg
Pkg.activate(".")
using LFTSampling
using LFTQuantumRotor
using BDIO
using Statistics
using ADerrors

length(ARGS) == 1 || error("Only one argument is expected! (Path to input file)")
isfile(ARGS[1]) || error("Path provided is not a file")

fname = ARGS[1]
dname = dirname(fname)

ens = LFTSampling.read_ensemble(fname, QuantumRotor)

iT = ens[1].params.iT

nconf = length(ens)

configs = zeros(Float64, nconf, iT)

for i in 1:nconf
    for t in 1:iT
        configs[i,t] = ens[i].phi[t]
    end
    print("$i \r")
end

configs = LFTQuantumRotor.Mod.(configs, 2pi)

Ct = similar(configs, nconf, 25)
meansphi = zeros(iT)

for i in 1:iT
    meansphi[i] = mean(configs[:,i])
end
Ct .= 0.0
for t in 1:25
    println(t)
    for it in 1:iT
        tup = mod1(it+t-1, iT)
        Ct[:,t] .+= (configs[:,it] .- meansphi[it]) .* (configs[:,tup] .- meansphi[tup]) / iT
    end
end

y = Vector{uwreal}(undef, 25)
for i in 1:25
    y[i] = uwreal(Ct[:,i], "test")
end

Ms = similar(y, 24)
for i in 1:length(Ms)
    Ms[i] = log(y[i]/y[i+1])
end

uwerr.(Ms)

open(joinpath(dname,"Meffs.txt"), "w") do file
    for row in eachrow(hcat(0:length(Ms)-1,ADerrors.value.(Ms), ADerrors.err.(Ms)))
        # Join the elements of the row into a string separated by commas
        row_string = join(row, ",")
        # Write the row string to the file, followed by a newline character
        write(file, row_string * "\n")
    end
end
