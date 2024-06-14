# Quantum Rotor
using Revise
import Pkg
Pkg.activate(".")
using LFTQuantumRotor
using LFTSampling
using Plots
using ADerrors
using BDIO
import Base: write
using Statistics

# cfile = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTQuantumRotor.jl/results/10-qnq0-vs-L/cfgs-run-I3.0-T8-th0.0_D2024-06-10-17-01-09.562/cfgs-run-I3.0-T8-th0.0_D2024-06-10-17-01-09.562.bdio"

# cfile = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTQuantumRotor.jl/results/10-qnq0-vs-L/cfgs-run-I3.0-T20-th0.0_D2024-06-10-18-21-47.835/cfgs-run-I3.0-T20-th0.0_D2024-06-10-18-21-47.835.bdio"

# cfile = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTQuantumRotor.jl/results/10-qnq0-vs-L/cfgs-run-I3.0-T120-th0.0_D2024-06-13-12-17-17.574/cfgs-run-I3.0-T120-th0.0_D2024-06-13-12-17-17.574.bdio"

cfile = ARGS[1]
dircfile = dirname(cfile)
q1q0path = joinpath(dircfile, "q1q0.bdio")
q1q0Q0path = joinpath(dircfile, "q1q0-Q0.bdio")
q1q0Q1path = joinpath(dircfile, "q1q0-Q1.bdio")
q2q0path = joinpath(dircfile, "q2q0.bdio")
q2q0Q0path = joinpath(dircfile, "q2q0-Q0.bdio")
q2q0Q1path = joinpath(dircfile, "q2q0-Q1.bdio")
phi1phi0path = joinpath(dircfile, "phi1phi0.bdio")
phi1phi0Q0path = joinpath(dircfile, "phi1phi0-Q0.bdio")
phi1phi0Q1path = joinpath(dircfile, "phi1phi0-Q1.bdio")
phi2phi0path = joinpath(dircfile, "phi2phi0.bdio")
phi2phi0Q0path = joinpath(dircfile, "phi2phi0-Q0.bdio")
phi2phi0Q1path = joinpath(dircfile, "phi2phi0-Q1.bdio")

qhistorypath = joinpath(dircfile, "Qhistory.png")

ens = LFTSampling.read_ensemble(cfile, QuantumRotor)


function qt(t, model)
    ut = mod1(t+1, model.params.iT)
    return LFTQuantumRotor.Mod(model.phi[ut] - model.phi[t], 2pi)
end

function qtmean(t, ens)
    return mean([qt(t, model) for model in ens])
end

function compute_qtq0(t, model, qtmeans)
    qtq0 = 0.0
    for tt in 1:model.params.iT
        utt = mod1(tt+t,model.params.iT)
        qtq0 += ((qt(tt,model)-qtmeans[tt]) * (qt(utt, model) - qtmeans[utt])) / model.params.iT
    end
    return qtq0
end

function global_qtq0(t, ens, ID; sector = nothing)
    qtmeans = [qtmean(it, ens) for it in 1:ens[1].params.iT]
    qtq0s = [compute_qtq0(t, model, qtmeans) for model in ens]
    if sector == nothing
        uwqtq0 = uwreal(qtq0s, ID)
    else
        Qs = LFTQuantumRotor.top_charge.(ens)
        uwqtq0 = reweight_Q(sector, Qs, qtq0s, ID)
    end
    uwqtq0c = uwqtq0
    return uwqtq0c
end

function compute_phit(t,model)
    return model.phi[t]
end

function phimean(t, ens)
    return mean([compute_phit(t, model) for model in ens])
end

function compute_phitphi0(t, model, phimeans)
    phitphi0 = 0.0
    for tt in 1:model.params.iT
        utt = mod1(tt+t,model.params.iT)
        phitphi0 += ((compute_phit(tt,model)-phimeans[tt]) * (compute_phit(utt, model) - phimeans[utt])) / model.params.iT
    end
    return phitphi0
end

function global_phitphi0(t, ens, ID; sector = nothing)
    phimeans = [phimean(it, ens) for it in 1:ens[1].params.iT]
    phitphi0s = [compute_phitphi0(t, model, phimeans) for model in ens]
    if sector == nothing
        uwphitphi0 = uwreal(phitphi0s, ID)
    else
        Qs = LFTQuantumRotor.top_charge.(ens)
        uwphitphi0 = reweight_Q(sector, Qs, phitphi0s, ID)
    end
    uwphitphi0c = uwphitphi0
    return uwphitphi0c
end


function write(pth::String, uwval::uwreal)
    fb = BDIO.BDIO_open(pth, "d", "Test file")
    ADerrors.write_uwreal(uwval, fb, 8)
    BDIO_close!(fb)
end

function read_uwv(pth::String)
  fb = BDIO_open(pth, "r")
  BDIO_seek!(fb)
  # Read observable
  b = read_uwreal(fb)
  BDIO_close!(fb)
  return b
end


function extract_rwQ(Q_sector, Qs)
    # return [isapprox(Q, Q_sector, atol=1e-5) ? 1.0 : 0.0 for Q in Qs]
    return [isapprox(Q, Q_sector, atol=1e-5) || isapprox(-Q, Q_sector, atol=1e-5) ? 1.0 : 0.0 for Q in Qs]
end

function wQ(Q_sector)
    Qs = LFTQuantumRotor.top_charge.(ens)
    return sum(extract_rwQ(Q_sector, Qs)) / length(Qs)
end


"""
Does reweighting to get value of MC observable `array` at topological sector
`Q`, where the topological charge list is given in `Qs`. Can fail if reweighting
is too aggressive.
"""
function reweight_Q(Q, Qs, array, ID)
    qdeltas = extract_rwQ(Q, Qs)
    uwv = uwreal(array .* qdeltas, ID)
    uwqdeltas = uwreal(qdeltas, ID)
    return uwv / uwqdeltas
end

"""
Just creates a new chain with the values of configurations with topological
charge `Q`. Would only be correct if samples withing a topological sector would
be uncorrelated.
"""
function extract_Q(Q, Qs, array, ID)
    qdeltas = extract_rwQ(Q, Qs)
    arrayQ = filter(!(x -> isapprox(x, 0.0, atol=1e-12)), array .* qdeltas)
    uwarrayQ = uwreal(arrayQ, ID)
    return uwarrayQ
end


Qs = LFTQuantumRotor.top_charge.(ens)
pl = plot(Qs)
savefig(pl, qhistorypath)

uwq1q0c = global_qtq0(1, ens, cfile, sector=nothing)
uwq1q0cQ0 = global_qtq0(1, ens, cfile, sector=0)
uwq1q0cQ1 = global_qtq0(1, ens, cfile, sector=1)

uwq2q0c = global_qtq0(2, ens, cfile, sector=nothing)
uwq2q0cQ0 = global_qtq0(2, ens, cfile, sector=0)
uwq2q0cQ1 = global_qtq0(2, ens, cfile, sector=1)

uwphi1phi0c = global_phitphi0(1, ens, cfile, sector=nothing)
uwphi1phi0cQ0 = global_phitphi0(1, ens, cfile, sector=0)
uwphi1phi0cQ1 = global_phitphi0(1, ens, cfile, sector=1)

uwphi2phi0c = global_phitphi0(2, ens, cfile, sector=nothing)
uwphi2phi0cQ0 = global_phitphi0(2, ens, cfile, sector=0)
uwphi2phi0cQ1 = global_phitphi0(2, ens, cfile, sector=1)

write(q1q0path, uwq1q0c)
write(q1q0Q0path, uwq1q0cQ0)
write(q1q0Q1path, uwq1q0cQ1)

write(q2q0path, uwq2q0c)
write(q2q0Q0path, uwq2q0cQ0)
write(q2q0Q1path, uwq2q0cQ1)

write(phi1phi0path, uwphi1phi0c)
write(phi1phi0Q0path, uwphi1phi0cQ0)
write(phi1phi0Q1path, uwphi1phi0cQ1)

write(phi2phi0path, uwphi2phi0c)
write(phi2phi0Q0path, uwphi2phi0cQ0)
write(phi2phi0Q1path, uwphi2phi0cQ1)

