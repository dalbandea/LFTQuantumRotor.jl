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
using ProgressBars
using DelimitedFiles

# cfile = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTQuantumRotor.jl/results/10-qnq0-vs-L/cfgs-run-I3.0-T8-th0.0_D2024-06-10-17-01-09.562/cfgs-run-I3.0-T8-th0.0_D2024-06-10-17-01-09.562.bdio"

# cfile = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTQuantumRotor.jl/results/10-qnq0-vs-L/cfgs-run-I3.0-T20-th0.0_D2024-06-10-18-21-47.835/cfgs-run-I3.0-T20-th0.0_D2024-06-10-18-21-47.835.bdio"

# cfile = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTQuantumRotor.jl/results/10-qnq0-vs-L/cfgs-run-I3.0-T120-th0.0_D2024-06-13-12-17-17.574/cfgs-run-I3.0-T120-th0.0_D2024-06-13-12-17-17.574.bdio"

# cfile = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTQuantumRotor.jl/results/10-qnq0-vs-L/500k-statistics/cfgs-run-I3.0-T20-th0.0_D2024-06-10-18-21-47.835/cfgs-run-I3.0-T20-th0.0_D2024-06-10-18-21-47.835.bdio"

# cfile = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTQuantumRotor.jl/results/10-qnq0-vs-L/500k-statistics/cfgs-run-I3.0-T250-th0.0_D2024-06-13-21-27-24.383/cfgs-run-I3.0-T250-th0.0_D2024-06-13-21-27-24.383.bdio"


# cfile = ARGS[1]
# dircfile = dirname(cfile)
# qtq0path = joinpath(dircfile, "Rqtq0.bdio")
# qtq0Q0path = joinpath(dircfile, "Rqtq0-Q0.bdio")
# qtq0Q1path = joinpath(dircfile, "Rqtq0-Q1.bdio")
# phitphi0path = joinpath(dircfile, "Rphitphi0.bdio")
# phitphi0Q0path = joinpath(dircfile, "Rphitphi0-Q0.bdio")
# phitphi0Q1path = joinpath(dircfile, "Rphitphi0-Q1.bdio")

# qhistorypath = joinpath(dircfile, "Qhistory.png")

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
    tt = 1
    utt = mod1(tt+t,model.params.iT)
    qtq0 += ((qt(tt,model)-qtmeans[tt]) * (qt(utt, model) - qtmeans[utt]))
    return qtq0
end

# function compute_qtq0(t, model, qtmeans)
#     qtq0 = 0.0
#     for tt in 1:model.params.iT
#         utt = mod1(tt+t,model.params.iT)
#         qtq0 += ((qt(tt,model)-qtmeans[tt]) * (qt(utt, model) - qtmeans[utt])) / model.params.iT
#     end
#     return qtq0
# end

function global_qtq0(t, ens, qtmeans, ID; sector = nothing)
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
    tt = 1
    utt = mod1(tt+t,model.params.iT)
    phitphi0 += ((compute_phit(tt,model)-phimeans[tt]) * (compute_phit(utt, model) - phimeans[utt]))
    return phitphi0
end

# function compute_phitphi0(t, model, phimeans)
#     phitphi0 = 0.0
#     for tt in 1:model.params.iT
#         utt = mod1(tt+t,model.params.iT)
#         phitphi0 += ((compute_phit(tt,model)-phimeans[tt]) * (compute_phit(utt, model) - phimeans[utt])) / model.params.iT
#     end
#     return phitphi0
# end

function global_phitphi0(t, ens, phimeans, ID; sector = nothing)
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

function f_R(R, uwsqtq0)
    T = length(uwsqtq0)
    R <= T/2 || error("R too big")
    res = uwsqtq0[1]
    for t in 1:R
        res += uwsqtq0[t+1] + uwsqtq0[T-t+1]
    end
    if R == T/2
        res -= uwsqtq0[R+1]
    end
    return res
end


qtmeans = [qtmean(it, ens) for it in ProgressBar(1:ens[1].params.iT)]

phimeans = [phimean(it, ens) for it in 1:ens[1].params.iT]

uwqtq0c = [global_qtq0(t, ens, qtmeans, cfile, sector=nothing) for t in ProgressBar(0:ens[1].params.iT-1)]

uwqtq0cQ0 = [global_qtq0(t, ens, qtmeans, cfile, sector=0) for t in ProgressBar(0:ens[1].params.iT-1)]

uwqtq0cQ1 = [global_qtq0(t, ens, qtmeans, cfile, sector=1) for t in ProgressBar(0:ens[1].params.iT-1)]

uwphitphi0c = [global_phitphi0(t, ens, phimeans, cfile, sector=nothing) for t in ProgressBar(0:ens[1].params.iT-1)]

uwphitphi0cQ0 = [global_phitphi0(t, ens, phimeans, cfile, sector=0) for t in 0:ens[1].params.iT-1]

uwphitphi0cQ1 = [global_phitphi0(t, ens, phimeans, cfile, sector=1) for t in 0:ens[1].params.iT-1]

Rs = collect(0:Int64(ens[1].params.iT/2))

fRsqtq0 = [f_R(R, uwqtq0c) for R in Rs]

fRsqtq0Q0 = [f_R(R, uwqtq0cQ0) for R in Rs]

fRsqtq0Q1 = [f_R(R, uwqtq0cQ1) for R in Rs]

fRsphitphi0 = [f_R(R, uwphitphi0c) for R in Rs]
fRsphitphi0Q0 = [f_R(R, uwphitphi0cQ0) for R in Rs]
fRsphitphi0Q1 = [f_R(R, uwphitphi0cQ1) for R in Rs]


uwerr.(fRsqtq0)
uwerr.(fRsqtq0Q0)
uwerr.(fRsqtq0Q1)

uwerr.(fRsphitphi0)
uwerr.(fRsphitphi0Q0)
uwerr.(fRsphitphi0Q1)

# write(qtq0path, uwqtq0c)
# write(qtq0Q0path, uwqtq0cQ0)
# write(qtq0Q1path, uwqtq0cQ1)

# write(phitphi0path, uwphitphi0c)
# write(phitphi0Q0path, uwphitphi0cQ0)
# write(phitphi0Q1path, uwphitphi0cQ1)

writedlm("Rqtq0.txt", hcat(Rs, ADerrors.value.(fRsqtq0), ADerrors.err.(fRsqtq0)), ',')

writedlm("Rqtq0Q0.txt", hcat(Rs, ADerrors.value.(fRsqtq0Q0), ADerrors.err.(fRsqtq0Q0)), ',')
writedlm("Rqtq0Q1.txt", hcat(Rs, ADerrors.value.(fRsqtq0Q1), ADerrors.err.(fRsqtq0Q1)), ',')


writedlm("Rphitphi0.txt", hcat(Rs, ADerrors.value.(fRsphitphi0), ADerrors.err.(fRsphitphi0)), ',')
writedlm("Rphitphi0Q0.txt", hcat(Rs, ADerrors.value.(fRsphitphi0Q0), ADerrors.err.(fRsphitphi0Q0)), ',')
writedlm("Rphitphi0Q1.txt", hcat(Rs, ADerrors.value.(fRsphitphi0Q1), ADerrors.err.(fRsphitphi0Q1)), ',')

plot(ADerrors.value.(fRsqtq0Q0), yerr=ADerrors.err.(fRsqtq0Q0))
