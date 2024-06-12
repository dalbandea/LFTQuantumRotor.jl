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

# cfile = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTQuantumRotor.jl/results/10-qnq0-vs-L/cfgs-run-I3.0-T8-th0.0_D2024-06-10-17-01-09.562/cfgs-run-I3.0-T8-th0.0_D2024-06-10-17-01-09.562.bdio"
# cfile = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTQuantumRotor.jl/results/10-qnq0-vs-L/cfgs-run-I3.0-T20-th0.0_D2024-06-10-18-21-47.835/cfgs-run-I3.0-T20-th0.0_D2024-06-10-18-21-47.835.bdio"
cfile = ARGS[1]
dircfile = dirname(cfile)
q1q0path = joinpath(dircfile, "q1q0.bdio")
q1q0Q0path = joinpath(dircfile, "q1q0-Q0.bdio")
q1q0Q1path = joinpath(dircfile, "q1q0-Q1.bdio")
qhistorypath = joinpath(dircfile, "Qhistory.pdf")

ens = LFTSampling.read_ensemble(cfile, QuantumRotor)

function compute_q1q0(model)
    q0 = LFTQuantumRotor.Mod(model.phi[2] - model.phi[1], 2pi)
    q1 = LFTQuantumRotor.Mod(model.phi[3] - model.phi[2], 2pi)
    return q1*q0
end

function compute_qt(t,model)
    return LFTQuantumRotor.Mod(model.phi[t+1] - model.phi[t], 2pi)
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

function global_q1q0(ens, ID)
    q1q0s = compute_q1q0.(ens)
    q0s = [compute_qt(1, model) for model in ens]
    q1s = [compute_qt(2, model) for model in ens]
    uwq1q0 = uwreal(q1q0s, cfile)
    uwq0 = uwreal(q0s, cfile)
    uwq1 = uwreal(q1s, cfile)
    uwq1q0c = uwq1q0 - uwq0*uwq1
    return uwq1q0c
end

function extract_rwQ(Q_sector, Qs)
    # return [isapprox(Q, Q_sector, atol=1e-5) ? 1.0 : 0.0 for Q in Qs]
    return [isapprox(Q, Q_sector, atol=1e-5) || isapprox(-Q, Q_sector, atol=1e-5) ? 1.0 : 0.0 for Q in Qs]
end

function q1q0_at_Q(Q, Qs, ens, ID)
    q1q0s = compute_q1q0.(ens)
    q0s = [compute_qt(1, model) for model in ens]
    q1s = [compute_qt(2, model) for model in ens]
    uwq1q0 = reweight_Q(Q, Qs, q1q0s, ID)
    uwq0 = reweight_Q(Q, Qs, q0s, ID)
    uwq1 = reweight_Q(Q, Qs, q1s, ID)
    uwq1q0c = uwq1q0 - uwq0*uwq1
    # uwq1q0 = extract_Q(Q, Qs, q1q0s, ID)
    # uwq0 = extract_Q(Q, Qs, q0s, ID)
    # uwq1 = extract_Q(Q, Qs, q1s, ID)
    # uwq1q0c = uwq1q0 - uwq0*uwq1
    return uwq1q0c
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

uwq1q0c = global_q1q0(ens, cfile)
uwq1q0cQ0 = q1q0_at_Q(0, Qs, ens, cfile)
uwq1q0cQ1 = q1q0_at_Q(1, Qs, ens, cfile)

# uwq1q0c = global_q1q0(ens, cfile*"t")
# uwq1q0cQ0 = q1q0_at_Q(0, Qs, ens, cfile*"--Q0")
# uwq1q0cQ1 = q1q0_at_Q(1, Qs, ens, cfile*"--Q1")

write(q1q0path, uwq1q0c)
write(q1q0Q0path, uwq1q0cQ0)
write(q1q0Q1path, uwq1q0cQ1)



# Trash

# Qs = LFTQuantumRotor.top_charge.(ens)
# pl = plot(Qs)

# uwq1q0c = global_q1q0(ens, cfile)
# uwq1q0cQ0 = q1q0_at_Q(0, Qs, ens, cfile)
# uwq1q0cQ1 = q1q0_at_Q(1, Qs, ens, cfile)
# uwq1q0cQm1 = q1q0_at_Q(-1, Qs, ens, cfile)

# uwq1q0cQ2 = q1q0_at_Q(2, Qs, ens, cfile)
# uwq1q0cQm2 = q1q0_at_Q(-2, Qs, ens, cfile)

# deltaq0 = extract_rwQ(0, Qs)
# deltaq1 = extract_rwQ(1, Qs)
# deltaqm1 = extract_rwQ(-1, Qs)

# q1q0Q0 = filter(!(x -> isapprox(x, 0.0, atol=1e-5)), compute_q1q0.(ens) .* deltaq0)

# q1q0Q1 = filter(!(x -> isapprox(x, 0.0, atol=1e-5)), compute_q1q0.(ens) .* deltaq1)
# q1q0Qm1 = filter(!(x -> isapprox(x, 0.0, atol=1e-5)), compute_q1q0.(ens) .* deltaqm1)

# uwq1q0Q0 = uwreal(q1q0Q0, cfile*"0")
# uwq1q0Q1 = uwreal(q1q0Q1, cfile*"1")
# uwq1q0Qm1 = uwreal(q1q0Qm1, cfile*"m1")

# uwq1q0c

# uwq1q0Q0

# uwq1q0cQ0

# uwq1q0Q1

# # Compute with right weights

# wz = 2* sum(deltaq1) +  sum(deltaq0)

# w1 = 2*sum(deltaq1)/wz
# w0 = sum(deltaq0)/wz

# w0 * uwq1q0Q0 + w1*uwq1q0Q1

# q1s = [compute_qt(2, model) for model in ens]

# uwq1 = reweight_Q(-1, Qs, q1s, "hi")

