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

# cfile = "/home/david/git/dalbandea/phd/codes/6-LFTs/LFTModels/LFTQuantumRotor.jl/results/11-qnq0-vs-L-OBC/cfgs-run-I3.0-T100-th0.0_D2024-06-15-20-10-41.733/cfgs-run-I3.0-T100-th0.0_D2024-06-15-20-10-41.733.bdio"

cfile = ARGS[1]

dircfile = dirname(cfile)
phi1phi0path = joinpath(dircfile, "phi1phi0.bdio")
phi2phi0path = joinpath(dircfile, "phi2phi0.bdio")

ens = LFTSampling.read_ensemble(cfile, QuantumRotor)


function compute_phit(t,model)
    return model.phi[t]
end

function phimean(t, ens)
    return mean([compute_phit(t, model) for model in ens])
end

# function compute_phitphi0(t, model, phimeans)
#     phitphi0 = 0.0
#     for tt in 1:model.params.iT
#         utt = mod1(tt+t,model.params.iT)
#         phitphi0 += ((compute_phit(tt,model)-phimeans[tt]) * (compute_phit(utt, model) - phimeans[utt])) / model.params.iT
#     end
#     return phitphi0
# end

function compute_phitphi0(t, t0, model, phimeans)
    phitphi0 = 0.0
    ut = t0+t
    phitphi0 += ((compute_phit(t0,model)-phimeans[t0]) * (compute_phit(ut, model) - phimeans[ut]))
    return phitphi0
end


function global_phitphi0(t, t0, ens, ID)
    phimeans = [phimean(it, ens) for it in 1:ens[1].params.iT]
    phitphi0s = [compute_phitphi0(t, t0, model, phimeans) for model in ens]
    uwphitphi0 = uwreal(phitphi0s, ID)
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


t0 = convert(Int64, ens[1].params.iT / 2)

uwphi1phi0c = global_phitphi0(1, t0, ens, cfile)
uwphi2phi0c = global_phitphi0(2, t0, ens, cfile)

uwerr(uwphi1phi0c)
uwerr(uwphi2phi0c)

write(phi1phi0path, uwphi1phi0c)
write(phi2phi0path, uwphi2phi0c)


