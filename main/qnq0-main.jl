# Quantum Rotor
using Revise
import Pkg
Pkg.activate(".")
using LFTQuantumRotor
using LFTSampling
using TOML
using Dates
using ProgressBars
using Logging

length(ARGS) == 1 || error("Only one argument is expected! (Path to input file)")
isfile(ARGS[1]) || error("Path provided is not a file")

infile = ARGS[1]
# infile = "main/infile.in"
pdata = TOML.parsefile(infile)

# Read model parameters

I = pdata["Model params"]["I"]
iT = pdata["Model params"]["T"]
BC = eval(Meta.parse(pdata["Model params"]["BC"]))
disc = eval(Meta.parse(pdata["Model params"]["disc"]))
theta = pdata["Model params"]["theta"]

# Read HMC parameters

tau = pdata["HMC params"]["tau"]
nsteps = pdata["HMC params"]["nsteps"]
ntherm = pdata["HMC params"]["ntherm"]
ntraj = pdata["HMC params"]["ntraj"]
discard = pdata["HMC params"]["discard"]
windings = pdata["HMC params"]["windings"]

# Working directory

wdir = pdata["Working directory"]["wdir"]

dt = Dates.now()
wdir_sufix = "_D"*Dates.format(dt, "yyyy-mm-dd-HH-MM-SS.ss")
fname = "cfgs-run-I$I-T$iT-th$theta"*wdir_sufix

fdir = joinpath(wdir, fname)
configfile = joinpath(fdir, fname*".bdio")
mkpath(fdir)
cp(infile, joinpath(fdir,splitpath(infile)[end]))


model = QuantumRotor(
                     I = I, 
                     iT = iT, 
                     BC = BC, 
                     disc = disc, 
                     theta = theta
                    )

# randomize!(model)

LFTQuantumRotor.coldstart!(model)

smplr = HMC(integrator = Leapfrog(tau, nsteps))
samplerws = LFTSampling.sampler(model, smplr)

# Logging.disable_logging(Logging.Debug)
# Logging.disable_logging(Logging.Info)

# Thermalization

for i in 1:ntherm
    @time sample!(model, samplerws, do_winding=windings)
    model.phi .= LFTQuantumRotor.Mod.(model.phi, 2pi)
end


# Run

ens = [deepcopy(model) for i in 1:ntraj]

LFTQuantumRotor.coldstart!.(ens)

@time for i in ProgressBar(1:ntraj)
    for j in 1:discard
        sample!(model, samplerws, do_winding=windings);
    end
    sample!(model, samplerws, do_winding=windings);
    model.phi .= LFTQuantumRotor.Mod.(model.phi, 2pi)
    ens[i].phi .= model.phi
    # save_cnfg(configfile, model)
end

save_ensemble(configfile, ens)

