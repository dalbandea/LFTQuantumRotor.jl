# Quantum Rotor
using Revise
import Pkg
Pkg.activate(".")
using LFTQuantumRotor
using LFTSampling

I = 1.0
iT = 10

model = QuantumRotor(I = I, iT = iT, BC = PeriodicBC, disc = CPAngleDifferenceDiscretization)

randomize!(model)

action(model)

top_charge(model)


function random_OBC_CP(I)
    xi = rand()
    return sqrt(2/I) * erfinv((2*xi-1) * erf(sqrt(I/2)*pi))
end

function PDF_OBC_CP(I)
    Z = sqrt(2pi/I) * erf(sqrt(I/2)*pi)
    return x -> exp(-I/2 * x^2) / Z
end

function logprob_OBC_CP(I)
    Z = sqrt(2pi/I) * erf(sqrt(I/2)*pi)
    return x -> -I/2 * x^2 - log(Z)
end

function PDF_configuration_OBC_CP(phi, I, T)
    factorized_PDF = PDF_OBC_CP(I)
    res = 1.0
    for i in 1:T
        res = res * factorized_PDF(phi[i])
    end
    return res
end

function logprob_configuration_OBC_CP(phi, I, T)
    factorized_logprob = logprob_OBC_CP(I)
    res = 0.0
    for i in 1:T
        res += factorized_logprob(phi[i])
    end
    return res
end



# Backup configuration
model_cp = deepcopy(model)

# Initial action
logp_i = -action(model)
logq_i = logprob_configuration_OBC_CP(model.phi, model.params.I, model.params.iT-1)

# Generate new configuration
model.phi .= [random_OBC_CP(1.0) for _ in 1:model.params.iT]

# Final action
logp_f = -action(model)
logq_f = logprob_configuration_OBC_CP(model.phi, model.params.I, model.params.iT-1)

# Accept-reject step
pacc = exp(logp_f - logp_i + logq_i - logq_f)



# Periodic Boundary Conditions
I = 4.0
iT = 400
model = QuantumRotor(I = I, iT = iT, BC = PeriodicBC, disc = CPAngleDifferenceDiscretization)
randomize!(model)

for _ in 1:1000
    # Backup configuration
    model_cp = deepcopy(model)
    # Initial action
    logp_i = -action(model)
    logq_i = logprob_configuration_OBC_CP(model.phi, model.params.I, model.params.iT-1)
    # Generate new configuration
    model.phi .= [random_OBC_CP(model.params.I) for _ in 1:model.params.iT]
    # Final action
    logp_f = -action(model)
    logq_f = logprob_configuration_OBC_CP(model.phi, model.params.I, model.params.iT-1)
    # Accept-reject step
    pacc = exp(logp_f - logp_i + logq_i - logq_f)
    println(pacc)
end

model.phi .= [random_OBC_CP(model.params.I) for _ in 1:model.params.iT]
top_charge(model)


# Ideal code: Periodic Boundary Conditions

I = 100.0
iT = 1000
model = QuantumRotor(I = I, iT = iT, BC = PeriodicBC, disc = CPAngleDifferenceDiscretization)
randomize!(model)
dist = QuantumRotorOBCDistribution(model.params.I)
smplr = MetropolisHastings(dist = dist)
samplerws = LFTsModule.sampler(model, smplr)

Qs = Vector{Float64}()
Ss = Vector{Float64}()

sample!(model, samplerws)

for _ in 1:10000
    sample!(model, samplerws)
    # push!(Qs, top_charge(model))
    # push!(Ss, action(model))
end

LFTsModule.acceptance(samplerws)

sample!(model, samplerws)

using ADerrors

uwchi = uwreal(Qs.^2/100, "test")

plot(Qs)

logq_i = log_prob(model, samplerws)
logp_i = -action(model)
model.phi[8] = sample(dist) 
logq_f = log_prob(model, samplerws)
logp_f = -action(model)
exp(logp_f - logp_i + logq_i - logq_f)

LFTsModule.update!(model, samplerws)

randomize!(model)


model.phi .= [sample(dist) for i in 1:1000]

sample(QuantumRotorOBCDistribution(10000000), 10)


# Normal Metropolis


I = 4.0
iT = 400
model = QuantumRotor(I = I, iT = iT, BC = PeriodicBC, disc = CPAngleDifferenceDiscretization)

randomize!(model)
top_charge(model)

smplr = Metropolis(weight=0.01)

samplerws = LFTsModule.sampler(model, smplr)

Qs = Vector{Float64}()
Ss = Vector{Float64}()

for i in 1:1000
    println(i)
    sample!(model, samplerws, do_winding = true)
    push!(Qs, top_charge(model))
    push!(Ss, action(model))
end
top_charge(model)

plot(Qs)


using ADerrors

uwchi = uwreal(Qs.^2/1000, "test4")
uwerr(uwchi)
uwchi


# Standard action, Open Boundary Conditions

I = 10.0
iT = 10000
model = QuantumRotor(I = I, iT = iT, BC = OpenBC, disc = StAngleDifferenceDiscretization)
model.phi .= sample(dist, iT) # gets stuck if initialized randomly
dist = QuantumRotorOBCDistribution(model.params.I)
smplr = MetropolisHastings(dist = dist)
samplerws = LFTsModule.sampler(model, smplr)
sample!(model, samplerws)

for _ in 1:10000
    sample!(model, samplerws)
    # push!(Qs, top_charge(model))
    # push!(Ss, action(model))
end

LFTsModule.acceptance(samplerws)


# Standard action, Periodic Boundary Conditions

I = 1000.0
iT = 1000000
model = QuantumRotor(I = I, iT = iT, BC = PeriodicBC, disc = StAngleDifferenceDiscretization)
model.phi .= sample(dist, iT) # gets stuck if initialized randomly
dist = QuantumRotorOBCDistribution(model.params.I)
smplr = MetropolisHastings(dist = dist)
samplerws = LFTsModule.sampler(model, smplr)
sample!(model, samplerws)

for _ in 1:1000
    sample!(model, samplerws)
    # push!(Qs, top_charge(model))
    # push!(Ss, action(model))
end

LFTsModule.acceptance(samplerws)


# Standard action OBC: HMC

I = 100
iT = 2
model = QuantumRotor(I = I, iT = iT, BC = OpenBC, disc = StAngleDifferenceDiscretization)
randomize!(model)

smplr = HMC(integrator = Leapfrog(1.0, 10))
samplerws = LFTSampling.sampler(model, smplr)

sample!(model, samplerws)

Ss = Vector{Float64}()
for i in 1:100000
    sample!(model, samplerws)
    S = action(model)
    # push!(Ss, S)
    push!(Ss, LFTsModule.Mod(model.phi[1],2pi)^2)
end

using ADerrors

id = "test19"

uws = uwreal(Ss, id)
uwerr(uws)
uws
taui(uws,id)
dtaui(uws,id)

cosphi = 1-uws/(I*(iT-1))
uwerr(cosphi)
cosphi
tau = taui(cosphi, id)
dtau = dtaui(cosphi, id)

using Plots

r = rho(uws,id)[1:200]
dr = drho(uws,id)[1:200]
plot(r, yerr=dr, seriestype=:scatter)


Is = Vector{Float64}()
taus = Vector{Float64}()
dtaus = Vector{Float64}()

push!(Is, I)
push!(taus, tau)
push!(dtaus, dtau)

plot(Is[4:end], taus[4:end], yerr=dtaus[4:end], seriestype=:scatter , yaxis=:log)


#

using ADerrors

Is = Vector{Float64}()
taus = Vector{Float64}()
dtaus = Vector{Float64}()

for I in 1:0.5:10
    iT = convert(Int64,I*100+1)
    # iT = 2
    model = QuantumRotor(I = I, iT = iT, BC = OpenBC, disc = StAngleDifferenceDiscretization)
    randomize!(model)
    smplr = HMC(integrator = Leapfrog(1.0, 30))
    samplerws = LFTsModule.sampler(model, smplr)
    Ss = Vector{Float64}()
    for i in 1:10000
        sample!(model, samplerws)
    end
    for i in 1:100000
        sample!(model, samplerws)
        S = action(model)
        push!(Ss, S)
    end
    id = "test15"
    uws = uwreal(Ss, id)
    uwerr(uws)
    uws
    tau = taui(uws, id)
    dtau = dtaui(uws, id)
    push!(Is, I)
    push!(taus, tau)
    push!(dtaus, dtau)
end

using Plots

plot(Is, taus, yerr=dtaus, seriestype=:scatter)
plot!(xlabel=L"I", ylabel=L"\tau_{S}", title="Quantum Rotor")



# Standard action PBC: HMC

I = 10
iT = I*100+1
model = QuantumRotor(I = I, iT = iT, BC = PeriodicBC, disc = StAngleDifferenceDiscretization)
# randomize!(model)
model.phi .= 0.0

smplr = HMC(integrator = Leapfrog(1.0, 75))

samplerws = LFTSampling.sampler(model, smplr)

sample!(model, samplerws, do_winding=false)

Ss = Vector{Float64}()
Qs = Vector{Float64}()
for i in 1:1000000
    sample!(model, samplerws, do_winding=false)
    S = action(model)
    Q = top_charge(model)
    push!(Ss, S)
    push!(Qs, Q)
    # push!(Ss, LFTsModule.Mod(model.phi[1],2pi)^2)
end

using ADerrors

id = "test41"

uws = uwreal(Ss, id)
uwerr(uws)
uws
taui(uws,id)
dtaui(uws,id)

cosphi = 1-uws/(I*(iT))
uwerr(cosphi)
cosphi

tau = taui(cosphi, id)
dtau = dtaui(cosphi, id)


uwchi = uwreal(Qs.^2/(iT), id)
uwerr(uwchi)
uwchi

uwchi = uwreal(Qs.^1, id)
uwerr(uwchi)
uwchi

taui(uwchi, id)
dtaui(uwchi, id)

using Plots

r = rho(uws,id)[1:200]
dr = drho(uws,id)[1:200]
plot(r, yerr=dr, seriestype=:scatter)


Is = Vector{Float64}()
taus = Vector{Float64}()
dtaus = Vector{Float64}()

push!(Is, I)
push!(taus, tau)
push!(dtaus, dtau)

plot(Is[4:end], taus[4:end], yerr=dtaus[4:end], seriestype=:scatter , yaxis=:log)


#

using ADerrors

Is = Vector{Float64}()
taus_S = Vector{Float64}()
dtaus_S = Vector{Float64}()
taus_chi = Vector{Float64}()
dtaus_chi = Vector{Float64}()
uwss = Vector{uwreal}()
uwchis = Vector{uwreal}()

for I in 1:10
    iT = I*100+1
    # iT = 10+1
    model = QuantumRotor(I = I, iT = iT, BC = PeriodicBC, disc = StAngleDifferenceDiscretization)
    randomize!(model)
    smplr = HMC(integrator = Leapfrog(1.0, 100))
    samplerws = LFTsModule.sampler(model, smplr)
    Ss = Vector{Float64}()
    Qs = Vector{Float64}()
    for i in 1:10000
        sample!(model, samplerws, do_winding=false)
    end
    for i in 1:100000
        sample!(model, samplerws, do_winding=false)
        S = action(model)
        Q = top_charge(model)
        push!(Ss, S)
        push!(Qs, Q)
    end
    id = "test7"
    uws = uwreal(Ss, id)
    uwerr(uws)
    uwchi = uwreal(Qs.^2/iT, id)
    uwerr(uwchi)
    push!(uwss, uws)
    push!(uwchis, uwchi)
    push!(Is, I)
    tau = taui(uws, id)
    dtau = dtaui(uws, id)
    push!(taus_S, tau)
    push!(dtaus_S, dtau)
    tau = taui(uwchi, id)
    dtau = dtaui(uwchi, id)
    push!(taus_chi, tau)
    push!(dtaus_chi, dtau)
end

using Plots


plot(Is, taus_S, yerr=dtaus_S, seriestype=:scatter)
plot!(xlabel=L"I", ylabel=L"\tau_{S}", title="Quantum Rotor, PBC, winding")

plot(Is, taus_chi, yerr=dtaus_chi, seriestype=:scatter)
plot!(xlabel=L"I", ylabel=L"\tau_{\chi_Q}", title="Quantum Rotor")

cosphis = 1 .- uwss ./ (Is .* (Is*100 .+ 1))

uwerr.(cosphis)
cosphis



# Write to file

using ADerrors


function main()
    for I in 7.5:0.5:10.0
        iT = convert(Int64, I*100+1)
        # iT = 10+1
        model = QuantumRotor(I = I, iT = iT, BC = PeriodicBC, disc = StAngleDifferenceDiscretization)
        # randomize!(model)
        model.phi .= 0.0
        smplr = HMC(integrator = Leapfrog(1.0, 75))
        samplerws = LFTsModule.sampler(model, smplr)
        Ss = Vector{Float64}()
        Qs = Vector{Float64}()
        datafile = "pbc-I$I-T$iT.txt"
        for i in 1:10000
            sample!(model, samplerws, do_winding=false)
        end
        for i in 1:1000000
            sample!(model, samplerws, do_winding=false)
            S = action(model)
            Q = top_charge(model)
            io_stat = open(datafile, "a")
            print(io_stat, S, ",", Q, "\n")
            close(io_stat)
        end
    end
end

main()


using Plots, DelimitedFiles

filterdir(dir::String, text::String) = filterdir(dir, [text])

function filterdir(dir::String, texts::Array{String})
    dirfiles = readdir(dir)
    occurrences = filter(s -> occursin(Regex("$(texts[1])"), s), dirfiles)
    for text in texts
        occurrences = filter(s -> occursin(Regex("$text"), s), occurrences)
    end
    return joinpath.(dir, occurrences)
end

datafiles = sort(filterdir(".", "pbc"), lt=natural)

id = "test53"

uwchis = Vector{uwreal}()
uwss = Vector{uwreal}()
for datafile in datafiles
    println("Reading $datafile...")
    data = readdlm(datafile, ',')
    uws = uwreal(data[:,1], datafile)
    uwerr(uws)
    # uwchi = uwreal(data[:,2], datafile)
    # uwerr(uwchi)
    push!(uwss, uws)
    # push!(uwchis, uwchi)
end


# tauchis = Vector{Float64}(undef, length(uwchis))
# dtauchis = Vector{Float64}(undef, length(uwchis))
tauss = Vector{Float64}(undef, length(uwss))
dtauss = Vector{Float64}(undef, length(uwss))
for i in eachindex(tauss)
    # tauchis[i] = taui(uwchis[i], datafiles[i]) 
    # dtauchis[i] = dtaui(uwchis[i], datafiles[i]) 
    tauss[i] = taui(uwss[i], datafiles[i]) 
    dtauss[i] = dtaui(uwss[i], datafiles[i]) 
end

Is = collect(1:0.5:7)

pl = plot(Is, tauss, yerr=dtauss, seriestype=:scatter)
plot!(pl, xlabel=L"I", ylabel=L"\tau_{S}", title="Quantum Rotor, PBC")

pl = plot(Is, tauchis, yerr=dtauchis, seriestype=:scatter, yscale=:log10)
plot!(pl, xlabel=L"I", ylabel=L"\tau_{Q}", title="Quantum Rotor, PBC")

id = "test56"
uws = uwreal(Ss, id)
uwerr(uws)
taui(uws, id)
dtaui(uws, id)

plot(Is, taus_S, yerr=dtaus_S, seriestype=:scatter)
plot!(xlabel=L"I", ylabel=L"\tau_{S}", title="Quantum Rotor, PBC, winding")

plot(Is, taus_chi, yerr=dtaus_chi, seriestype=:scatter)
plot!(xlabel=L"I", ylabel=L"\tau_{\chi_Q}", title="Quantum Rotor")

cosphis = 1 .- uwss ./ (Is .* (Is*100 .+ 1))

uwerr.(cosphis)
cosphis



# Trivializing Map Standard action PBC: HMC

I = 3.0
iT = I*100+1
model = QuantumRotor(I = I, iT = iT, BC = PeriodicBC, disc = StAngleDifferenceTrivMapDiscretization)
# randomize!(model)
model.phi .= 0.0

smplr = HMC(integrator = Leapfrog(6.75, 25))
samplerws = LFTsModule.sampler(model, smplr)

sample!(model, samplerws, do_winding=false)

Ss = Vector{Float64}()
Qs = Vector{Float64}()
for i in 1:1000000
    sample!(model, samplerws, do_winding=false)
    S = action(model)
    model.phi .= model.phi ./ sqrt(10*model.params.I)
    Q = top_charge(model)
    model.phi .= model.phi .* sqrt(10*model.params.I)
    push!(Ss, S)
    push!(Qs, Q)
    # push!(Ss, LFTsModule.Mod(model.phi[1],2pi)^2)
end

acc = Ss .!= circshift(Ss, -1)
mean(acc)

using ADerrors

id = "test43"

uws = uwreal(Ss, id)
uwerr(uws)
uws
taui(uws,id)
dtaui(uws,id)

cosphi = 1-uws/(I*(iT))
uwerr(cosphi)
cosphi

tau = taui(cosphi, id)
dtau = dtaui(cosphi, id)


uwchi = uwreal(Qs.^2/(iT), id)
uwerr(uwchi)
uwchi

uwchi = uwreal(Qs.^1, id)
uwerr(uwchi)
uwchi

taui(uwchi, id)
dtaui(uwchi, id)

using Plots

r = rho(uws,id)[1:200]
dr = drho(uws,id)[1:200]
plot(r, yerr=dr, seriestype=:scatter)


Is = Vector{Float64}()
taus = Vector{Float64}()
dtaus = Vector{Float64}()

push!(Is, I)
push!(taus, tau)
push!(dtaus, dtau)

plot(Is[4:end], taus[4:end], yerr=dtaus[4:end], seriestype=:scatter , yaxis=:log)

