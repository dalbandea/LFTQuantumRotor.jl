function sample!(qrws::QuantumRotor, mws::AbstractMetropolis; do_winding = false)
    sweep!(qrws, mws)
    do_winding && winding_step!(qrws)
    return nothing
end

function sweep!(qrws::QuantumRotor, mws::AbstractMetropolis)
    for t in 1:qrws.params.iT
        update_variable!(qrws, mws, t)
        # println(t)
    end
    return nothing
end

function update_variable!(qrws::QuantumRotor, mws::AbstractMetropolis, t::Int64)

    #change in phi
    delta = mws.params.weight * randn()

    sini = daction_t(qrws, t)

    qrws.phi[t] += delta

    sfin = daction_t(qrws, t)

    dS = sfin - sini

    accepted = accept_reject!(dS)

    mws.params.nsampled += 1
    if accepted == true
        mws.params.naccepted += 1
    else
        qrws.phi[t] -= delta
    end

    return nothing
end

daction_t(qrws::QuantumRotor, t::Int64) = daction_t(qrws, qrws.params.disc, t)
function daction_t(qrws::QuantumRotor, disc::Type{CPAngleDifferenceDiscretization}, t::Int64)
    ds = zero(qrws.PRC)
    if t < qrws.params.iT
        ds = qrws.params.I/2 * (action_t(qrws, t) + boundary_action(qrws))
    end
    return ds
end

function winding_step!(qrws::QuantumRotor)
    ws_cp = deepcopy(qrws)

    sini = action(qrws)

    r = rand()

    if r > 0.5
        winding!(qrws)
    else
        antiwinding!(qrws)
    end

    sfin = action(qrws)

    ds = sfin - sini

    metropolis_accept_reject!(qrws, ws_cp, ds)

    return nothing
end
