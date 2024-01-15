import BDIO: BDIO_write!, BDIO_read
import LFTSampling: save_cnfg_header, read_cnfg_info

function BDIO.BDIO_write!(fb::BDIO.BDIOstream, qrws::LFTQuantumRotor.QuantumRotor) 
    BDIO.BDIO_write!(fb, qrws.phi)
    return nothing
end

function BDIO.BDIO_read(fb::BDIO.BDIOstream, qrws::LFTQuantumRotor.QuantumRotor) 
    BDIO.BDIO_read(fb, qrws.phi)
    return nothing
end

function save_cnfg_header(fb::BDIO.BDIOstream, qrws::QuantumRotor)
    BDIO.BDIO_write!(fb, [qrws.params.I])
    BDIO.BDIO_write!(fb, [convert(Int32, qrws.params.iT)])
    BDIO_write!(fb, string(qrws.params.BC)*"\0")
    BDIO_write!(fb, string(qrws.params.disc)*"\0")
    BDIO_write!(fb, string(typeof(qrws.params.theta))*"\0")
    BDIO.BDIO_write!(fb, [qrws.params.theta])
    return nothing
end

function read_cnfg_info(fname::String, ::Type{LFTQuantumRotor.QuantumRotor}; modul::Module = LFTQuantumRotor)

    fb = BDIO.BDIO_open(fname, "r")

    while BDIO.BDIO_get_uinfo(fb) != 1
        BDIO.BDIO_seek!(fb)
    end

    ifoo    = Vector{Float64}(undef, 1)
    BDIO.BDIO_read(fb, ifoo)
    I    = ifoo[1]
    ifoo    = Vector{Int32}(undef, 1)
    BDIO.BDIO_read(fb, ifoo)
    iT      = convert(Int64, ifoo[1])
    BC      = eval(Meta.parse(BDIO.BDIO_read_str(fb)))
    disc    = eval(Meta.parse(BDIO.BDIO_read_str(fb)))
    thtp    = Base.eval(modul, Meta.parse(BDIO.BDIO_read_str(fb)))
    thfoo   = [zero(thtp)]
    BDIO.BDIO_read(fb, thfoo)
    theta   = thfoo[1]

    model = LFTQuantumRotor.QuantumRotor(
                         Float64, 
                         thtp,
                         I = I, 
                         iT = iT, 
                         BC = BC, 
                         disc = disc,
                         theta = theta
                        )

    return fb, model
end
