# Depends on EvGen.jl, LinAlg.jl and Numeric.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Mechanism to apply a cut to generated events"
module EvCut

using ..EvGen: electron_momentum, Event, OUTGOING_COUNT, outgoing_momenta,
               min_photon_energy
using ..LinAlg: X, Z, E
using ..Numeric: Float
using LinearAlgebra: dot

export EventCut, keep_event


"Cuts on generated events"
struct EventCut
    "Cut on maximum cosine of (beam, photons) angle"
    a_cut::Float
    
    "Cut on maximum cosine of (photon, photon) angle"
    b_cut::Float
    
    "Cut on minimum photon energy"
    e_min::Float
    
    "Cut on minimum cosine of (beam, normal to the photon plane) angle"
    sin_cut::Float
end


"Decide whether a generated event passes the cut or should be rejected"
function keep_event(cut::EventCut, event::Event)::Bool
    # Check if the outgoing photons pass the energy cut
    if min_photon_energy(event) < cut.e_min
        return false
    end

    # Get the incoming electron 4-momentum and outgoing photon 4-momenta
    p_el = electron_momentum(event)
    ps_out = outgoing_momenta(event)

    # Check if the (beam, photon) angles pass the cut
    cos_nums = ps_out[:, X:Z] * p_el[X:Z]
    cos_denoms = ps_out[E] * p_el[E]
    for (num, denom) ∈ zip(cos_nums, cos_denoms)
        if abs(num) > cut.a_cut * denom
            return false
        end
    end

    # Check if the (photon, photon) angles pass the cut
    for ph1=1:OUTGOING_COUNT, ph2=ph1+1:OUTGOING_COUNT
        p_ph1 = ps_out[ph1, :]
        p_ph2 = ps_out[ph2, :]
        cos_num = dot(p_ph1[X:Z], p_ph2[X:Z])
        cos_denom = p_ph1[E] * p_ph2[E]
        if cos_num > cut.b_cut * cos_denom
            return false
        end
    end

    # TODO: Not implemented yet
    throw(AssertionError("Not implemented yet"))
end

end