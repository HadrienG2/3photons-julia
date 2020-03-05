# Depends on Errors.jl, EvData.jl, LinAlg.jl and Numeric.jl being include-d
# beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Mechanism to apply a cut to generated events"
module EvCut

using ..Errors: @enforce
using ..EvData: Event, INCOMING_E₋, NUM_OUTGOING, OUTGOING, min_photon_energy
using ..LinAlg: XYZ, E
using ..Numeric: Float
using LinearAlgebra: ×, ⋅, norm

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
    p_e₋ = event[INCOMING_E₋, :]
    ps_out = event[OUTGOING, :]

    # Check if the (beam, photon) angles pass the cut
    cos_nums = ps_out[:, XYZ] * p_e₋[XYZ]
    cos_denoms = ps_out[:, E] * p_e₋[E]
    for (num, denom) ∈ zip(cos_nums, cos_denoms)
        if abs(num) > cut.a_cut * denom
            return false
        end
    end

    # Check if the (photon, photon) angles pass the cut
    for ph1=1:NUM_OUTGOING, ph2=ph1+1:NUM_OUTGOING
        p_ph1 = ps_out[ph1, :]
        p_ph2 = ps_out[ph2, :]
        cos_num = p_ph1[XYZ] ⋅ p_ph2[XYZ]
        cos_denom = p_ph1[E] * p_ph2[E]
        if cos_num > cut.b_cut * cos_denom
            return false
        end
    end

    # Compute a vector which is normal to the outgoing photon plane
    # This notion is only valid because we have three output photons
    @enforce (NUM_OUTGOING == 3) "This code assumes 3 outgoing particles"
    n_ppp = ps_out[1, XYZ] × ps_out[2, XYZ]

    # Compute the cosine of the angle between the beam and this vector
    cos_num = p_e₋[XYZ] ⋅ n_ppp
    cos_denom = p_e₋[E] * norm(n_ppp)

    # Check if the (beam, normal to photon plane) angle passes the cut
    if abs(cos_num) < cut.sin_cut * cos_denom
        return false
    end

    # If all checks passed, we're good
    true
end

end