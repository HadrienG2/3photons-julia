# Depends on EvGen.jl and Numeric.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Mechanism to apply a cut to generated events"
module EvCut

using ..EvGen: Event
using ..Numeric: Float

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
    # TODO: Not implemented yet
    throw(AssertionError("Not implemented yet"))
end

end