"Mechanism to apply a cut to generated events"
module EvCut

include("numeric.jl")

using .Numeric: Float

export EventCut


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

end