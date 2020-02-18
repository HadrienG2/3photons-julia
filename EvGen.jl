# Depends on Errors.jl, Numeric.jl and Random.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Event generation and storage"
module EvGen

using ..Errors: @enforce
using ..Numeric: Float
using ..Random: RandomGenerator
using StaticArrays: SMatrix, @SMatrix

export EventGenerator, generate_event!


"Number of incoming particles"
const INCOMING_COUNT = 2

"Number of outgoing particles (replaces original INP)"
const OUTGOING_COUNT = 3

"Number of particles in an event"
const PARTICLE_COUNT = INCOMING_COUNT + OUTGOING_COUNT


"Generator of ee -> ppp events"
struct EventGenerator
    "Total center-of-mass energy of the collision"
    e_tot::Float

    "Weight of generated events"
    event_weight::Float

    "Incoming electron and positron momenta"
    incoming_momenta::SMatrix{INCOMING_COUNT, 4, Float}
end


"""
Storage for ee -> ppp event data

Encapsulates a vector of incoming and outgoing 4-momenta.
"""
const Event = SMatrix{PARTICLE_COUNT, 4, Float}


# === CONSTRUCTION ===

"""
Initialize event generation for a center-of-mass energy of e_tot.

Combines former functionality of ppp constructor and IBEGIN-based lazy
initialization from the original C++ 3photons code.
"""
function EventGenerator(e_tot::Real)
    # Check on the number of particles. The check for N<101 is gone since unlike
    # the original RAMBO, we don't use arrays of hardcoded size.
    @enforce (OUTGOING_COUNT > 1)

    # As currently written, this code only works for two incoming particles
    @enforce (INCOMING_COUNT == 2)

    # Compute some numerical constants. Replaces the lazy initialization from
    # the original RAMBO code with something less branchy.
    println("IBegin")
    z = (OUTGOING_COUNT-1) * log(π/2)  # Replaces Z[INP-1] in the original code
    for k = 2:OUTGOING_COUNT-1
        z -= 2 * log(k - 1)
    end
    z -= log(OUTGOING_COUNT - 1)

    # NOTE: The check on total energy is gone, because we only generate massless
    #       photons and so the total energy will always be enough.
    #       Counting of nonzero masses is also gone because it was unused.

    # All generated events will have the same weight: pre-compute it
    ln_weight = (2 * OUTGOING_COUNT - 4) * log(e_tot) + z
    @enforce (ln_weight >= -180 && ln_weight < 174)
    event_weight = exp(ln_weight)

    # Compute the incoming particle momenta
    half_e_tot = e_tot / 2
    incoming_momenta = @SMatrix [ -half_e_tot 0. 0. half_e_tot ;
                                   half_e_tot 0. 0. half_e_tot ]

    # Construct and return the output data structure
    EventGenerator(
        e_tot,
        event_weight,
        incoming_momenta,
    )
end


# === EVENT GENERATION ===

"""
Use a highly specialized version of the RAMBO (RAndom Momenta
Beautifully Organized) algorithm from S.D. Ellis, R. Kleiss and W.J.
Stirling to generate the 4-momenta of the three outgoing photons.

All events have the same weight, it can be queried via evgen.event_weight.

The 4-momenta of output photons are sorted by decreasing energy.
"""
function generate_event!(rng::RandomGenerator, evgen::EventGenerator)::Event
    # TODO: Not implemented yet
    throw(AssertionError("Not implemented yet"))
end

end