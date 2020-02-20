# Depends on Errors.jl, LinAlg.jl, Numeric.jl and Random.jl being include-d
# beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Event generation and storage"
module EvGen

using ..Errors: @enforce
using ..LinAlg: X, Y, Z, E
using ..Numeric: Float
using ..Random: RandomGenerator, random!
using LinearAlgebra: dot
using StaticArrays: SMatrix, SVector, @SMatrix, @SVector

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
Generate massless outgoing 4-momenta in infinite phase space

The output momenta are provided as a matrix where rows are 4-momentum
components (Px, Py, Pz, E) and columns are particles.
"""
function generate_event_raw!(rng::RandomGenerator)::SMatrix{4, OUTGOING_COUNT, Float}
    # Must be generalized to support OUTGOING_COUNT != 3
    @enforce (OUTGOING_COUNT == 3) "This function assumes 3 outgoing particles"

    # This implementation targets maximal reproducibility with respect to the
    # original 3photons program, at the expense of performance.
    #
    # TODO: Provide and expose a port of Rust version's faster-evgen impl'

    # Generate the basic random parameters of the particles
    # (This code is convoluted because it replicates the RNG call order of the
    # original 3photons program, which itself isn't so straightforward)
    cos_theta_idx = 1
    phi_idx = 2
    exp_min_e_idx = 3
    params = @SMatrix [
        if coord == cos_theta_idx
            2 * random!(rng) - 1
        elseif coord == phi_idx
            2π * random!(rng)
        elseif coord == exp_min_e_idx
            random!(rng) * random!(rng)
        else
            throw(AssertionError("Unexpected coordinate"))
        end
        for coord=1:3, _par=1:OUTGOING_COUNT
    ]
    cos_theta = params[cos_theta_idx, :]
    phi = params[phi_idx, :]
    exp_min_e = params[exp_min_e_idx, :]

    # Compute the outgoing momenta
    cos_phi = map(cos, phi)
    sin_phi = map(sin, phi)
    sin_theta = map(c -> sqrt(1 - c*c), cos_theta)
    energy = map(e_me -> -log(e_me + eps(Float)), exp_min_e)
    @SMatrix [
        energy[par] *
            if coord == X
                sin_theta[par] * sin_phi[par]
            elseif coord == Y
                sin_theta[par] * cos_phi[par]
            elseif coord == Z
                cos_theta[par]
            elseif coord == E
                1
            else
                throw(AssertionError("Unexpected coordinate"))
            end
        for coord=1:4, par=1:OUTGOING_COUNT
    ]
end


"""
Storage for ee -> ppp event data

Encapsulates a vector of incoming and outgoing 4-momenta.
"""
const Event = SMatrix{PARTICLE_COUNT, 4, Float}


"""
Use a highly specialized version of the RAMBO (RAndom Momenta
Beautifully Organized) algorithm from S.D. Ellis, R. Kleiss and W.J.
Stirling to generate the 4-momenta of the three outgoing photons.

All events have the same weight, it can be queried via evgen.event_weight.

The 4-momenta of output photons are sorted by decreasing energy.
"""
function generate_event!(rng::RandomGenerator, evgen::EventGenerator)::Event
    # Generate massless outgoing momenta in infinite phase space
    q = generate_event_raw!(rng)

    # Calculate the parameters of the conformal transformation
    r = @SVector [ sum(q[coord, :]) for coord=1:4 ]
    r_norm_2 = r[E]^2 - dot(r[X:Z], r[X:Z])  # FIXME: No squared norm?
    alpha = evgen.e_tot / r_norm_2
    r_norm = sqrt(r_norm_2)
    beta = 1 / (r_norm + r[E])

    # TODO: Not implemented yet
    throw(AssertionError("Not implemented yet"))
end

end