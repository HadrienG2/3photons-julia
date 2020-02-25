# Depends on Errors.jl, LinAlg.jl, Numeric.jl and Random.jl being include-d
# beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Event generation and storage"
module EvGen

using ..Errors: @enforce
using ..LinAlg: X, Y, Z, XYZ, E
using ..Numeric: Float
using ..Random: RandomGenerator, random!
using LinearAlgebra: ⋅
using StaticArrays: MMatrix, MVector, SMatrix, SVector, @SMatrix, @SVector

export EventGenerator, electron_momentum, generate_event!, outgoing_momenta,
       min_photon_energy


"Number of incoming particles"
const NUM_INCOMING = 2

"Number of outgoing particles (replaces original INP)"
const NUM_OUTGOING = 3

"Number of particles in an event"
const NUM_PARTICLES = NUM_INCOMING + NUM_OUTGOING

"Number of possible spin values of the outgoing particles"
const NUM_SPINS = 2


"Generator of ee -> ppp events"
struct EventGenerator
    "Total center-of-mass energy of the collision"
    e_tot::Float

    "Weight of generated events"
    event_weight::Float

    # FIXME: Need to specify SMatrix length to avoid type instability?
    "Incoming electron and positron momenta"
    incoming_momenta::SMatrix{NUM_INCOMING, 4, Float, NUM_INCOMING*4}
end


# === CONSTRUCTION ===

"""
Initialize event generation for a center-of-mass energy of e_tot.

Combines former functionality of ppp constructor and IBEGIN-based lazy
initialization from the original C++ 3photons code.
"""
function EventGenerator(e_tot::Float; jit_warmup::Bool=false)
    # Check on the number of particles. The check for N<101 is gone since unlike
    # the original RAMBO, we don't use arrays of hardcoded size.
    @enforce (NUM_OUTGOING > 1)

    # As currently written, this code only works for two incoming particles
    @enforce (NUM_INCOMING == 2)

    # Compute some numerical constants. Replaces the lazy initialization from
    # the original RAMBO code with something less branchy.
    if !jit_warmup
        println("IBegin")
    end
    z = (NUM_OUTGOING-1) * log(π/2)  # Replaces Z[INP-1] in the original code
    for k = 2:NUM_OUTGOING-1
        z -= 2 * log(k - 1)
    end
    z -= log(NUM_OUTGOING - 1)

    # NOTE: The check on total energy is gone, because we only generate massless
    #       photons and so the total energy will always be enough.
    #       Counting of nonzero masses is also gone because it was unused.

    # All generated events will have the same weight: pre-compute it
    ln_weight = (2 * NUM_OUTGOING - 4) * log(e_tot) + z
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


# === GENERATED EVENTS ===

"Row of the incoming electron in the event data matrix"
const INCOMING_E₋ = 1

"Row of the incoming positron in the event data matrix"
const INCOMING_E₊ = 2


"""
Storage for ee -> ppp event data

Encapsulates a vector of incoming and outgoing 4-momenta.
"""
const Event = SMatrix{NUM_PARTICLES, 4, Float}


"Extract the electron 4-momentum"
function electron_momentum(event::Event)::SVector{4, Float}
    event[INCOMING_E₋, :]
end


"Access the outgoing 4-momenta"
function outgoing_momenta(event::Event)::SMatrix{NUM_OUTGOING, 4, Float}
    OUTGOING = SVector{NUM_OUTGOING}(NUM_INCOMING+1:NUM_PARTICLES)
    event[OUTGOING, :]
end


"Minimal outgoing photon energy"
function min_photon_energy(event::Event)::Float
    # Use the fact that photons are sorted by decreasing energy
    #
    # FIXME: Revise this code once option to disable sorting returns.
    #
    outgoing_momenta(event)[NUM_OUTGOING, E]
end


# === EVENT GENERATION ===

"""
Generate massless outgoing 4-momenta in infinite phase space

The output momenta are provided as a matrix where rows are 4-momentum
components (Px, Py, Pz, E) and columns are particles.
"""
function generate_event_raw!(rng::RandomGenerator)::SMatrix{4, NUM_OUTGOING, Float}
    # Must be generalized to support NUM_OUTGOING != 3
    @enforce (NUM_OUTGOING == 3) "This function assumes 3 outgoing particles"

    # This implementation targets maximal reproducibility with respect to the
    # original 3photons program, at the expense of performance.
    #
    # FIXME: Provide and expose a port of Rust version's faster-evgen impl'

    # Generate the basic random parameters of the particles
    # (This code is convoluted because it replicates the RNG call order of the
    # original 3photons program, which itself isn't so straightforward)
    cos_θ_idx = 1
    φ_idx = 2
    exp_min_e_idx = 3
    params = @SMatrix [
        if coord == cos_θ_idx
            2*random!(rng) - 1
        elseif coord == φ_idx
            2π * random!(rng)
        elseif coord == exp_min_e_idx
            random!(rng) * random!(rng)
        else
            throw(AssertionError("Unexpected coordinate"))
        end
        for coord=1:3, _par=1:NUM_OUTGOING
    ]
    cos_θ = params[cos_θ_idx, :]
    φ = params[φ_idx, :]
    exp_min_e = params[exp_min_e_idx, :]

    # Compute the outgoing momenta
    # NOTE: Unlike Rust, Julia needs manual sincos optimization here
    sincos_φ = sincos.(φ)
    sin_θ = sqrt.(1 .- cos_θ.^2)
    energy = -log.(exp_min_e .+ eps(Float))
    @SMatrix [
        energy[par] *
            if coord == X
                sin_θ[par] * sincos_φ[par][1]
            elseif coord == Y
                sin_θ[par] * sincos_φ[par][2]
            elseif coord == Z
                cos_θ[par]
            elseif coord == E
                1
            else
                throw(AssertionError("Unexpected coordinate"))
            end
        for coord=1:4, par=1:NUM_OUTGOING
    ]
end


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
    r_norm² = r[E]^2 - r[XYZ]⋅r[XYZ]  # FIXME: No squared norm func?
    α = evgen.e_tot / r_norm²
    r_norm = √r_norm²
    β = 1 / (r_norm + r[E])

    # Perform the conformal transformation from Q's to output 4-momenta
    #
    # FIXME: Switch to SVector in no-photon-sorting mode, once added, as that
    #        will rid us of a few remaining heap allocations.
    #
    tr_q = transpose(q)
    rq = tr_q[:, XYZ] * r[XYZ]
    p_e = MVector{NUM_OUTGOING}(α * (r[E] * tr_q[:, E] - rq))
    b_rq_e = β * rq - tr_q[:, E]
    p_xyz = MMatrix{NUM_OUTGOING, 3}(
        α * (r_norm * tr_q[:, XYZ] + b_rq_e * transpose(r[XYZ]))
    )

    # Sort the output 4-momenta in order of decreasing energy
    #
    # FIXME: Make this optional, as in Rust version
    #
    for par1=1:NUM_OUTGOING, par2=par1+1:NUM_OUTGOING
        if p_e[par2] > p_e[par1]
            p_e[par1], p_e[par2] = p_e[par2], p_e[par1]
            p_xyz[par1, :], p_xyz[par2, :] = p_xyz[par2, :], p_xyz[par1, :]
        end
    end

    # Build the final event: incoming momenta + output 4-momenta
    #
    # FIXME: Discuss with StaticArrays devs why this is 2x faster than an
    #        idiomatic `vcat(evgen.incoming_momenta, hcat(p_xyz, p_e))`...
    #
    res = zeros(MMatrix{NUM_PARTICLES, 4})
    INCOMING = SVector{NUM_INCOMING}(1:NUM_INCOMING)
    OUTGOING = SVector{NUM_OUTGOING}(NUM_INCOMING+1:NUM_INCOMING+NUM_OUTGOING)
    res[INCOMING, :] = evgen.incoming_momenta
    res[OUTGOING, XYZ] = p_xyz
    res[OUTGOING, E] = p_e
    SMatrix(res)
end

end