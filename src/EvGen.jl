# Depends on Errors.jl, EvData.jl, LinAlg.jl, Numeric.jl and Random.jl being
# include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Event generation facilities"
module EvGen

using ..Errors: @enforce
using ..EvData: Event, INCOMING, OUTGOING, NUM_INCOMING, NUM_OUTGOING,
                NUM_PARTICLES
using ..LinAlg: X, Y, Z, XYZ, E
using ..Numeric: Float
using ..Random: RandomGenerator, random!
using LinearAlgebra: ⋅
using StaticArrays: MMatrix, MVector, SMatrix, SVector, @SMatrix, @SVector

export EventGenerator, generate_event!


# FIXME: How to express that NO_PHOTON_SORTING should be a Bool? Merely typing
#        NO_PHOTON_SORTING<:Bool doesn't work, calling constructor with a
#        Val(true) argument for a no_photon_sorting leads to a method mismatch.
"Generator of ee -> ppp events"
struct EventGenerator{NO_PHOTON_SORTING}
    "Total center-of-mass energy of the collision"
    e_tot::Float

    "Weight of generated events"
    event_weight::Float

    # FIXME: Need to specify SMatrix length to avoid type instability?
    "Incoming electron and positron momenta"
    incoming_momenta::SMatrix{NUM_INCOMING, 4, Float, NUM_INCOMING*4}

    "Disable sorting of photons by energy"
    no_photon_sorting::Val{NO_PHOTON_SORTING}
end


# === CONSTRUCTION ===

"""
Initialize event generation for a center-of-mass energy of e_tot.

Combines former functionality of ppp constructor and IBEGIN-based lazy
initialization from the original C++ 3photons code.
"""
function EventGenerator(
    e_tot::Float,
    no_photon_sorting::Val{NO_PHOTON_SORTING};
    jit_warmup::Bool=false
) where NO_PHOTON_SORTING
    # Check on the number of particles. The check for N<101 is gone since unlike
    # the original RAMBO, we don't use arrays of hardcoded size.
    @enforce (NUM_OUTGOING > 1) "At least one particle should be coming out"

    # As currently written, this code only works for two incoming particles
    @enforce (NUM_INCOMING == 2) "This code assumes two incoming particles"

    # Compute some numerical constants. Replaces the lazy initialization from
    # the original RAMBO code with something less branchy.
    if !jit_warmup
        println("IBegin")
    end
    z = (NUM_OUTGOING-1) * log(π/2)  # Replaces Z[INP-1] in the original code
    for k = 2:NUM_OUTGOING-1
        z -= 2 * log(k-1)
    end
    z -= log(NUM_OUTGOING-1)

    # NOTE: The check on total energy is gone, because we only generate massless
    #       photons and so the total energy will always be enough.
    #       Counting of nonzero masses is also gone because it was unused.

    # All generated events will have the same weight: pre-compute it
    ln_weight = (2*NUM_OUTGOING - 4) * log(e_tot) + z
    @enforce (ln_weight ≥ -180 && ln_weight ≤ 174)
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
        no_photon_sorting
    )
end


# === EVENT GENERATION ===

"""
Generate massless outgoing 4-momenta in infinite phase space

The output momenta are provided as a matrix where rows are 4-momentum
components (Px, Py, Pz, E) and columns are particles.
"""
function generate_event_raw!(rng::RandomGenerator)::SMatrix{4, NUM_OUTGOING, Float}
    # This implementation targets maximal reproducibility with respect to the
    # original 3photons program, at the expense of performance.
    #
    # FIXME: Provide and expose a port of Rust version's faster-evgen mode

    # Generate the basic random parameters of the particles
    # (This code is convoluted because it replicates the RNG call order of the
    # original 3photons program, which itself isn't so straightforward)
    cos_θ_idx = 1
    φ_idx = 2
    exp_min_e_idx = 3
    params = @SMatrix [
        if row == cos_θ_idx
            2*random!(rng) - 1
        elseif row == φ_idx
            2π * random!(rng)
        elseif row == exp_min_e_idx
            random!(rng) * random!(rng)
        else
            throw(AssertionError("Unexpected parameter"))
        end
        for row=1:3, _part=1:NUM_OUTGOING
    ]
    cos_θ = params[cos_θ_idx, :]
    φ = params[φ_idx, :]
    exp_min_e = params[exp_min_e_idx, :]

    # Compute the outgoing momenta
    # NOTE: Unlike Rust, Julia benefits from manual sincos optimization here
    sincos_φ = sincos.(φ)
    sin_θ = sqrt.(1 .- cos_θ.^2)
    energy = -log.(nextfloat.(exp_min_e))
    @SMatrix [
        energy[part] *
            if coord == X
                sin_θ[part] * sincos_φ[part][1]
            elseif coord == Y
                sin_θ[part] * sincos_φ[part][2]
            elseif coord == Z
                cos_θ[part]
            elseif coord == E
                1
            else
                throw(AssertionError("Unexpected coordinate"))
            end
        for coord=1:4, part=1:NUM_OUTGOING
    ]
end


"""
Use a highly specialized version of the RAMBO (RAndom Momenta
Beautifully Organized) algorithm from S.D. Ellis, R. Kleiss and W.J.
Stirling to generate the 4-momenta of the three outgoing photons.

All events have the same weight, it can be queried via evgen.event_weight.

The 4-momenta of output photons are sorted by decreasing energy.
"""
function generate_event!(
    rng::RandomGenerator,
    evgen::EventGenerator{NO_PHOTON_SORTING}
)::Event{NO_PHOTON_SORTING} where NO_PHOTON_SORTING
    # Generate massless outgoing momenta in infinite phase space
    q = generate_event_raw!(rng)

    # Calculate the parameters of the conformal transformation
    r = @SVector [ sum(q[coord, :]) for coord=1:4 ]
    r_norm² = r[E]^2 - r[XYZ]⋅r[XYZ]  # FIXME: No squared norm func?
    α = evgen.e_tot / r_norm²
    r_norm = √r_norm²
    β = 1 / (r_norm + r[E])

    # Perform the conformal transformation from Q's to output 4-momenta
    tr_q = transpose(q)
    rq = tr_q[:, XYZ] * r[XYZ]
    p_e_svector = α * (r[E] * tr_q[:, E] - rq)
    p_e = NO_PHOTON_SORTING ? p_e_svector : MVector(p_e_svector)
    b_rq_e = β * rq - tr_q[:, E]
    p_xyz_smatrix = α * (r_norm * tr_q[:, XYZ] + b_rq_e * transpose(r[XYZ]))
    p_xyz = NO_PHOTON_SORTING ? p_xyz_smatrix : MMatrix(p_xyz_smatrix)

    # Sort the output 4-momenta in order of decreasing energy
    if !NO_PHOTON_SORTING
        for par1=1:NUM_OUTGOING-1, par2=par1+1:NUM_OUTGOING
            if p_e[par2] > p_e[par1]
                p_e[par1], p_e[par2] = p_e[par2], p_e[par1]
                p_xyz[par1, :], p_xyz[par2, :] = p_xyz[par2, :], p_xyz[par1, :]
            end
        end
    end

    # Build the final event data: incoming momenta + output 4-momenta
    #
    # FIXME: Discuss with StaticArrays devs why this is 2x faster than the
    #        straightforward `vcat(evgen.incoming_momenta, hcat(p_xyz, p_e))`...
    #
    ps = zeros(MMatrix{NUM_PARTICLES, 4})
    ps[INCOMING, :] = evgen.incoming_momenta
    ps[OUTGOING, XYZ] = p_xyz
    ps[OUTGOING, E] = p_e
    Event(
        SMatrix(ps),
        evgen.no_photon_sorting
    )
end

end