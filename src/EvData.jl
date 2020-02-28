# Depends on LinAlg.jl and Numeric.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Event definition and storage"
module EvData

using ..LinAlg: E
using ..Numeric: Float
using StaticArrays: SMatrix, SVector

export Event, INCOMING, INCOMING_E₋, INCOMING_E₊, min_photon_energy,
       NUM_INCOMING, NUM_OUTGOING, NUM_PARTICLES, NUM_SPINS


# === EVENT TYPE DEFINITION ===

"Number of incoming particles"
const NUM_INCOMING = 2

"Number of outgoing particles (replaces original INP)"
const NUM_OUTGOING = 3

"Number of particles in an event"
const NUM_PARTICLES = NUM_INCOMING + NUM_OUTGOING

"Number of possible spin values of the outgoing particles"
const NUM_SPINS = 2


# === EVENT DATA STORAGE ===

# FIXME: Need to specify SMatrix length to avoid type instability in structs?
"""
Storage for ee -> ppp event data

This is a matrix of incoming and outgoing 4-momenta.
"""
const Event = SMatrix{NUM_PARTICLES, 4, Float, NUM_PARTICLES*4}

"Row of the incoming electron in the event data matrix"
const INCOMING_E₋ = 1

"Row of the incoming positron in the event data matrix"
const INCOMING_E₊ = 2

"Rows of the incoming particles in the event data matrix"
const INCOMING = SVector{NUM_INCOMING}(1:NUM_INCOMING)

"Rows of the outgoing particles in the event data matrix"
const OUTGOING = SVector{NUM_OUTGOING}(NUM_INCOMING+1:NUM_PARTICLES)


"Minimal outgoing photon energy"
function min_photon_energy(event::Event)::Float
    # Access photon momenta
    outgoing_momenta = event[OUTGOING, :]

    # Use the fact that photons are sorted by decreasing energy
    #
    # FIXME: Revise this code once option to disable sorting returns.
    #
    outgoing_momenta[NUM_OUTGOING, E]
end

end