# Depends on Errors.jl, EvGen.jl, LinAlg.jl and Numeric.jl being include-d
# beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Facilities for computing spinor products"
module Spinor

using ..Errors: @enforce
using ..EvGen: Event, NUM_INCOMING, NUM_OUTGOING, NUM_PARTICLES, NUM_SPINS
using ..LinAlg: E, X, Y, Z
using ..Numeric: Float
using StaticArrays: SMatrix, SVector, @SMatrix, @SVector

export NUM_HELICITIES, SpinorProducts


# PhotonHelicities declaration and SpinorProducts methods are very specific to
# our simulation's e+e- -> ppp configuration.
@enforce (NUM_INCOMING == 2) && (NUM_OUTGOING == 3) && (NUM_SPINS == 2) """
This code assumes two incoming particles and three outgoing photons
"""


"Massless 4-momenta spinor inner products"
const SpinorProducts = SMatrix{NUM_PARTICLES, NUM_PARTICLES, Complex{Float}}

"Build spinor products from previously generated particle 4-momenta"
function SpinorProducts(event::Event)
    # Compute the spinor products (method from M. Mangano and S. Parke)
    # NOTE: No fancy √ syntax here, for some reason it doesn't compose with '.'
    xx = sqrt.(event[:, E] + event[:, Z])
    fx = @SVector [
        if xx[par] > eps(Float)
            complex(event[par, X], event[par, Y]) / xx[par]
        else
            Complex(sqrt(2 * event[par, E]))
        end
        for par=1:NUM_PARTICLES
    ]

    # Fill up the Gram matrix
    # TODO: Can we leverage antisymmetry + zero diagonal better?
    @SMatrix [
        fx[i] * xx[j] - fx[j] * xx[i]
        for i=1:NUM_PARTICLES, j=1:NUM_PARTICLES
    ]
end


"Output photon helicities (M is - and P is +)"
@enum PhotonHelicities begin
    MMM = 0b000
    MMP = 0b001
    MPM = 0b010
    MPP = 0b011
    PMM = 0b100
    PMP = 0b101
    PPM = 0b110
    PPP = 0b111
end

"Number of photon helicity configurations"
const NUM_HELICITIES = NUM_OUTGOING^NUM_SPINS

end