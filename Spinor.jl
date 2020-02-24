# Depends on Errors.jl and EvGen.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Facilities for computing spinor products"
module Spinor

using ..Errors: @enforce
using ..EvGen: NUM_OUTGOING, NUM_SPINS


# PhotonHelicities declaration is specific to the case of 3 output photons
@enforce (NUM_OUTGOING == 3) """
This code assumes three outgoing particles with spin +/-1
"""

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