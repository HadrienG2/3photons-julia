# Depends on Coupling.jl, Errors.jl, EvGen.jl and Numeric.jl being include-d
# beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Intermediary results of the computation for one event"
module ResCont

using ..Coupling: Couplings
using ..Errors: @enforce
using ..EvGen: Event, NUM_OUTGOING
using ..Numeric: Float
using ..Spinor: 𝛼_amp, 𝛽₊_amp, 𝛽₋_amp, MMM, NUM_HELICITIES, PhotonHelicities,
                PPP, SpinorProducts
using StaticArrays: SMatrix, SVector

export NUM_RESULTS, ResultContibution, ResultVector


"Number of results (matrix elements)"
const NUM_RESULTS = 5

"Storage for matrix elements"
const ResultVector = SVector{NUM_RESULTS, Float}


"Array of square matrix elements contribution with detail of helicities"
const ResultContribution = SMatrix{NUM_RESULTS, NUM_HELICITIES, Float}


"Construct the matrix element from the spinor products"
function ResultContribution(couplings::Couplings, event::Event)
    # This code is very specific to the current problem definition
    @enforce (NUM_OUTGOING == 3) "This code assumes 3 outgoing photons"
    @enforce (NUM_RESULTS == 5) "This code computes 5 results"

    # Compute spinor inner products
    spinor = SpinorProducts(event)

    # Compute the helicity amplitudes, formerly known as a, b_p and b_m, for
    # each possible output spin configuration
    hels = PhotonHelicities.(SVector{NUM_HELICITIES}(0b000:0b111))
    𝛼_amps = map(hel-> 𝛼_amp(spinor, hel) * couplings.g_𝛼, hels)
    𝛽₊_amps = map(hel -> 𝛽₊_amp(spinor, hel) * couplings.g_𝛽₊, hels)
    𝛽₋_amps = map(hel -> 𝛽₋_amp(spinor, hel) * couplings.g_𝛽₋, hels)

    # TODO: Not implemented yet
    throw(AssertionError("Not implemented yet"))
end

end