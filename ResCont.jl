# Depends on Coupling.jl, EvGen.jl and Numeric.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Intermediary results of the computation for one event"
module ResCont

using ..Coupling: Couplings
using ..EvGen: Event
using ..Numeric: Float
using ..Spinor: NUM_HELICITIES
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
    # TODO: Not implemented yet
    throw(AssertionError("Not implemented yet"))
end

end