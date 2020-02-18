# Depends on numeric.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Intermediary results of the computation for one event"
module ResCont

using ..Numeric: Float
using StaticArrays: SVector

export NUM_RESULTS, ResultVector


"Number of results (matrix elements)"
const NUM_RESULTS = 5

"Storage for matrix elements"
const ResultVector = SVector{NUM_RESULTS, Float}

end