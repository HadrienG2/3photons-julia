# Depends on Coupling.jl, Errors.jl, EvData.jl and Numeric.jl being include-d
# beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Intermediary results of the computation for one event"
module ResCont

using ..Coupling: Couplings
using ..Errors: @enforce
using ..EvData: Event, NUM_OUTGOING
using ..Numeric: Float
using ..Spinor: 𝛼_amp, 𝛽₊_amp, 𝛽₋_amp, MMM, NUM_HELICITIES, PhotonHelicities,
                PPP, SpinorProducts
using StaticArrays: SMatrix, SVector, @SMatrix, @SVector

export A, B₊, B₋, I_MX, m²_sums, NUM_RESULTS, R_MX, ResultContibution,
       ResultVector


"Number of results (matrix elements)"
const NUM_RESULTS = 5

"Storage for matrix elements"
const ResultVector = SVector{NUM_RESULTS, Float}

"Index of the electromagnetic matrix element"
const A = 1

"Index of the positive electroweak matrix element"
const B₊ = 2

"Index of the negative electroweak matrix element"
const B₋ = 3

"Index of the real part of the mixed matrix element"
const R_MX = 4

"Index of the imaginary part of the mixed matrix element"
const I_MX = 5


# FIXME: Need to specify SMatrix length to avoid type instability?
"Array of square matrix elements contribution with detail of helicities"
const ResultContribution = SMatrix{NUM_RESULTS, NUM_HELICITIES, Float, NUM_RESULTS*NUM_HELICITIES}


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
    mixed_amps = 2(𝛼_amps .* conj(𝛽₊_amps))

    # Compute the matrix elements
    @SMatrix [
        if contrib == A
            abs2(𝛼_amps[hel])
        elseif contrib == B₊
            abs2(𝛽₊_amps[hel])
        elseif contrib == B₋
            abs2(𝛽₋_amps[hel])
        elseif contrib == R_MX
            real(mixed_amps[hel])
        elseif contrib == I_MX
            imag(mixed_amps[hel])
        else
            throw(AssertionError("Unexpected contribution"))
        end
        for contrib=1:NUM_RESULTS, hel=1:NUM_HELICITIES
    ]
end


"Compute the sums of the squared matrix elements for each contribution"
function m²_sums(rc::ResultContribution)::ResultVector
    @SVector [ sum(rc[i, :]) for i=1:NUM_RESULTS ]
end

end