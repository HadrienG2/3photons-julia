# Depends on Coupling.jl, Errors.jl, EvData.jl and Numeric.jl being include-d
# beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Manipulation of matrix elements"
module MatElems

using ..Coupling: Couplings
using ..Errors: @enforce
using ..EvData: Event, NUM_OUTGOING
using ..Numeric: Float
using ..Spinor: 𝛼_amp, 𝛽₊_amp, 𝛽₋_amp, HELICITIES, NUM_HELICITIES,
                SpinorProducts
using StaticArrays: SMatrix, SVector, @SMatrix, @SVector

export A, B₊, B₋, I_MX, m²_sums, MEsContibutions, MEsVector, NUM_MAT_ELEMS, R_MX


# === MATRIX ELEMENTS ===

"Number of matrix elements"
const NUM_MAT_ELEMS = 5

"Storage for per-matrix element data"
const MEsVector = SVector{NUM_MAT_ELEMS, Float}

"Index of the electromagnetic element"
const A = 1

"Index of the positive electroweak element"
const B₊ = 2

"Index of the negative electroweak element"
const B₋ = 3

"Index of the real part of the mixed element"
const R_MX = 4

"Index of the imaginary part of the mixed element"
const I_MX = 5


# === PER-HELICITY CONTRIBUTIONS TO MATRIX ELEMENTS ===

# FIXME: Need to specify SMatrix length to avoid type instability?
"Array of squared matrix elements contributions, with detail of helicities"
const MEsContributions = SMatrix{NUM_MAT_ELEMS, NUM_HELICITIES, Float, NUM_MAT_ELEMS*NUM_HELICITIES}


"Construct the matrix element contributions from the event data"
function MEsContributions(couplings::Couplings, event::Event)
    # This code is very specific to the current problem definition
    @enforce (NUM_MAT_ELEMS == 5) "This code assumes 5 matrix elements"

    # Compute spinor inner products
    spinor = SpinorProducts(event)

    # Compute the helicity amplitudes, formerly known as a, b_p and b_m, for
    # each possible output spin configuration
    𝛼_amps = map(hel-> 𝛼_amp(spinor, hel) * couplings.g_𝛼, HELICITIES)
    𝛽₊_amps = map(hel -> 𝛽₊_amp(spinor, hel) * couplings.g_𝛽₊, HELICITIES)
    𝛽₋_amps = map(hel -> 𝛽₋_amp(spinor, hel) * couplings.g_𝛽₋, HELICITIES)
    mixed_amps = 2(𝛼_amps .* conj(𝛽₊_amps))

    # Compute the matrix elements
    @SMatrix [
        if elem == A
            abs2(𝛼_amps[hel])
        elseif elem == B₊
            abs2(𝛽₊_amps[hel])
        elseif elem == B₋
            abs2(𝛽₋_amps[hel])
        elseif elem == R_MX
            real(mixed_amps[hel])
        elseif elem == I_MX
            imag(mixed_amps[hel])
        else
            throw(AssertionError("Unexpected contribution"))
        end
        for elem=1:NUM_MAT_ELEMS, hel=1:NUM_HELICITIES
    ]
end


"Compute the sums of the squared matrix elements for each contribution"
function m²_sums(rc::MEsContributions)::MEsVector
    @SVector [ sum(rc[elem, :]) for elem=1:NUM_MAT_ELEMS ]
end

end