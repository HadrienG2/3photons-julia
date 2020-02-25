# Depends on Errors.jl, EvGen.jl, LinAlg.jl and Numeric.jl being include-d
# beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Facilities for computing spinor products"
module Spinor

using ..Errors: @enforce
using ..EvGen: Event, INCOMING_E₋, INCOMING_E₊, NUM_INCOMING, NUM_OUTGOING,
               NUM_PARTICLES, NUM_SPINS
using ..LinAlg: E, X, Y, Z
using ..Numeric: Float
using StaticArrays: SMatrix, SVector, @SMatrix, @SVector

export 𝛼_amp, 𝛽₊_amp, 𝛽₋_amp, NUM_HELICITIES, SpinorProducts


# PhotonHelicities declaration and SpinorProducts methods are very specific to
# our simulation's e+e- -> ppp configuration.
@enforce (NUM_INCOMING == 2) && (NUM_OUTGOING == 3) && (NUM_SPINS == 2) """
This code assumes two incoming particles and three outgoing photons
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
const NUM_HELICITIES = NUM_SPINS^NUM_OUTGOING


"Massless 4-momenta spinor inner products"
const SpinorProducts = SMatrix{NUM_PARTICLES, NUM_PARTICLES, Complex{Float}}


# === CONSTRUCTION ===

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


# === GRAM MATRIX ACCESSORS ===

function s(sx::SpinorProducts, i::Int, j::Int)::Complex{Float}
    sx[i, j]
end

function t(sx::SpinorProducts, i::Int, j::Int)::Complex{Float}
    -conj(sx[i, j])
end


# === AMPLITUDE COMPUTATIONS ===

# Internal computations from specific helicity configurations + indices

const E₋ = INCOMING_E₋
const E₊ = INCOMING_E₊

"Standard amplitude for helicities ++-"
function 𝛼_ppm(sx::SpinorProducts, k1::Int, k2::Int, k3::Int)::Complex{Float}
    -√8 * s(sx, E₋, E₊) * s(sx, E₋, k3)^2 /
        (s(sx, E₋, k1) * s(sx, E₋, k2) * s(sx, E₊, k1) * s(sx, E₊, k2))
end

"Standard amplitude for helicities +--"
function 𝛼_pmm(sx::SpinorProducts, k1::Int, k2::Int, k3::Int)::Complex{Float}
    -√8 * t(sx, E₋, E₊) * t(sx, E₊, k1)^2 /
        (t(sx, E₊, k2) * t(sx, E₊, k3) * t(sx, E₋, k2) * t(sx, E₋, k3))
end

"Anomalous amplitude 𝛽₊ for helicities ++-"
function 𝛽₊_ppm(sx::SpinorProducts, k1::Int, k2::Int, k3::Int)::Complex{Float}
    -√8 * t(sx, E₋, E₊) * (t(sx, k1, k2) * s(sx, k3, E₋))^2
end

"Anomalous amplitude 𝛽₊ for helicities +--"
function 𝛽₊_pmm(sx::SpinorProducts, k1::Int, k2::Int, k3::Int)::Complex{Float}
    -√8 * s(sx, E₋, E₊) * (t(sx, k1, E₊) * s(sx, k2, k3))^2
end

"Anomalous amplitude 𝛽₋ for helicities +++"
function 𝛽₋_ppp(sx::SpinorProducts, k1::Int, k2::Int, k3::Int)::Complex{Float}
    -√8 * s(sx, E₋, E₊) * ((t(sx, k1, k2) * t(sx, k3, E₊))^2 +
                           (t(sx, k1, k3) * t(sx, k2, E₊))^2 +
                           (t(sx, k2, k3) * t(sx, k1, E₊))^2)
end

"Anomalous amplitude 𝛽₋ for helicities ---"
function 𝛽₋_mmm(sx::SpinorProducts, k1::Int, k2::Int, k3::Int)::Complex{Float}
    -√8 * t(sx, E₋, E₊) * ((s(sx, k1, E₋) * s(sx, k2, k3))^2 +
                           (s(sx, k2, E₋) * s(sx, k1, k3))^2 +
                           (s(sx, k3, E₋) * s(sx, k1, k2))^2)
end


# External interface from helicities enum

"Standard amplitude for given photon helicities"
function 𝛼_amp(sx::SpinorProducts, helicities::PhotonHelicities)::Complex{Float}
    if helicities == MMM
        Complex(0.0)
    elseif helicities == MMP
        𝛼_pmm(sx, 5, 3, 4)
    elseif helicities == MPM
        𝛼_pmm(sx, 4, 3, 5)
    elseif helicities == MPP
        𝛼_ppm(sx, 4, 5, 3)
    elseif helicities == PMM
        𝛼_pmm(sx, 3, 4, 5)
    elseif helicities == PMP
        𝛼_ppm(sx, 3, 5, 4)
    elseif helicities == PPM
        𝛼_ppm(sx, 3, 4, 5)
    elseif helicities == PPP
        Complex(0.0)
    else
        throw(AssertionError("Unexpected helicities"))
    end
end

"Anomalous amplitude 𝛽₊ for given photon helicities"
function 𝛽₊_amp(sx::SpinorProducts, helicities::PhotonHelicities)::Complex{Float}
    if helicities == MMM
        Complex(0.0)
    elseif helicities == MMP
        𝛽₊_pmm(sx, 5, 3, 4)
    elseif helicities == MPM
        𝛽₊_pmm(sx, 4, 3, 5)
    elseif helicities == MPP
        𝛽₊_ppm(sx, 4, 5, 3)
    elseif helicities == PMM
        𝛽₊_pmm(sx, 3, 4, 5)
    elseif helicities == PMP
        𝛽₊_ppm(sx, 3, 5, 4)
    elseif helicities == PPM
        𝛽₊_ppm(sx, 3, 4, 5)
    elseif helicities == PPP
        Complex(0.0)
    else
        throw(AssertionError("Unexpected helicities"))
    end
end

"Anomalous amplitude 𝛽₋ for given photon helicities"
function 𝛽₋_amp(sx::SpinorProducts, helicities::PhotonHelicities)::Complex{Float}
    if helicities == MMM
        𝛽₋_mmm(sx, 3, 4, 5)
    elseif helicities == MMP
        Complex(0.0)
    elseif helicities == MPM
        Complex(0.0)
    elseif helicities == MPP
        Complex(0.0)
    elseif helicities == PMM
        Complex(0.0)
    elseif helicities == PMP
        Complex(0.0)
    elseif helicities == PPM
        Complex(0.0)
    elseif helicities == PPP
        𝛽₋_ppp(sx, 3, 4, 5)
    else
        throw(AssertionError("Unexpected helicities"))
    end
end

end