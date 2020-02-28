# Depends on Config.jl, Errors.jl, EvData.jl, Numeric.jl, ResAcc.jl and
# ResCont.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"""
This module contains everything that is needed to compute, store, and analyze
the final results: differential cross-section, sum & variance.
"""
module ResFin

using ..Config: Configuration
using ..Errors: @enforce
using ..EvData: NUM_SPINS
using ..Numeric: Float
using ..ResAcc: ResultsAccumulator
using ..ResCont: A, B₊, B₋, I_MX, NUM_RESULTS, R_MX
using LinearAlgebra: norm
using Printf: @printf
using StaticArrays: MMatrix, SMatrix, SVector, @MMatrix, @SMatrix, @SVector

export FinalResults, print_eric, print_fawzi, SP₋, SP₊


# FIXME: Need to specify SMatrix length to avoid type instability?
"""
Matrix of per-spin result contributions

Rows are spins, columns are result contributions (in the ResCont.jl sense)
"""
const PerSpinResults = SMatrix{NUM_SPINS, NUM_RESULTS, Float, NUM_SPINS*NUM_RESULTS};

"Index of negative spin data"
const SP₋ = 1

"Index of positive spin data"
const SP₊ = 2


"Final results of the simulation"
struct FinalResults
    "Number of integrated events"
    selected_events::UInt

    "Cross-section for each spin"
    spm²::PerSpinResults

    "Variance for each spin"
    vars::PerSpinResults

    "Total cross-section"
    σ::Float

    "Relative precision"
    prec::Float

    "Total variance"
    variance::Float

    "Beta minimum (???)"
    𝛽_min::Float

    "Statistical significance B+(pb-1/2) (???)"
    ss₊::Float

    "Incertitide associated with ss₊"
    inc_ss₊::Float

    "Statistical significance B-(pb-1/2) (???)"
    ss₋::Float

    "Incertitude associated with ss₋"
    inc_ss₋::Float

    # FIXME: How does one express "should be private" in Julia?
    "Configuration of the simulation"
    _cfg::Configuration
end


"Turn integrated simulation data into finalized results"
function FinalResults(acc::ResultsAccumulator)::FinalResults
    # This code depends on some aspects of the problem definition
    @enforce (NUM_SPINS == 2) "This code currently assumes 2 spins"
    @enforce (NUM_RESULTS == 5) "This code currently assumes 5 matrix elements"

    # Simulation configuration shorthands
    cfg = acc.cfg
    n_ev = cfg.num_events

    # Compute the relative uncertainties for one spin
    #
    # FIXME: Why doesn't (v_spm², v_var) ∈ zip(acc.spm², acc.vars) work?
    #
    vars_1spin = @SVector [
        begin
            v_spm² = acc.spm²[i]
            v_var = (acc.vars[i] - v_spm²^2 / n_ev) / (n_ev - 1)
            √(v_var / n_ev) / abs(v_spm² / n_ev)
        end
        for i=1:NUM_RESULTS
    ]

    # Copy for the opposite spin
    spm² = @MMatrix [
        acc.spm²[res]
        for _spin=1:NUM_SPINS, res=1:NUM_RESULTS
    ]
    vars = @SMatrix [
        vars_1spin[res]
        for _spin=1:NUM_SPINS, res=1:NUM_RESULTS
    ]

    # Electroweak polarisations factors for the 𝛽₊/𝛽₋ anomalous contribution
    pol₊ = -2 * cfg.sin²_w
    pol₋ = 1 + pol₊

    # Take polarisations into account
    spm²[SP₋, B₊:B₋] *= pol₋^2
    spm²[SP₊, B₊:B₋] *= pol₊^2
    spm²[SP₋, R_MX:I_MX] *= pol₋
    spm²[SP₊, R_MX:I_MX] *= pol₊

    # Flux factor (=1/2s for 2 initial massless particles)
    flux = 1 / (2 * cfg.e_tot^2)

    # Apply physical coefficients and Z⁰ propagator to each spin
    spm² *= acc.fact_com * flux * acc.norm_weight
    gm_Z⁰ = cfg.g_Z⁰ * cfg.m_Z⁰
    spm²[:, B₊:I_MX] *= acc.propag / gm_Z⁰
    spm²[:, B₊:B₋] /= gm_Z⁰
    spm²[:, R_MX] *= acc.ecart_pic

    # Compute other parts of the result
    𝛽_min = √(sum(spm²[:, A]) / sum(spm²[:, B₊]))

    ss_denom = sum(spm²[:, A])
    ss_norm = 1 / (2 * √ss_denom)

    ss₊ = sum(spm²[:, B₊]) * ss_norm
    ss₋ = sum(spm²[:, B₋]) * ss_norm

    inc_ss_common = norm(spm²[:, A] .* vars[:, A]) / (2 * abs(ss_denom))

    inc_ss₊ = norm(spm²[:, B₊] .* vars[:, B₊]) / abs(sum(spm²[:, B₊])) +
              inc_ss_common
    inc_ss₋ = norm(spm²[:, B₋] .* vars[:, B₋]) / abs(sum(spm²[:, B₋])) +
              inc_ss_common

    variance = (acc.variance - acc.σ^2 / n_ev) / (n_ev - 1)
    prec = √(variance / n_ev) / abs(acc.σ / n_ev)
    σ = acc.σ * flux

    # Return the final results
    #
    # FIXME: Isn't there any way to say which field we are talking about?
    #
    FinalResults(
        acc.selected_events,
        spm²,
        vars,
        σ,
        prec,
        variance,
        𝛽_min,
        ss₊,
        inc_ss₊,
        ss₋,
        inc_ss₋,
        cfg,
    )
end


"Display results using Eric's (???) parametrization"
function print_eric(results::FinalResults)
    @enforce (NUM_SPINS == 2) "This code assumes particles of spin +/-1"
    @enforce (NUM_RESULTS == 5) "This code assumes a 5-results configuration"

    cfg = results._cfg
    spm² = results.spm²

    µ_th = cfg.br_e₊_e₋ * cfg.convers / (8 * 9 * 5 * π^2 * cfg.m_Z⁰ * cfg.g_Z⁰)
    λ₀₍₋₎ = (spm²[SP₋, B₋] - spm²[SP₋, B₊]) / 2
    λ₀₍₊₎ = (spm²[SP₊, B₋] - spm²[SP₊, B₊]) / 2
    µ₀₍₋₎ = (spm²[SP₋, B₋] + spm²[SP₋, B₊]) / 2
    µ₀₍₊₎ = (spm²[SP₊, B₋] + spm²[SP₊, B₊]) / 2
    µ_num = (spm²[SP₋, B₊] + spm²[SP₋, B₋] + spm²[SP₊, B₊] + spm²[SP₊, B₋]) / 4

    println()
    println("       :        -          +")
    @printf("sigma0  : %.6f | %.6f\n", spm²[SP₋, A]/2, spm²[SP₊, A]/2)
    @printf("alpha0  : %.5e | %.4e\n", spm²[SP₋, I_MX]/2, spm²[SP₊, I_MX]/2)
    @printf("beta0   : %.0f | %.0f\n", spm²[SP₋, R_MX]/2, spm²[SP₊, R_MX]/2)
    @printf("lambda0 : %.4f | %.4f\n", λ₀₍₋₎, λ₀₍₊₎)
    @printf("mu0     : %.4f | %.5f\n", µ₀₍₋₎, µ₀₍₊₎)
    @printf("mu/lamb : %.5f | %.5f\n", µ₀₍₋₎/λ₀₍₋₎, µ₀₍₊₎/λ₀₍₊₎)
    @printf("mu (num): %.4f\n", µ_num)
    @printf("rapport : %.6f\n", µ_num/µ_th)
    @printf("mu (th) : %.4f\n", µ_th)
end


"""
Display Fawzi's (???) analytical results and compare them to the Monte Carlo
results that we have computed
"""
function print_fawzi(results::FinalResults)
    @enforce (NUM_RESULTS == 5) "This code assumes a 5-results configuration"

    cfg = results._cfg
    ev_cut = cfg.event_cut
    spm² = results.spm²
    vars = results.vars

    mre = cfg.m_Z⁰ / cfg.e_tot
    gre = cfg.g_Z⁰ * cfg.m_Z⁰ / cfg.e_tot^2
    x = 1 - mre^2
    sdz = complex(x, -gre) / (x^2 + gre^2)
    δ = (1 - ev_cut.b_cut) / 2
    ε = 2 * ev_cut.e_min / cfg.e_tot
    bra = cfg.m_Z⁰ / (3 * 6 * π^3 * 16 * 120)
    σ = 12π / cfg.m_Z⁰^2 * cfg.br_e₊_e₋ * cfg.g_Z⁰ * bra / cfg.e_tot^2 *
        (cfg.e_tot / cfg.m_Z⁰)^8 * abs2(sdz) * cfg.convers

    f₁ = 1 - 15ε^4 - 9/7 * (1 - 70ε^4)δ^2 + 6/7 * (1 + 70ε^4)δ^3
    g₁ = 1 - 30ε^4 - 9/7 * (1 - 70ε^4)δ - 90ε^4 * δ^2 - 1/7 * (1 - 420ε^4)δ^3
    g₂ = 1 - 25ε^4 - 6/7 * (1 - 70ε^4)δ - 3/7 * (1 + 210ε^4)δ^2 -
         8/21 * (1 - 52.5ε^4)δ^3
    g₃ = 1 - 195/11 * ε^4 - 18/77 * (1 - 7ε^4)δ - 9/11 * (9/7 - 70ε^4)δ^2 -
         8/11 * (1 - 105/11 * ε^4)δ^3

    ff = f₁ * (1 - ev_cut.sin_cut^3)
    gg = g₁ - 27/16 * g₂ * ev_cut.sin_cut + 11/16 * g₃ * ev_cut.sin_cut^3

    σ₊ = σ * (ff + 2gg)
    σ₋ = σ₊ + 2σ * gg

    mc₊ = sum(spm²[:, B₊]) / 4
    mc₋ = sum(spm²[:, B₋]) / 4
    incr₊ = norm(spm²[:, B₊] .* vars[:, B₊]) / abs(sum(spm²[:, B₊]))
    incr₋ = norm(spm²[:, B₋] .* vars[:, B₋]) / abs(sum(spm²[:, B₋]))

    println()
    println("s (pb) :   Sig_cut_Th    Sig_Th      Rapport")
    println("       :   Sig_Num")
    println("       :   Ecart_relatif  Incertitude")
    println()
    @printf("s+(pb) : %.5f | %.5f | %.6f\n", σ₊, 3σ, σ₊/(3σ))
    @printf("       : %.5f\n", mc₊)
    @printf("       : %.6f | %.8f | %.2f\n", mc₊/σ₊-1, incr₊, (mc₊/σ₊-1)/incr₊)
    println()
    @printf("s-(pb) : %.5f | %.4f | %.6f\n", σ₋, 5σ, σ₋/(5σ))
    @printf("       : %.5f\n", mc₋)
    @printf("       : %.6f | %.9f | %.2f\n", mc₋/σ₋-1, incr₋, (mc₋/σ₋-1)/incr₋)
    println()
end

end