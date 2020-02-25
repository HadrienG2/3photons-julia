# Depends on Config.jl, Errors.jl, EvGen.jl, Numeric.jl and ResCont.jl being
# include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"""
This module contains everything that is needed to compute, store, and analyze
the final results: differential cross-section, sum & variance.
"""
module ResFin

using ..Config: Configuration
using ..Errors: @enforce
using ..EvGen: NUM_OUTGOING, NUM_SPINS
using ..Numeric: Float
using ..ResCont: A, B₊, B₋, I_MX, m²_sums, NUM_RESULTS, ResultContribution,
                 ResultVector, R_MX
using LinearAlgebra: ⋅
using StaticArrays: MMatrix, SMatrix, @MMatrix, @SMatrix

export integrate_contrib!, finalize_results, merge_results!, ResultsBuilder,
       SP₋, SP₊


# === RESULTS ACCUMULATION ===

"""
This struct will accumulate intermediary results during integration, and
ultimately compute the final results (see FinalResults below).
"""
mutable struct ResultsBuilder
    # === RESULT ACCUMULATORS ===

    "Number of integrated events"
    selected_events::UInt

    "Accumulated cross-section for each contribution"
    spm²::ResultVector

    "Accumulated variance for each contribution"
    vars::ResultVector

    "Impact of each contribution on the cross-section"
    σ_contribs::ResultVector

    "Accumulated total cross-section"
    σ::Float

    "Accumulated total variance"
    variance::Float

    # === PHYSICAL CONSTANTS (CACHED FOR FINALIZATION) ===

    "Configuration of the simulation"
    cfg::Configuration

    """
    Common factor, non-averaged over spins=1
                     /(symmetry factor)
                     *(conversion factor GeV^-2->pb)
    To average over spins, add :
                     /(number of incoming helicities)
    """
    fact_com::Float

    "Event weight, with total phase space normalization"
    norm_weight::Float

    "Z° propagator"
    propag::Float

    "??? (Ask Vincent Lafage)"
    ecart_pic::Float
end


"Prepare for results integration"
function ResultsBuilder(cfg::Configuration, event_weight::Float)
    # Common factor (see definition and remarks above)
    fact_com = 1 / 6 * cfg.convers
    gzr = cfg.g_Z⁰ / cfg.m_Z⁰

    # Sum over polarisations factors
    p_aa = 2.
    p_ab = 1 - 4 * cfg.sin²_w
    p_bb = p_ab + 8 * cfg.sin²_w^2

    # Homogeneity coefficient
    c_aa = fact_com * p_aa
    c_ab = fact_com * p_ab / cfg.m_Z⁰^2
    c_bb = fact_com * p_bb / cfg.m_Z⁰^4

    # Switch to dimensionless variable
    ζ = (cfg.e_tot / cfg.m_Z⁰)^2
    ecart_pic = (ζ - 1) / gzr
    propag = 1 / (1 + ecart_pic^2)

    # Apply total phase space normalization to the event weight
    norm = (2π)^(4 - 3*NUM_OUTGOING) / cfg.num_events
    # NOTE: This replaces the original WTEV, previously reset every event
    norm_weight = event_weight * norm

    # Compute how much each result contribution adds to the cross-section.
    # Again, this avoids duplicate work in the integration loop.
    com_contrib = norm_weight / 4
    aa_contrib = com_contrib * c_aa
    bb_contrib = com_contrib * c_bb * propag / gzr^2
    ab_contrib = com_contrib * c_ab * 2 * cfg.𝛽₊ * propag / gzr
    σ_contribs = ResultVector(
        aa_contrib,
        bb_contrib * cfg.𝛽₊^2,
        bb_contrib * cfg.𝛽₋^2,
        ab_contrib * ecart_pic,
        -ab_contrib,
    )

    # Return a complete results builder
    #
    # FIXME: Isn't there any way to say which field we are talking about?
    #
    ResultsBuilder(
        0,                   # selected_events
        zeros(NUM_RESULTS),  # spm²
        zeros(NUM_RESULTS),  # vars
        σ_contribs,
        0,                   # σ
        0,                   # variance

        cfg,
        fact_com,
        norm_weight,
        propag,
        ecart_pic,
    )
end


"Integrate one intermediary result into the simulation results"
function integrate_contrib!(builder::ResultsBuilder, result::ResultContribution)
    builder.selected_events += 1
    spm²_dif = m²_sums(result)
    builder.spm² += spm²_dif
    builder.vars += spm²_dif.^2
    weight = spm²_dif ⋅ builder.σ_contribs
    builder.σ += weight
    builder.variance += weight^2
end


"Integrate pre-aggregated results from another ResultsBuilder"
function merge_results!(dest::ResultsBuilder, src::ResultsBuilder)
    dest.selected_events += src.selected_events
    dest.spm² += src.spm²
    dest.vars += src.vars
    dest.σ += src.σ
    dest.variance += src.variance
end


# === FINAL RESULTS ===

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
function finalize_results(builder::ResultsBuilder)::FinalResults
    # This code depends on some aspects of the problem definition
    @enforce (NUM_SPINS == 2) "This code currently assumes 2 spins"
    @enforce (NUM_RESULTS == 5) "This code currently assumes 5 matrix elements"

    # Simulation configuration shorthands
    cfg = builder.cfg
    n_ev = cfg.num_events

    # Compute the relative uncertainties for one spin
    for (v_spm², v_var) ∈ zip(builder.spm², builder.vars)
        v_var = (v_var - v_spm²^2 / n_ev) / (n_ev - 1)
        v_var = √(v_var / n_ev) / abs(v_spm² / n_ev)
    end

    # Copy for the opposite spin
    spm² = @MMatrix [
        builder.spm²[res]
        for _spin=1:NUM_SPINS, res=1:NUM_RESULTS
    ]
    vars = @SMatrix [
        builder.vars[res]
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
    spm² *= builder.fact_com * flux * builder.norm_weight
    gm_Z⁰ = cfg.g_Z⁰ * cfg.m_Z⁰
    spm²[:, B₊:I_MX] *= builder.propag / gm_Z⁰
    spm²[:, B₊:B₋] /= gm_Z⁰
    spm²[:, R_MX] *= builder.ecart_pic

    # Compute other parts of the result
    𝛽_min = √((spm²[SP₋, A] + spm²[SP₊, A]) / (spm²[SP₋, B₊] + spm²[SP₊, B₊]))

    ss_denom = spm²[SP₋, A] + spm²[SP₊, A]
    ss_norm = 1 / (2 + √ss_denom)

    ss₊ = (spm²[SP₋, B₊] + spm²[SP₊, B₊]) * ss_norm
    ss₋ = (spm²[SP₋, B₋] + spm²[SP₊, B₋]) * ss_norm

    inc_ss_common =
        √((spm²[SP₋, A] * vars[SP₋, A])^2 + (spm²[SP₊, A] * vars[SP₊, A])^2) /
            (2 * abs(ss_denom))

    inc_ss₊ =
        √((spm²[SP₋, B₊] * vars[SP₋, B₊])^2 + (spm²[SP₊, B₊] * vars[SP₊, B₊])^2) /
            abs(spm²[SP₋, B₊] + spm²[SP₊, B₊]) +
        inc_ss_common
    inc_ss₋ =
        √((spm²[SP₋, B₋] * vars[SP₋, B₋])^2 + (spm²[SP₊, B₋] * vars[SP₊, B₋])^2) /
            abs(spm²[SP₋, B₋] + spm²[SP₊, B₋]) +
        inc_ss_common

    variance = (builder.variance - builder.σ^2 / n_ev) / (n_ev - 1)
    prec = √(variance / n_ev) / abs(builder.σ / n_ev)
    σ = builder.σ * flux

    # Return the final results
    #
    # FIXME: Isn't there any way to say which field we are talking about?
    #
    FinalResults(
        builder.selected_events,
        spm²,
        vars,
        σ,
        variance,
        𝛽_min,
        prec,
        ss₊,
        inc_ss₊,
        ss₋,
        inc_ss₋,
        cfg,
    )
end

end