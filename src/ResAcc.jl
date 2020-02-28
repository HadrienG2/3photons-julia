# Depends on Config.jl, Errors.jl, EvData.jl, MatElems.jl and Numeric.jl being
# include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Integration of simulation results across generated events"
module ResAcc

using ..Config: Configuration
using ..Errors: @enforce
using ..EvData: NUM_OUTGOING
using ..MatElems: m²_sums, MEsContributions, MEsVector, NUM_MAT_ELEMS
using ..Numeric: Float
using LinearAlgebra: ⋅

export integrate_event!, merge_results!, ResultsAccumulator


"""
This struct will accumulate intermediary results during integration, and
ultimately be used to compute the final results (see ResFin module).
"""
mutable struct ResultsAccumulator
    # === RESULT ACCUMULATORS ===

    "Number of integrated events"
    selected_events::UInt

    "Accumulated cross-section for each contribution"
    spm²::MEsVector

    "Accumulated variance for each contribution"
    vars::MEsVector

    "Impact of each contribution on the cross-section"
    σ_contribs::MEsVector

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
function ResultsAccumulator(cfg::Configuration, event_weight::Float)
    # This code depends on some aspects of the problem definition
    @enforce (NUM_MAT_ELEMS == 5) "This code assumes 5 matrix elements"

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
    σ_contribs = MEsVector(
        aa_contrib,              # A
        bb_contrib * cfg.𝛽₊^2,   # B₊
        bb_contrib * cfg.𝛽₋^2,   # B₋
        ab_contrib * ecart_pic,  # R_MX
        -ab_contrib,             # I_MX
    )

    # Return a complete results accumulator
    #
    # FIXME: Isn't there any way to say which field we are talking about?
    #
    ResultsAccumulator(
        0,                 # selected_events
        zeros(MEsVector),  # spm²
        zeros(MEsVector),  # vars
        σ_contribs,
        0,                 # σ
        0,                 # variance

        cfg,
        fact_com,
        norm_weight,
        propag,
        ecart_pic,
    )
end


"Integrate an event's contribution into the simulation results"
function integrate_event!(acc::ResultsAccumulator, result::MEsContributions)
    acc.selected_events += 1
    spm²_dif = m²_sums(result)
    acc.spm² += spm²_dif
    acc.vars += spm²_dif.^2
    weight = spm²_dif ⋅ acc.σ_contribs
    acc.σ += weight
    acc.variance += weight^2
end


"Integrate pre-aggregated results from another ResultsAccumulator"
function merge_results!(dest::ResultsAccumulator, src::ResultsAccumulator)
    dest.selected_events += src.selected_events
    dest.spm² += src.spm²
    dest.vars += src.vars
    dest.σ += src.σ
    dest.variance += src.variance
end

end