# Depends on Config.jl, EvGen.jl, Numeric.jl and ResCont.jl being include-d
# beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"""
This module contains everything that is needed to compute, store, and analyze
the final results: differential cross-section, sum & variance.
"""
module ResFin

using ..Config: Configuration
using ..EvGen: OUTGOING_COUNT
using ..Numeric: Float
using ..ResCont: NUM_RESULTS, ResultVector

export ResultsBuilder


"""
This struct will accumulate intermediary results during integration, and
ultimately compute the final results (see FinalResults below).
"""
struct ResultsBuilder
    # === RESULT ACCUMULATORS ===

    "Number of integrated events"
    selected_events::UInt

    "Accumulated cross-section for each contribution"
    spm2::ResultVector

    "Accumulated variance for each contribution"
    vars::ResultVector

    "Impact of each contribution on the cross-section"
    sigma_contribs::ResultVector

    "Accumulated total cross-section"
    sigma::Float

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
    gzr = cfg.g_z0 / cfg.m_z0

    # Sum over polarisations factors
    p_aa = 2.
    p_ab = 1 - 4 * cfg.sin2_w
    p_bb = p_ab + 8 * cfg.sin2_w^2

    # Homogeneity coefficient
    c_aa = fact_com * p_aa
    c_ab = fact_com * p_ab / cfg.m_z0^2
    c_bb = fact_com * p_bb / cfg.m_z0^4

    # Switch to dimensionless variable
    dzeta = (cfg.e_tot / cfg.m_z0)^2
    ecart_pic = (dzeta - 1) / gzr
    propag = 1 / (1 + ecart_pic^2)

    # Apply total phase space normalization to the event weight
    norm = (2π)^(4 - 3*OUTGOING_COUNT) / cfg.num_events
    # NOTE: This replaces the original WTEV, previously reset every event
    norm_weight = event_weight * norm

    # Compute how much each result contribution adds to the cross-section.
    # Again, this avoids duplicate work in the integration loop.
    com_contrib = norm_weight / 4
    aa_contrib = com_contrib * c_aa
    bb_contrib = com_contrib * c_bb * propag / gzr^2
    ab_contrib = com_contrib * c_ab * 2 * cfg.beta_plus * propag / gzr
    sigma_contribs = ResultVector(
        aa_contrib,
        bb_contrib * cfg.beta_plus^2,
        bb_contrib * cfg.beta_minus^2,
        ab_contrib * ecart_pic,
        -ab_contrib,
    )

    # Return a complete results builder
    #
    # FIXME: Isn't there any way to say which field we are talking about?
    #
    ResultsBuilder(
        0,                   # selected_events
        zeros(NUM_RESULTS),  # spm2
        zeros(NUM_RESULTS),  # vars
        sigma_contribs,
        0,                   # sigma
        0,                   # variance

        cfg,
        fact_com,
        norm_weight,
        propag,
        ecart_pic,
    )
end

end