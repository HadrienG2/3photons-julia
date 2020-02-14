# This file depends on config.jl and numeric.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Physical couplings used for result computations"
module Coupling

using ..Config: Configuration
using ..Numeric: Float

export Couplings


"Set of physical couplings"
struct Couplings
    "Standard Model contribution electromagnetic coupling √(4𝜋𝛼)³"
    g_a::Float

    "𝛽₊ anomalous contribution electroweak coupling"
    g_bp::Float

    "𝛽₋ anomalous contribution electroweak coupling"
    g_bm::Float
end


"Fill in the parameters using data from the configuration file"
function Couplings(cfg::Configuration)
    e2 = 4 * π * cfg.alpha
    e2_z = 4 * π * cfg.alpha_z
    cos2_w = 1. - cfg.sin2_w
    g_beta = -sqrt(e2_z / (4 * cos2_w * cfg.sin2_w)) / cfg.m_z0 ^ 4

    # FIXME: Isn't there any way to say which field we are talking about?
    Couplings(
        -(sqrt(e2)^3),  # g_a
        g_beta,         # g_bp
        g_beta,         # g_bm
    )
end

end