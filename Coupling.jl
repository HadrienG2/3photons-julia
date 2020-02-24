# Depends on Config.jl and Numeric.jl being include-d beforehand
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
    g_𝛼::Float

    "𝛽₊ anomalous contribution electroweak coupling"
    g_𝛽₊::Float

    "𝛽₋ anomalous contribution electroweak coupling"
    g_𝛽₋::Float
end


"Fill in the parameters using data from the configuration file"
function Couplings(cfg::Configuration)
    e² = 4π * cfg.𝛼
    e²_Z = 4π * cfg.𝛼_Z
    cos²_w = 1. - cfg.sin²_w
    g_𝛽 = -√(e²_Z / (4 * cos²_w * cfg.sin²_w)) / cfg.m_Z⁰^4

    # FIXME: Isn't there any way to say which field we are talking about?
    Couplings(
        -√e²^3,  # g_𝛼
        g_𝛽,     # g_𝛽₊
        g_𝛽,     # g_𝛽₋
    )
end

end