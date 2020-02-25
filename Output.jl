# Depends on Config.jl, Numeric.jl and ResFin.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"""
This module is in charge of outputting the final simulation results to the
standard output and various files
"""
module Output

using ..Config: Configuration
using ..Numeric: Float
using ..ResFin: FinalResults

export dump_results


"Output the simulation results to the console and to disk"
function dump_results(cfg::Configuration,
                      res_fin::FinalResults,
                      elapsed_secs::Float)
    # TODO: Finish translating the program
    throw(AssertionError("Not implemented yet"))
end

end