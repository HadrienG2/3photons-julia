# Depends on Config.jl, Numeric.jl and ResFin.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"""
This module is in charge of outputting the final simulation results to the
standard output and various files
"""
module Output

import Dates

using ..Config: Configuration
using ..Numeric: Float
using ..ResFin: FinalResults, print_eric, print_fawzi
using Printf: @sprintf

export dump_results


"Output the simulation results to the console and to disk"
function dump_results(cfg::Configuration,
                      res_fin::FinalResults,
                      elapsed_secs::Float)
    # Print out some final results on stdout
    print_eric(res_fin)
    print_fawzi(res_fin)

    # Compute a timestamp of when the run ended
    current_time = Dates.now()
    timestamp = Dates.format(current_time, "dd-uuu-yy   HH:MM:SS")

    # So, apparently, Julia has println but not writeln because... reasons?
    # Well, that's a good occasion to enforce 3photon-style formatting anyway.
    writeln(file, str) = write(file, " $(str)\n")
    label(key) = @sprintf("%-31s: ", key)
    writeln(file, key, value) = writeln(file, label(key)*string(value))

    # Write execution timings to a file
    open("res.times", "w") do tim_file
        # Write a timestamp of when the run ended
        writeln(tim_file, timestamp)

        # Write program performance stats
        writeln(tim_file, "---------------------------------------------")
        writeln(tim_file, "Temps ecoule", "???")
        writeln(tim_file, "Temps ecoule utilisateur", elapsed_secs)
        writeln(tim_file, "Temps ecoule systeme", "???")
        secs_per_ev = elapsed_secs / cfg.num_events
        writeln(tim_file, "Temps ecoule par evenement", secs_per_ev)
    end

    # Write main results file. Try to mimick the original C++ format as well as
    # possible to ease comparisons, even where it makes little sense.
    open("res.data", "w") do dat_file
        # TODO: Finish translating the program
        throw(AssertionError("Not implemented yet"))
    end
    
    # Append the results of this run to a cumulative file
    #
    # NOTE: This part is completely broken in the C++ version, I did my best to
    #       fix it in this version.
    open("pil.mc", "a") do cum_dat_file
        # TODO: Finish translating the program
        throw(AssertionError("Not implemented yet"))
    end
end

end