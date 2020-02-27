# Depends on Config.jl, Numeric.jl, ResCont.jl and ResFin.jl being include-d
# beforehand
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
using ..ResCont: NUM_RESULTS
using ..ResFin: FinalResults, NUM_SPINS, print_eric, print_fawzi
using Printf: @sprintf

export dump_results


"Output the simulation results to the console and to disk"
function dump_results(cfg::Configuration,
                      res::FinalResults,
                      elapsed_secs::Float)
    # Print out some final results on stdout
    print_eric(res)
    print_fawzi(res)

    # Compute a timestamp of when the run ended
    current_time = Dates.now()
    timestamp = Dates.format(current_time, "dd-uuu-yy   HH:MM:SS")

    # Number of significant digits in file output
    #
    # Must print one less than the actual machine type precision to match the
    # output of the C++ version of 3photons.
    #
    sig_digits = floor(Int, -log10(eps(Float))) - 1

    # So, apparently, Julia has println but not writeln because... reasons?
    # Well, that's a good occasion to enforce 3photon-style formatting anyway.
    writeln(file) = write(file, "\n")
    writeln(file, str) = write(file, " $(str)\n")
    label(key) = @sprintf("%-31s: ", key)
    format(value) = string(value)
    # FIXME: Should honor sig_digits here, but @sprintf only supports static
    #        format string and doesn't support .* for externally controlled
    #        precision. Tried Formatting, but it doesn't even support %g...
    format(value::AbstractFloat) = @sprintf("%.14g", value)
    writeln(file, key, value) = writeln(file, label(key)*format(value))

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
        # Shorthand
        ev_cut = cfg.event_cut

        # Write the results to the file
        writeln(dat_file, "Nombre d'evenements", cfg.num_events)
        writeln(dat_file, "... apres coupure", res.selected_events)
        writeln(dat_file, "energie dans le CdM      (GeV)", cfg.e_tot)
        writeln(dat_file, "coupure / cos(photon,faisceau)", ev_cut.a_cut)
        writeln(dat_file, "coupure / cos(photon,photon)", ev_cut.b_cut)
        writeln(dat_file, "coupure / sin(normale,faisceau)", ev_cut.sin_cut)
        writeln(dat_file, "coupure sur l'energie    (GeV)", ev_cut.e_min)
        writeln(dat_file, "1/(constante de structure fine)", 1/cfg.𝛼)
        writeln(dat_file, "1/(structure fine au pic)", 1/cfg.𝛼_Z)
        writeln(dat_file, "facteur de conversion GeV-2/pb", cfg.convers)
        writeln(dat_file, "Masse du Z0              (GeV)", cfg.m_Z⁰)
        writeln(dat_file, "Largeur du Z0            (GeV)", cfg.g_Z⁰)
        writeln(dat_file, "Sinus^2 Theta Weinberg", cfg.sin²_w)
        writeln(dat_file, "Taux de branchement Z--->e+e-", cfg.br_e₊_e₋)
        writeln(dat_file, "Beta plus", cfg.𝛽₊)
        writeln(dat_file, "Beta moins", cfg.𝛽₋)
        writeln(dat_file, "---------------------------------------------")
        writeln(dat_file, "Section Efficace          (pb)", res.σ)
        writeln(dat_file, "Ecart-Type                (pb)", res.σ*res.prec)
        # Work around opinion divergence between Rust and Julia's %g logic
        # FIXME: Should honor sig_digits precision here, but see above.
        prec_str = @sprintf("%.13e", res.prec)
        writeln(dat_file, "Precision Relative", prec_str)
        writeln(dat_file, "---------------------------------------------")
        writeln(dat_file, "Beta minimum", res.𝛽_min)
        writeln(dat_file, "Stat. Significance  B+(pb-1/2)", res.ss₊)
        writeln(dat_file, "Incert. Stat. Sign. B+(pb-1/2)", res.ss₊*res.inc_ss₊)
        writeln(dat_file, "Stat. Significance  B-(pb-1/2)", res.ss₋)
        writeln(dat_file, "Incert. Stat. Sign. B-(pb-1/2)", res.ss₋*res.inc_ss₋)

        # Write more results (nature and purpose unclear in C++ code...)
        writeln(dat_file)
        decimals = min(sig_digits-1, 7)
        for sp=1:NUM_SPINS
            for k=1:NUM_RESULTS
                writeln(
                    dat_file,
                    # FIXME: Should honor decimals precision here, but see above
                    @sprintf("%2d%3d%15.7e%15.7e%15.7e",
                             sp,
                             k,
                             res.spm²[sp, k],
                             abs(res.spm²[sp, k]) * res.vars[sp, k],
                             res.vars[sp, k])
                )
            end
            writeln(dat_file)
        end
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