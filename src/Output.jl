# Depends on Config.jl, Errors.jl, EvData.jl, MatElems.jl, Numeric.jl and
# ResFin.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"""
This module is in charge of outputting the final simulation results to the
standard output and various files
"""
module Output

import Dates

using ..Config: Configuration
using ..Errors: @enforce
using ..EvData: NUM_SPINS
using ..MatElems: A, B₊, B₋, NUM_MAT_ELEMS, R_MX
using ..Numeric: Float
using ..ResFin: FinalResults, print_eric, print_fawzi
using LinearAlgebra: norm
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

    # Facilities for replicating 3photons' output styling
    writeln_3p(file) = write(file, "\n")
    writeln_3p(file, str) = write(file, " $str\n")
    label_3p(key) = @sprintf("%-31s: ", key)
    format_3p(val) = string(val)
    # FIXME: Should honor sig_digits here, but @sprintf only supports static
    #        format string and doesn't support .* for externally controlled
    #        precision. Tried Formatting, but it doesn't even support %g...
    format_3p(val::AbstractFloat) = @sprintf("%.14g", val)
    writeln_3p(file, key, val) = writeln_3p(file, label_3p(key)*format_3p(val))

    # Write execution timings to a file
    open("res.times", "w") do tim_file
        # Write a timestamp of when the run ended
        writeln_3p(tim_file, timestamp)

        # Write program performance stats
        writeln_3p(tim_file, "---------------------------------------------")
        writeln_3p(tim_file, "Temps ecoule", "???")
        writeln_3p(tim_file, "Temps ecoule utilisateur", elapsed_secs)
        writeln_3p(tim_file, "Temps ecoule systeme", "???")
        secs_per_ev = elapsed_secs / cfg.num_events
        writeln_3p(tim_file, "Temps ecoule par evenement", secs_per_ev)
    end

    # Write main results file. Try to mimick the original C++ format as well as
    # possible to ease comparisons, even where it makes little sense.
    open("res.data", "w") do dat_file
        # Shorthand
        ev_cut = cfg.event_cut

        # Write the results to the file
        writeln_3p(dat_file, "Nombre d'evenements", cfg.num_events)
        writeln_3p(dat_file, "... apres coupure", res.selected_events)
        writeln_3p(dat_file, "energie dans le CdM      (GeV)", cfg.e_tot)
        writeln_3p(dat_file, "coupure / cos(photon,faisceau)", ev_cut.a_cut)
        writeln_3p(dat_file, "coupure / cos(photon,photon)", ev_cut.b_cut)
        writeln_3p(dat_file, "coupure / sin(normale,faisceau)", ev_cut.sin_cut)
        writeln_3p(dat_file, "coupure sur l'energie    (GeV)", ev_cut.e_min)
        writeln_3p(dat_file, "1/(constante de structure fine)", 1/cfg.𝛼)
        writeln_3p(dat_file, "1/(structure fine au pic)", 1/cfg.𝛼_Z)
        writeln_3p(dat_file, "facteur de conversion GeV-2/pb", cfg.convers)
        writeln_3p(dat_file, "Masse du Z0              (GeV)", cfg.m_Z⁰)
        writeln_3p(dat_file, "Largeur du Z0            (GeV)", cfg.g_Z⁰)
        writeln_3p(dat_file, "Sinus^2 Theta Weinberg", cfg.sin²_w)
        writeln_3p(dat_file, "Taux de branchement Z--->e+e-", cfg.br_e₊_e₋)
        writeln_3p(dat_file, "Beta plus", cfg.𝛽₊)
        writeln_3p(dat_file, "Beta moins", cfg.𝛽₋)
        writeln_3p(dat_file, "---------------------------------------------")
        writeln_3p(dat_file, "Section Efficace          (pb)", res.σ)
        writeln_3p(dat_file, "Ecart-Type                (pb)", res.σ*res.prec)
        # Work around opinion divergence between Rust and Julia's %g logic
        # FIXME: Should honor sig_digits precision here, but see above.
        prec_str = @sprintf("%.13e", res.prec)
        writeln_3p(dat_file, "Precision Relative", prec_str)
        writeln_3p(dat_file, "---------------------------------------------")
        writeln_3p(dat_file, "Beta minimum", res.𝛽_min)
        writeln_3p(dat_file, "Stat. Significance  B+(pb-1/2)", res.ss₊)
        incert_ss₊ = res.ss₊*res.inc_ss₊
        writeln_3p(dat_file, "Incert. Stat. Sign. B+(pb-1/2)", incert_ss₊)
        writeln_3p(dat_file, "Stat. Significance  B-(pb-1/2)", res.ss₋)
        incert_ss₋ = res.ss₋*res.inc_ss₋
        writeln_3p(dat_file, "Incert. Stat. Sign. B-(pb-1/2)", incert_ss₋)

        # Write more results (nature and purpose unclear in C++ code...)
        writeln_3p(dat_file)
        decimals = min(sig_digits-1, 7)
        for sp=1:NUM_SPINS
            for elem=1:NUM_MAT_ELEMS
                writeln_3p(
                    dat_file,
                    # FIXME: Should honor decimals precision here, but see above
                    @sprintf("%2d%3d%15.7e%15.7e%15.7e",
                             sp,
                             elem,
                             res.spm²[sp, elem],
                             abs(res.spm²[sp, elem]) * res.vars[sp, elem],
                             res.vars[sp, elem])
                )
            end
            writeln_3p(dat_file)
        end
        for elem=1:NUM_MAT_ELEMS
            tmp₁ = sum(res.spm²[:, elem])
            tmp₂ = norm(res.spm²[:, elem] .* res.vars[:, elem])
            writeln_3p(
                dat_file,
                # FIXME: Should honor decimals precision here, but see above
                @sprintf("  %3d%15.7e%15.7e%15.7e",
                         elem,
                         tmp₁/4,
                         tmp₂/4,
                         tmp₂/abs(tmp₁))
            )
        end
    end
    
    # Append the results of this run to a cumulative file
    #
    # NOTE: This part is completely broken in the C++ version, I did my best to
    #       fix it in this version.
    open("pil.mc", "a") do cum_dat_file
        @enforce (NUM_MAT_ELEMS == 5) "This code assumes 5 matrix elements"

        println(cum_dat_file, timestamp)
        res₁ = sum(res.spm²[:, A])
        res₂ = sum(res.spm²[:, B₊]) * cfg.𝛽₊^2
        res₃ = sum(res.spm²[:, B₋]) * cfg.𝛽₋^2
        res₄ = sum(res.spm²[:, R_MX]) * cfg.𝛽₊
        avg = (res₁ + res₂ + res₃ + res₄) / 4
        println(cum_dat_file, "$(cfg.e_tot) $(res₁/4) $(res₂/4) $(res₃/4) "*
                              "$(res₄/4) $avg $(res.σ)")
    end
end

end