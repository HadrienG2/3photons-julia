# 3 photons: A simple Monte Carlo simulation
#
#
# # Introduction (for the physicist)
#
# This small computational program computes cross-section for the particle
# physics process electron + positron gives three photons (e⁺e⁻ → 𝛾𝛾𝛾).
#
# It distinguishes a classical Standard Model contribution, of purely Quantum
# ElectroDynamic origin and an hypothetic, beyond the Standard Model, New
# Physics contribution, phenomenologically described by two effective
# operators.
#
# It was designed in the LEP era, so these new interactions occurs between the
# Z⁰ boson and the three photons.
#
# The effective operator can be related to specific models, among which
# magnetic monopoles that run in a four points loop. The two operators exhibit
# different (???)
#
#
# # Introduction (for the numerical guy)
#
# The physicist want to compute a (multidimensional) integral, so we chose a
# Monte Carlo algorithm
#
#
# # Introduction (for the computer guy)
#
# this program started in a purely procedural style:
#
# * read in parameters and initialise counters
# * loop over (random) event,
#     * determining their geometrical and energy configuration,
#     * their phase space weight,
#     * their transition probability for each polarisation/helicity
#       configuration, depending on coupling strength
#     * sum it up
# * then display / store the result.
#
# The use of common (for the original Fortran) or struct (in C) or record
# types (in Ada) or classes (in C++) illustrates an object oriented design.
#
# The fact that we can plug each phase's output as the input of the next phase
# lend to a functionnal approach.

# FIXME: The above is really a program-wide doc, figure out if Julia supports
#        exposing that and if so how it should be expressed.


# FIXME: This manual dependency tracking is incredibly ugly, and I'm very
#        surprised that I apparently need to engage in it. Investigate if Julia
#        truly doesn't provide any better way to organize source code.
include("Errors.jl")      # No dependency
include("LinAlg.jl")      # No dependency
include("Numeric.jl")     # No dependency
include("EvData.jl")      # Depends on: LinAlg.jl, Numeric.jl
include("Random.jl")      # Depends on: Errors.jl, Numeric.jl
include("EvCut.jl")       # Depends on: Errors.jl, EvData.jl, LinAlg.jl,
                          #             Numeric.jl
include("EvGen.jl")       # Depends on: Errors.jl, EvData.jl, LinAlg.jl,
                          #             Numeric.jl, Random.jl
include("Spinor.jl")      # Depends on: Errors.jl, EvData.jl, LinAlg.jl,
                          #             Numeric.jl
include("Config.jl")      # Depends on: Errors.jl, EvCut.jl, Numeric.jl
include("Coupling.jl")    # Depends on: Config.jl, Numeric.jl
include("MatElems.jl")    # Depends on: Coupling.jl, Errors.jl, EvData.jl,
                          #             Numeric.jl
include("ResAcc.jl")      # Depends on: Config.jl, Errors.jl, EvData.jl,
                          #             Numeric.jl, MatElems.jl
include("ResFin.jl")      # Depends on: Config.jl, Errors.jl, EvData.jl,
                          #             Numeric.jl, ResAcc.jl, MatElems.jl
include("Output.jl")      # Depends on: Config.jl, EvData.jl, Numeric.jl,
                          #             MatElems.jl, ResFin.jl
include("Scheduling.jl")  # Depends on: Errors.jl, Random.jl, ResAcc.jl,
                          #             ResFin.jl


"Artificial module introduced as a performance optimization"
module MainModule

using ..Config: Configuration
using ..Coupling: Couplings
using ..EvCut: keep_event
using ..EvGen: EventGenerator, generate_event!
using ..MatElems: MEsContributions
using ..Output: dump_results
using ..Random: RandomGenerator
using ..ResAcc: integrate_event!, ResultsAccumulator
using ..Scheduling: run_simulation

export main


# FIXME: How to specify that NO_PHOTON_SORTING should be a Bool?
"Artificial function introduced as a performance optimization"
function main(
    no_photon_sorting::Val{NO_PHOTON_SORTING};
    jit_warmup::Bool=false,
) where NO_PHOTON_SORTING
    # === CONFIGURATION READOUT ===

    # Load the configuration from its file
    # TODO: Optimize this after optimizing simulation body
    cfg = try
        Configuration("valeurs"; jit_warmup=jit_warmup)
    catch e
        println("Failed to load the configuration")
        rethrow(e)
    end

    # === SIMULATION INITIALIZATION ===

    # NOTE: Unlike the C++ version, we start the clock after configuration I/O,
    #       to avoid IO-induced timing fluctuations.
    start_time_s = time()

    # NOTE: Removed final particle mass array. Since we are simulating photons,
    #       we know that all masses will be zero.

    # NOTE: Deleted the original WTEV value. In the C++ code, it was overwritten
    #       by the first RAMBO call w/o having ever been read!

    # Compute physical couplings
    couplings = Couplings(cfg)

    # Initialize the event generator
    evgen = EventGenerator(cfg.e_tot,
                           no_photon_sorting;
                           jit_warmup=jit_warmup)

    # === SIMULATION EXECUTION ===

    """
    This kernel simulates a number of events, given an initial random number
    generator state, and return the accumulated intermediary results.
    """
    function simulate_events(num_events::UInt,
                             rng::RandomGenerator)::ResultsAccumulator
        # Setup a results accumulator
        res_acc = ResultsAccumulator(cfg, evgen.event_weight)

        # Simulate the requested number of events
        for _ = 1:num_events
            # Generate an event
            event = generate_event!(rng, evgen)

            # If the event passes the cut...
            if keep_event(cfg.event_cut, event)
                # Compute the total weight, including matrix elements
                me_contribs = MEsContributions(couplings, event)

                # NOTE: The original code would display the result here

                # Integrate the event's contribution into the results
                integrate_event!(res_acc, me_contribs)

                # NOTE: The FORTRAN code would fill histograms here
            end
        end

        # Return the accumulated results
        res_acc
    end

    # Run the simulation
    # TODO: Optimize this first
    num_events = jit_warmup ? UInt(1) : cfg.num_events
    println("DEBUG: Running simulation body... ")
    @time result = run_simulation(num_events, simulate_events)

    # NOTE: This is where the FORTRAN code would normalize histograms

    # === RESULTS DISPLAY AND STORAGE ===

    # Measure how much time has elapsed
    elapsed_secs = time() - start_time_s

    # Send the results to the standard output and to disk
    # TODO: Optimize this after optimizing simulation body
    if !jit_warmup
        try
            dump_results(cfg, result, elapsed_secs)
        catch e
            println("Failed to output the results")
            rethrow(e)
        end
    end
end

end


using Profile

# Collect "compilation flags" from CLI args
no_photon_sorting = Val("--no-photon-sorting" in ARGS)

# FIXME: The need for such ceremony to get optimal performance from the main
#        function and non-noisy profiles is arguably a Julia issue...
MainModule.main(no_photon_sorting; jit_warmup=true)
Profile.clear_malloc_data()
# DEBUG: Hide this while dumping code, bring it back after analysis
@profile @time MainModule.main(no_photon_sorting)
Profile.print(mincount=120, noisefloor=1.0)  # Ignore <2%run_sim + low vs parent