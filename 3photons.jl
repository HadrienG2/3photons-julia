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
include("Random.jl")      # Used, depends on: Errors.jl, Numeric.jl
include("EvGen.jl")       # Used, depends on: Errors.jl, LinAlg.jl, Numeric.jl,
                          #                   Random.jl
include("Spinor.jl")      # Depends on: Errors.jl, EvGen.jl, LinAlg.jl,
                          #             Numeric.jl
include("EvCut.jl")       # Used, depends on: Errors.jl, EvGen.jl, LinAlg.jl,
                          #                   Numeric.jl
include("Config.jl")      # Used, depends on: Errors.jl, EvCut.jl, Numeric.jl
include("Coupling.jl")    # Used, depends on: Config.jl, Numeric.jl
include("ResCont.jl")     # Used, depends on: Coupling.jl, Errors.jl, EvGen.jl,
                          #                   Numeric.jl
include("ResFin.jl")      # Used, depends on: Config.jl, Errors.jl, EvGen.jl,
                          #                   Numeric.jl, ResCont.jl
include("Output.jl")      # Used, depends on: Config.jl, Numeric.jl, ResFin.jl
include("Scheduling.jl")  # Used, depends on: Errors.jl, Random.jl, ResFin.jl


"Artificial module introduced as a performance optimization"
module MainModule

using ..Config: Configuration
using ..Coupling: Couplings
using ..EvCut: keep_event
using ..EvGen: EventGenerator, generate_event!
using ..Output: dump_results
using ..Random: RandomGenerator
using ..ResCont: ResultContribution
using ..ResFin: integrate_contrib!, ResultsBuilder
using ..Scheduling: run_simulation

export main


# TODO: After translating, turns this into more idiomatic Julia (e.g. more
#       fancy unicode variable names and broadcasting, more genericity...)

"Artificial function introduced as a performance optimization"
function main(;jit_warmup::Bool=false)
    # === CONFIGURATION READOUT ===

    # Load the configuration from its file
    cfg = Configuration("valeurs"; jit_warmup=jit_warmup)


    # === SIMULATION INITIALIZATION ===

    # Record when the simulation started
    #
    # NOTE: Unlike the C++ version, we do this after loading the configuration
    #       file, which reduces IO-induced timing fluctuations.
    #
    start_time_s = time()

    # NOTE: Removed final particle mass array. Since we are simulating photons,
    #       we know that all masses will be zero.

    # NOTE: Deleted the original WTEV value. In the C++ code, it was overwritten
    #       by the first RAMBO call w/o having ever been read!

    # Compute physical couplings
    couplings = Couplings(cfg)

    # Initialize the event generator
    evgen = EventGenerator(cfg.e_tot; jit_warmup=jit_warmup)


    # === SIMULATION EXECUTION ===

    """
    This kernel simulates a number of events, given an initial random number
    generator state, and return the accumulated intermediary results.
    """
    function simulate_events(num_events::UInt, rng::RandomGenerator)::ResultsBuilder
        # Setup a results accumulator
        res_builder = ResultsBuilder(cfg, evgen.event_weight)

        # Simulate the requested number of events
        #
        # DEBUG: If that seems stuck, try @time for allocations profiling...
        #
        for _ = 1:num_events
            # Generate an event
            event = generate_event!(rng, evgen)

            # If the event passes the cut, compute the total weight (incl.
            # matrix elements) and integrate it into the final results.
            if keep_event(cfg.event_cut, event)
                res_contrib = ResultContribution(couplings, event)
                # NOTE: The original code would display the result here
                integrate_contrib!(res_builder, res_contrib)
                # NOTE: The FORTRAN code would fill histograms here
            end
        end

        # Return the accumulated results
        res_builder
    end

    # Run the simulation
    num_events = jit_warmup ? UInt(1) : cfg.num_events
    result = run_simulation(num_events, simulate_events)

    # NOTE: This is where the FORTRAN code would normalize histograms


    # === RESULTS DISPLAY AND STORAGE ===

    # Measure how much time has elapsed
    elapsed_secs = time() - start_time_s

    # Send the results to the standard output and to disk
    if !jit_warmup
        dump_results(cfg, result, elapsed_secs)
    end
end

end


using Profile

# FIXME: The need for such ceremony to get optimal performance from the main
#        function is arguably a Julia bug...
MainModule.main(jit_warmup=true)  # Get JIT out of our profile
@profile @time MainModule.main()
Profile.print(mincount=130, noisefloor=1.0)  # Ignore <2% + low wrt parent