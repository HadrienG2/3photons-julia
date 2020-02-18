# FIXME: This manual dependency tracking is incredibly ugly, and I'm very
#        surprised that I apparently need to engage in it. Investigate if Julia
#        truly doesn't provide any better way to organize source code.
include("errors.jl")    # No dependency
include("numeric.jl")   # No dependency
include("evcut.jl")     # Depends on: numeric.jl
include("rescont.jl")   # Depends on: numeric.jl
include("config.jl")    # Used, depends on: errors.jl, evcut.jl, numeric.jl
include("coupling.jl")  # Used, depends on: config.jl, numeric.jl
include("event.jl")     # Used, depends on: errors.jl, numeric.jl
include("random.jl")    # Used, no dependency
include("resfin.jl")    # Used, depends on: config.jl, event.jl, numeric.jl,
                        #                   rescont.jl


using .Config: Configuration
using .Coupling: Couplings
using .Event: EventGenerator
using .Random: RandomGenerator
using .ResFin: ResultsBuilder


# === CONFIGURATION READOUT ===

# Load the configuration from its file
cfg = Configuration("valeurs")


# === SIMULATION INITIALIZATION ===

# Record when the simulation started
#
# NOTE: Unlike the C++ version, we do this after loading the configuration file,
#       which reduces IO-induced timing fluctuations.
#
start_time_s = time()

# NOTE: Removed final particle mass array. Since we are simulating photons, we
#       know that all masses will be zero.

# NOTE: Deleted the original WTEV value. In the C++ code, it was overwritten by
#       the first RAMBO call w/o having ever been read!

# Compute physical couplings
couplings = Couplings(cfg)

# Initialize the event generator
evgen = EventGenerator(cfg.e_tot)


# === SIMULATION EXECUTION ===

"""
This kernel simulates a number of events, given an initial random number
generator state, and return the accumulated intermediary results.
"""
function simulate_events(num_events::UInt, rng::RandomGenerator)::ResultsBuilder
    # Setup a results accumulator
    res_builder = ResultsBuilder(cfg, evgen.event_weight)

    # Simulate the requested number of events
    for _ = 1:num_events
        # TODO: Not implemented yet
        throw(AssertionError("Not implemented yet"))
    end

    # Return the accumulated results
    res_builder
end

# TODO: Replace with actual execution once we have that
simulate_events(cfg.num_events, RandomGenerator())


# TODO: Finish translating the program
#
# TODO: After translating, turns this into more idiomatic Julia (e.g. unicode
#       variable names, more genericity...
throw(AssertionError("Not implemented yet"))