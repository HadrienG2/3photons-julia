# FIXME: This manual dependency tracking is incredibly ugly, and I'm very
#        surprised that I apparently need to engage in it. Investigate if Julia
#        truly doesn't provide any better way to organize source code.
include("errors.jl")    # No dependency
include("numeric.jl")   # No dependency
include("evcut.jl")     # Depends on: numeric.jl
include("config.jl")    # Used, depends on: errors.jl, evcut.jl, numeric.jl
include("event.jl")     # Used, depends on: errors.jl, numeric.jl
include("coupling.jl")  # Used, depends on: config.jl, numeric.jl


using .Config: Configuration
using .Coupling: Couplings
using .Event: EventGenerator


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

# TODO: Finish translating the program
#
# TODO: After translating, turns this into more idiomatic Julia (e.g. unicode
#       variable names, more genericity...