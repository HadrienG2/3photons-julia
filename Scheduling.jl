# No include dependency
#
# FIXME: Isn't there a way to spell this out in code???


"""
This module takes care of scheduling the simulation work, encapsulating use of
multiple threads and anything else that will come in the future
"""
module Scheduling

"""
Size of the simulated event batches

Simulated events are grouped in batches of a certain size in order to reduce
accumulation error and achieve perfect reproducibility between sequential
and parallel runs of the simulation.

This constant may need to be tuned in the future if CPUs become faster or
synchronization overhead changes. But the rate of such change is expected to
be low enough for hard-coding of this constant to be reasonable.
"""
const EVENT_BATCH_SIZE = 10000


# FIXME: Port multi-threaded simulation schedule from Rust version and expose
#        support for it somehow.
#
# FIXME: Figure out how to cleanly express the Julia equivalent of an
#        `impl Send + Sync + Fn(usize, RandomGenerator) -> ResultsBuilder`
#        static typing contract on simulate_events.
#
"""
Run the simulation in the manner that was configured.

Takes as parameters the total number of events to be simulated, and a
simulation kernel that simulates a certain number of events given an initial
random number generator state.

Returns the finalized simulation results
"""
function run_simulation(num_events::UInt, simulate_events::Function)::FinalResults
    # TODO: Finish translating the program
    throw(AssertionError("Not implemented yet"))
end

end