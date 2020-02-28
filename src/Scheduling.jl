# Depends on Errors.jl, Random.jl, ResAcc.jl and ResFin.jl being include-d
# beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"""
This module takes care of scheduling the simulation work, encapsulating use of
multiple threads and anything else that will come in the future
"""
module Scheduling

using ..Errors: @enforce
using ..Random: RandomGenerator
using ..ResAcc: merge_results!, ResultsAccumulator
using ..ResFin: FinalResults

export run_simulation


"""
Size of the simulated event batches

Simulated events are grouped in batches of a certain size in order to reduce
accumulation error and achieve perfect reproducibility between sequential
and parallel runs of the simulation.

This constant may need to be tuned in the future if CPUs become faster or
synchronization overhead changes. But the rate of such change is expected to
be low enough for hard-coding of this constant to be reasonable.
"""
const EVENT_BATCH_SIZE = UInt(10_000)


# FIXME: Extract sequential back-end in a separate file
#
# FIXME: Figure out how to cleanly express the Julia equivalent of an
#        `impl Send + Sync + Fn(usize, RandomGenerator) -> ResultsAccumulator`
#        static typing contract on simulate_events.
#
"""
Simulate events in sequential mode

We use batched logic even in sequential mode, in order to achieve
reproducibility with respect to multi-threaded runs.

Note that this is anyways generally a good thing to do when accumulating
lots of results, as otherwise the accumulator will eventually grow much
larger than the accumulated values and numerical accumulation errors
will start to blow up.
"""
function run_simulation_seq(num_events::UInt,
                            rng::RandomGenerator,
                            simulate_events::Function)::ResultsAccumulator
    # Some double-checking cannot hurt...
    @enforce (num_events > 0) "Must simulate at least one event"

    # Initialize the accumulator with the first batch of events
    first_batch_size = min(EVENT_BATCH_SIZE, num_events)
    num_events -= first_batch_size
    res_acc = simulate_events(first_batch_size, rng)

    # Simulate and integrate complete batches of events (if any)
    num_full_batches = num_events / EVENT_BATCH_SIZE
    for _ = 1:num_full_batches
        merge_results!(res_acc, simulate_events(EVENT_BATCH_SIZE, rng))
    end
    num_events %= EVENT_BATCH_SIZE

    # Integrate the remaining events
    merge_results!(res_acc, simulate_events(num_events, rng))

    # Return the final accumulated results
    res_acc
end


# FIXME: Figure out how to cleanly express the Julia equivalent of an
#        `impl Send + Sync + Fn(usize, RandomGenerator) -> ResultsAccumulator`
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
    # Check that the user is being reasonable (should have already been checked
    # at configuration time, but bugs can happen...)
    @enforce (num_events > 0) "Must simulate at least one event"

    # Initialize the random number generator
    rng = RandomGenerator()

    # Integrate simulation results in sequential mode
    #
    # FIXME: Port multi-threaded simulation schedule from Rust version and
    #        expose support for it.
    #
    res_acc = run_simulation_seq(num_events, rng, simulate_events)

    # Finalize the results
    FinalResults(res_acc)
end

end