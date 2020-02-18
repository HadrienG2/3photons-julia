# No include dependency
#
# FIXME: Isn't there a way to spell this out in code???


"""
Random number generation module. Uses either a port of 3photon's ranf random
number generator or the "rand" crate that is the Rust standard for RNGs.
"""
module Random

using StaticArrays: MVector

export RandomGenerator


# Generated random numbers will have a granularity of 1/MODULO
const Integer = Int32
const MODULO = 1000 * 1000 * 1000
const INV_MODULO = 1 / MODULO


"Random number generator, from Knuth's ranf (in Seminumerical Algorithm)"
struct RanfGenerator
    "Seed which this generator was initialized with"
    seed::Integer

    "Current set of random numbers, maps to IA in original code"
    numbers::MVector{56, Integer}

    "Index of the current random number, maps to MCALL in original code"
    index::UInt
end


# === RANDOM NUMBER GENERATION ===

"""
Generate 55 new random numbers between 0 and 1/FMODUL
This roughly maps to the IRN55 method in the original code.
"""
function reset!(rng::RanfGenerator)
    for i = 2:25
        rng.numbers[i] -= rng.numbers[i + 31]
        if rng.numbers[i] < 0
            rng.numbers[i] += MODULO
        end
    end
    for i = 26:56
        rng.numbers[i] -= rng.numbers[i - 24]
        if rng.numbers[i] < 0
            rng.numbers[i] += MODULO
        end
    end
end


# === CONSTRUCTION ===

"""
Create a new generator with an arbitrary seed.
This roughly maps to the IN55 method in the original code.
"""
function seeded_new(seed::Integer)::RanfGenerator
    # Start by zero-initializing the generator state
    #
    # FIXME: Isn't there any way to say which field we are talking about?
    #
    result = RanfGenerator(
        seed,
        zeros(56),  # numbers
        56,         # index
    )

    # Run the IN55 initialization algorithm
    result.numbers[56] = seed
    j = seed
    k = 1
    for i = 1:54
        ii = (21 * i) % 55 + 1
        result.numbers[ii] = k
        k = j - k
        if k < 0
            k += MODULO
        end
        j = result.numbers[ii]
    end

    # Warm up the sequence a bit
    for _ = 1:10
        reset!(result)
    end

    # Return the initialized generator
    result
end


"Create a new generator, with state faithful to original 3photons code"
function RanfGenerator()
    # TODO: Would be nice to figure out the seed constraints of seeded_new and
    #       publicize that interface too.
    seeded_new(Integer(234612947))
end


# === CHOICE OF RNG ===

# FIXME: Extract ranf into its a submodule and allow switching between ranf and
#        xoshiro algorithms like in the Rust version
"Random Number Generator implementation in use"
const RandomGenerator = RanfGenerator

end