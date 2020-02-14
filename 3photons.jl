# NOTE: Run with julia --project=. 3photons.jl to load proper packages
#
# FIXME: There has to be a better way. If there isn't, replace usage of
#        IterTools with eager collection usage, it's not that important.

# TODO: After translating, turns this into more idiomatic Julia (e.g. unicode
#       variable names, more genericity...


include("config.jl")

using .Config: Configuration


# Load the configuration from its file
cfg = Configuration("valeurs")
