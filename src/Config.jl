# Depends on Errors.jl, EvCut.jl and Numeric.jl being include-d beforehand
#
# FIXME: Isn't there a way to spell this out in code???


"Mechanism for loading and sharing the simulation configuration"
module Config

import IterTools

using ..Errors: @enforce
using ..EvCut: EventCut
using ..Numeric: Float
using Printf: @printf

export Configuration


"Simulation configuration"
struct Configuration
    "Number of events to be simulated"
    num_events::UInt
    
    "Collision energy at center of mass (GeV)"
    e_tot::Float
    
    "Cuts on the angles and energies of generated photons"
    event_cut::EventCut
    
    "Fine structure constant"
    𝛼::Float
    
    "Fine structure constant at the Z peak"
    𝛼_Z::Float
    
    "Conversion factor from GeV^(-2) to pb"
    convers::Float
    
    "Z⁰ boson mass (GeV)"
    m_Z⁰::Float
    
    "Z⁰ boson width (GeV)"
    g_Z⁰::Float
    
    "Square sine of Weinberg's Theta"
    sin²_w::Float
    
    "Branching factor from Z to e+/e-"
    br_e₊_e₋::Float
    
    "Beta + (???)"
    𝛽₊::Float
    
    "Beta - (???)"
    𝛽₋::Float
    
    "Number of histogram bins (UNUSED)"
    n_bin::Int32
    
    "Whether intermediary results should be displayed (UNUSED)"
    impr::Bool
    
    "Whether results should be plotted in a histogram (UNUSED)"
    plot::Bool
end


"Display the configuration, following formatting of the original 3photons code"
function Base.show(io::IO, c::Configuration)
    println(io, "ITOT           : ", c.num_events)
    println(io, "ETOT           : ", c.e_tot)
    println(io, "oCutpar.ACUT   : ", c.event_cut.a_cut)
    println(io, "oCutpar.BCUT   : ", c.event_cut.b_cut)
    println(io, "oCutpar.EMIN   : ", c.event_cut.e_min)
    @printf(io, "oCutpar.SINCUT : %.0f\n", c.event_cut.sin_cut)
    println(io, "ALPHA          : ", c.𝛼)
    println(io, "ALPHAZ         : ", c.𝛼_Z)
    @printf(io, "CONVERS        : %.0f\n", c.convers)
    println(io, "oParam.MZ0     : ", c.m_Z⁰)
    println(io, "oParam.GZ0     : ", c.g_Z⁰)
    println(io, "SIN2W          : ", c.sin²_w)
    println(io, "BREPEM         : ", c.br_e₊_e₋)
    @printf(io, "BETAPLUS       : %.0f\n", c.𝛽₊)
    @printf(io, "BETAMOINS      : %.0f\n", c.𝛽₋)
    println(io, "NBIN           : ", c.n_bin)
    println(io, "oParam.IMPR    : ", c.impr)
    println(io, "PLOT           : ", c.plot)
end


"Constructor that loads the simulation configuration from a file"
function Configuration(file_name::AbstractString; jit_warmup::Bool=false)
    # Open the config file
    open(file_name) do config_file
        # Parse a config file entry, given its type
        function parse_entry(ty::DataType, field::AbstractString)
            # 3photons's config file uses an unusual syntax for booleans
            if ty == Bool
                # Special-case this boolean syntax
                if field == ".true."
                    true
                elseif field ==".false."
                    false
                else
                    throw(DomainError(field, "Not a 3photons-style boolean"))
                end
            else
                # Use julia's standard string parser for anything else
                parse(ty, field)
            end
        end

        # Iterate over the config file's lines, extracting the first chunk
        # of non-whitespace on each line (if any) and skipping empty lines.
        line_it = eachline(config_file)
        split_first_it = IterTools.imap(s -> split(s; limit=2), line_it)
        non_empty_it = Iterators.filter(a -> length(a) > 0, split_first_it)
        first_field_it = IterTools.imap(a -> a[1], non_empty_it)

        # Fetch and decode the next config file entry, given its type
        next = iterate(first_field_it)
        function next_entry!(ty::DataType)
            (field, state) = next
            next = iterate(first_field_it, state)
            parse_entry(ty, field)
        end
        
        # Decode the configuration
        #
        # FIXME: Isn't there any way to say which field we are talking about?
        #
        config = Configuration(
            next_entry!(UInt),       # num_events
            next_entry!(Float),      # e_tot
            EventCut(
                next_entry!(Float),  # a_cut
                next_entry!(Float),  # b_cut
                next_entry!(Float),  # e_min
                next_entry!(Float),  # sin_cut
            ),
            next_entry!(Float),      # 𝛼
            next_entry!(Float),      # 𝛼_Z
            next_entry!(Float),      # convers
            next_entry!(Float),      # m_Z⁰
            next_entry!(Float),      # g_Z⁰
            next_entry!(Float),      # sin²_w
            next_entry!(Float),      # br_e₊_e₋
            next_entry!(Float),      # 𝛽₊
            next_entry!(Float),      # 𝛽₋
            next_entry!(Int32),      # n_bin
            next_entry!(Bool),       # impr
            next_entry!(Bool),       # plot
        )

        # Display it like the C++ version would (this eases comparisons)
        if !jit_warmup
            print(config)
        end

        # A sensible simulation must run for at least one event
        @enforce (config.num_events > 0) "Invalid event count"

        # We don't support the original code's PAW-based plotting features, so
        # we make sure that they weren't enabled.
        #
        # FIXME: Consider enabling some plotting in the Julia version later.
        #
        @enforce (!config.plot) "Plotting is not supported"

        # We don't support the initial code's debugging feature which displays
        # all intermediary results during sampling. Such a feature should be set
        # up at build time to avoid run-time costs.
        #
        # FIXME: Revise this design decision for the Julia version. There is no
        #        compile time versus run-time distinction here.
        #
        @enforce (!config.impr) """
        Individual result printing is not supported.
        This debugging feature has a run-time performance cost even when unused.
        It should be implemented at compile-time instead.
        """

        # We're done checking the conifguration and can return it
        config
    end
end

end