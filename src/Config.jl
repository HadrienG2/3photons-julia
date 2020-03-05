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
function Base.show(io::IO, cfg::Configuration)
    println(io, "ITOT           : ", cfg.num_events)
    println(io, "ETOT           : ", cfg.e_tot)
    println(io, "oCutpar.ACUT   : ", cfg.event_cut.a_cut)
    println(io, "oCutpar.BCUT   : ", cfg.event_cut.b_cut)
    println(io, "oCutpar.EMIN   : ", cfg.event_cut.e_min)
    @printf(io, "oCutpar.SINCUT : %.0f\n", cfg.event_cut.sin_cut)
    println(io, "ALPHA          : ", cfg.𝛼)
    println(io, "ALPHAZ         : ", cfg.𝛼_Z)
    @printf(io, "CONVERS        : %.0f\n", cfg.convers)
    println(io, "oParam.MZ0     : ", cfg.m_Z⁰)
    println(io, "oParam.GZ0     : ", cfg.g_Z⁰)
    println(io, "SIN2W          : ", cfg.sin²_w)
    println(io, "BREPEM         : ", cfg.br_e₊_e₋)
    @printf(io, "BETAPLUS       : %.0f\n", cfg.𝛽₊)
    @printf(io, "BETAMOINS      : %.0f\n", cfg.𝛽₋)
    println(io, "NBIN           : ", cfg.n_bin)
    println(io, "oParam.IMPR    : ", cfg.impr)
    println(io, "PLOT           : ", cfg.plot)
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

        # Fetch and decode the next config file entry, given its name & type
        next = iterate(first_field_it)
        function next_entry!(name::AbstractString, ty::DataType)
            field = try
                (field, state) = next
                next = iterate(first_field_it, state)
                field
            catch e
                println("Missing configuration of $name")
                rethrow(e)
            end
            try
                parse_entry(ty, field)
            catch e
                println("Could not parse configuration of $name")
                rethrow(e)
            end
        end
        
        # Decode the configuration
        #
        # FIXME: Isn't there any way to say which field we are talking about?
        #
        config = Configuration(
            next_entry!("num_events", UInt),
            next_entry!("e_tot", Float),
            EventCut(
                next_entry!("a_cut", Float),
                next_entry!("b_cut", Float),
                next_entry!("e_min", Float),
                next_entry!("sin_cut", Float),
            ),
            next_entry!("𝛼", Float),
            next_entry!("𝛼_Z", Float),
            next_entry!("convers", Float),
            next_entry!("m_Z⁰", Float),
            next_entry!("g_Z⁰", Float),
            next_entry!("sin²_w", Float),
            next_entry!("br_e₊_e₋", Float),
            next_entry!("𝛽₊", Float),
            next_entry!("𝛽₋", Float),
            next_entry!("n_bin", Int32),
            next_entry!("impr", Bool),
            next_entry!("plot", Bool),
        )

        # Display it like the C++ version would (this eases comparisons)
        if !jit_warmup
            print(config)
        end

        # A sensible simulation must run for at least one event
        @enforce (config.num_events > 0) "Please simulate at least one event"

        # We don't support the original code's PAW-based plotting features, so
        # we make sure that they weren't enabled.
        #
        # TODO: Consider bringing back plotting in the Julia version later, that
        #       should be easier than in the Rust version.
        #
        @enforce (!config.plot) "Plotting is not supported by this version"

        # We don't support the initial code's debugging feature which displays
        # all intermediary results during sampling. Such a feature should be set
        # up at build time to avoid run-time costs.
        #
        # TODO: Revise this Rust version design decision for the Julia version.
        #       With JITting, compile/run time tradeoffs work differently.
        #
        @enforce (!config.impr) """
        Individual result printing is not supported.
        This debugging feature has a run-time performance cost even when unused,
        so it should be implemented at compile-time instead.
        """

        # We're done checking the conifguration and can return it
        config
    end
end

end