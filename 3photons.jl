# Run with julia --project=. 3photons.jl
# FIXME: There has to be a better way
import IterTools

# TODO: Break this into multiple modules
# TODO: After translating, turns this into more idiomatic Julia (e.g. unicode
#       variable names, more genericity...

"Floating-point type used throughout the simulation."
const Float = Float64

"Cuts on generated events"
struct EventCut
    "Cut on maximum cosine of (beam, photons) angle"
    a_cut::Float
    
    "Cut on maximum cosine of (photon, photon) angle"
    b_cut::Float
    
    "Cut on minimum photon energy"
    e_min::Float
    
    "Cut on minimum cosine of (beam, normal to the photon plane) angle"
    sin_cut::Float
end

"Simulation configuration"
struct Configuration
    "Number of events to be simulated"
    num_events::UInt
    
    "Collision energy at center of mass (GeV)"
    e_tot::Float
    
    "Cuts on the angles and energies of generated photons"
    event_cut::EventCut
    
    "Fine structure constant"
    alpha::Float
    
    "Fine structure constant at the Z peak"
    alpha_z::Float
    
    "Conversion factor from GeV^(-2) to pb"
    convers::Float
    
    "Z⁰ boson mass (GeV)"
    m_z0::Float
    
    "Z⁰ boson width (GeV)"
    g_z0::Float
    
    "Square sine of Weinberg's Theta"
    sin2_w::Float
    
    "Branching factor from Z to e+/e-"
    br_ep_em::Float
    
    "Beta + (???)"
    beta_plus::Float
    
    "Beta - (???)"
    beta_minus::Float
    
    "Number of histogram bins (UNUSED)"
    n_bin::Int32
    
    "Whether intermediary results should be displayed (UNUSED)"
    impr::Bool
    
    "Whether results should be plotted in a histogram (UNUSED)"
    plot::Bool
end

"Constructor that loads the simulation configuration from a file"
function Configuration(file_name::AbstractString)
    open(file_name) do config_file
        # Iterate over the config file's lines, extracting the first chunk
        # of non-whitespace on each line (if any) and skipping empty lines.
        line_it = eachline(config_file)
        split_first_it = IterTools.imap(s -> split(s; limit=2), line_it)
        non_empty_it = Iterators.filter(a -> length(a) > 0, split_first_it)
        first_field_it = IterTools.imap(a -> a[1], non_empty_it)
        
        # DEBUG printout
        println("DEBUG: Here are the raw configuration items")
        for line in first_field_it
            println(line)
        end
    end
end

# DEBUG
config = Configuration("valeurs")
