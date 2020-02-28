# No include dependency
#
# FIXME: Isn't there a way to spell this out in code???


"Basic numerical concepts used throughout the program"
module Numeric

export ComplexF, Float


# FIXME: Make choice of floating-point type configurable
"Floating-point type used throughout the simulation."
const Float = Float64

"Complex number type with Float coordinates"
const ComplexF = Complex{Float}

end