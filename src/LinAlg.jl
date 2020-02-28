# No include dependency
#
# FIXME: Isn't there a way to spell this out in code???


"Common linear algebra declarations, e.g. 4-momentum handling logic"
module LinAlg

using StaticArrays: SVector

export E, X, XYZ, Y, Z


"X impulsion of a 4-momentum"
const X = 1

"Y impulsion of a 4-momentum"
const Y = 2

"Z impulsion of a 4-momentum"
const Z = 3

"Energy of a 4-momentum"
const E = 4

# FIXME: It's sad that efficient StaticArray slicing requires this weirdness
"Spatial part of a 4-momentum"
const XYZ = SVector{3}(X:Z)

end