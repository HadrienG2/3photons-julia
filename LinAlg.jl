# No include dependency
#
# FIXME: Isn't there a way to spell this out in code???


"Common linear algebra declarations, e.g. 4-momentum handling logic"
module LinAlg

using StaticArrays: SVector

export X, XYZ, Y, Z, E


"X impulsion of a 4-momentum"
const X = 1

"Y impulsion of a 4-momentum"
const Y = 2

"Z impulsion of a 4-momentum"
const Z = 3

"Energy of a 4-momentum"
const E = 4

"Helper for efficiently slicing the spatial part of a 4-momentum"
const XYZ = SVector{3}(X:Z)

end