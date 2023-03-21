module PolyDAQP
using DAQP, LinearAlgebra

include("types.jl")
export Polyhedron

include("core.jl")
export slice,iscontained,isfeasible,sum,+,⊕

include("center.jl")
export center

include("eliminate.jl")
export eliminate

include("minrep.jl")
export minrep
include("vrep.jl")

include("boundingbox.jl")
export bounding_box 

include("plot.jl")
export pplot

end # module
