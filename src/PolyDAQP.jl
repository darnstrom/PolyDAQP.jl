module PolyDAQP
using DAQP, LinearAlgebra

include("core.jl")
export slice,iscontained,isfeasible,sum,+,⊕

include("center.jl")
export center

include("eliminate.jl")
export center

include("minrep.jl")
export minrep
include("vrep.jl")

include("boundingbox.jl")
export bounding_box 

include("types.jl")
export Polyhedron

include("plot.jl")
export pplot

end # module
