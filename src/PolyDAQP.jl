module PolyDAQP
using DAQP, LinearAlgebra

export bounding_box

include("types.jl")
export Polyhedron

include("core.jl")
export contains,isempty,isfulldim,isfeasible,∈,project,proj,proj⊥,intersect

include("center.jl")
export center

include("eliminate.jl")
export eliminate

include("minrep.jl")
export minrep
include("vrep.jl")

include("boundingbox.jl")
#export bounding_box 

include("plot.jl")
export pplot

include("slice.jl")
export slice

include("mappings.jl")
export sum,+,⊕,*

end # module
