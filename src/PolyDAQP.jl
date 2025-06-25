module PolyDAQP
using DAQPBase, LinearAlgebra
const DAQP = DAQPBase

# Types
include("types.jl")
export Polyhedron

# Methods
include("boundingbox.jl")
export bounding_box

include("center.jl")
export center

include("core.jl")
export project,proj,proj⊥,intersect


include("eliminate.jl")
export eliminate

include("mappings.jl")
export sum,+,⊕,*

include("minrep.jl")
export minrep

include("plot.jl")
export pplot

include("predicate.jl")
export contains,isempty,isfulldim,∈

include("slice.jl")
export slice

include("vrep.jl")

include("normalize.jl")
export normalize!

end # module
