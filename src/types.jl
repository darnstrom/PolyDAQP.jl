"""
Polyedron defined by ``A x ≤ b``
"""
mutable struct Polyhedron 
    A::Matrix{Float64}
    b::Vector{Float64}
end
