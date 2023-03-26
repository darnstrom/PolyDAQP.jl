## Projection
function project(x::Vector{<:Real}, p::Polyhedron)
    return project(x,p.A,p.b)
end
function project(x::Vector{<:Real},A::Matrix{<:Real}, b::Vector{<:Real})
    m = length(b)
    xp,_,_,_ = DAQP.quadprog(DAQP.QPj(zeros(0,0),x,A,b,fill(-1e30,m),zeros(Int32,m);A_rowmaj=true))
    return xp
end
"""
    proj(p,x)
The projection of a vector `x` onto the Polyhedron `p`
"""
proj(p::Polyhedron,x::Vector{<:Real}) = project(x,p)
"""
    proj⊥(p,x)
The projection of a vector `x` onto the orthogonal complement of `p`
"""
proj⊥(p::Polyhedron,x::Vector{<:Real}) = x-project(x,p)
## Intersection
"""
    p ∩ q
Compute the intersection of polyhedra `p` and `q`
"""
function Base.intersect(p::Polyhedron,q::Polyhedron)
    return Polyhedron(hcat(p.A,q.A), vcat(p.b,q.b))
end
## Size 
Base.:size(p::Polyhedron) = size(p.A)
