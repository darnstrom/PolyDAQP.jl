## Bounding box
"""
    ub,lb = bounding_box(p)
Return `ub`,`lb` that define the bounding box of the Polyhedron `p`
(``p \\subseteq \\{x: lb ≤ x ≤ ub\\}``)
"""
function bounding_box(p::Polyhedron ;sense=nothing)
    return bounding_box(p.A,p.b;sense)
end

function bounding_box(A::Matrix{<:Real}, b::Vector{<:Real};sense=nothing)
    n,m = size(A)
    bl = fill(-1e30,m)
    sense = isnothing(sense) ? zeros(Int32,m) : sense
    ub,lb = zeros(n),zeros(n)
    for i = 1:n 
        f = zeros(n)
        for (fi,bound) in [(1,lb),(-1,ub)]
            f[i] = fi 
            x,fval,exitflag,info = DAQP.linprog(f,A,b,bl,sense;A_rowmaj=true)
            bound[i] = (exitflag==1) ? x[i] : -fi*Inf 
        end
    end
    return ub,lb 
end
