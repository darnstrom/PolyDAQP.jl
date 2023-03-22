## Naive vertex enumeration for 2D
function vrep_2d(p::Polyhedron;tol=1e-10)
    return vrep_2d(p.A,p.b;tol)
end

function vrep_2d(A::Matrix{<:Real},b::Vector{<:Real};tol = 1e-10)
    n,m = size(A)
    n != 2 && error("vrep currently supported for two-dimensional polyhedra") 
    vs = Vector{Float64}[]
    m==0 && return vs
    for i = 1:m
        for j = i+1:m 
            L = lu!([A[:,i]'; A[:,j]'],check=false)
            r = [b[i];b[j]]
            !issuccess(L) && continue
            v = ldiv!(L,r)
            if(contains(A,b,v;tol))
                push!(vs,v)
            end
        end
    end
    return convexhull(vs)
end

function convexhull(vs)
    nv = length(vs)
    nv == 0 && return vs
    clockwise(p, q, r) = 0 < ((r[2] - p[2]) * (q[1] - p[1]) - (r[1] - p[1]) * (q[2] - p[2]))
    # Gift wrapping
    hull = [argmin(first.(vs))] # Start with left-most point 

    while length(hull) == 1 || hull[1] != hull[end] 
        next = hull[end] % nv + 1 # select the preceeding point as next candidate 
        for j = 1:nv # compare with other points 
            if clockwise(vs[hull[end]], vs[next], vs[j]) # If "more to the left", switch
                next = j
            end
        end
        push!(hull, next)
    end
    return vs[hull]
end

