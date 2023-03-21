## Naive vertex enumeration for 2D
function vrep_2d(p::Polyhedron)
    return vrep_2d(p.A,p.b)
end
function vrep_2d(A::Matrix{<:Real},b::Vector{<:Real})
    n,m = size(A)
    n != 2 && error("vrep currently supported for two-dimensional polyhedra") 
    vs = Vector{Float64}[]
    m==0 && return vs
    visited,i = falses(m),1
    visited[i] = true
    for out = 1:m
        for j = 1:m 
            visited[j] && continue
            L = lu!([A[:,i]'; A[:,j]'],check=false)
            r = [b[i];b[j]]
            !issuccess(L) && continue
            v = ldiv!(L,r)
            if(contains(A,b,v))
                push!(vs,v)
                i, visited[j] = j, true
                break
            end
        end
    end
    # The last point has to connect to the first facet 
    L = lu!([A[:,i]'; A[:,1]'],check=false)
    r = [b[i];b[1]]
    v = ldiv!(L,r)
    if(contains(A,b,v))
        push!(vs,v)
    else
        @warn "Could not close cycle"
    end
    all(visited) || @warn "Missing vertices"
    return vs
end
