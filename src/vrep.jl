## Naive vertex enumeration for 2D
function vrep_2d(A,b)
    n,m = size(A)
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
            if(iscontained(v,A,b))
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
    if(iscontained(v,A,b))
        push!(vs,v)
    else
        @warn "Could not close cycle"
    end
    all(visited) || @warn "Missing vertices"
    return vs
end
