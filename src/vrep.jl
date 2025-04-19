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
    if(isempty(vs)) return vs;
    c = sum(vs)/length(vs) 
    return sort(vs, by= x->atan(x[2]-c[2],x[1]-c[1]))
end
