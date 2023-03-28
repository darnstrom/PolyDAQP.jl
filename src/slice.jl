## Slice
function slice(p::Polyhedron,ids::Vector{<:Integer};values=nothing, tol_zero=1e-10)
    return Polyhedron(slice(p.A,p.b,ids;values,tol_zero)...);
end
function slice(A::Matrix{<:Real}, b::Vector{<:Real},ids::Vector{<:Integer};values = nothing,tol_zero=1e-10)
    n,m = size(A)
    ids_keep = setdiff(collect(1:n),ids)
    values = isnothing(values) ? zeros(length(ids)) : values
    rows = Int[]
    scalings = Float64[]
    bout = Float64[]
    for i = 1:m 
        anorm = norm(view(A,ids_keep,i))
        bi = b[i]-A[ids,i]'*values;
        if(anorm > tol_zero)
            push!(rows,i)
            push!(scalings,1/anorm)
            push!(bout,bi/anorm)
        elseif(bi < -tol_zero) # infeasible or 0 == 0 
            return zeros(n-length(ids),0),zeros(0) 
        end

    end
    return A[ids_keep,rows]*Diagonal(scalings), bout 
end
