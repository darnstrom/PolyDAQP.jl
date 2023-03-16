## Slice
function slice(A,b,ids;values = nothing)
    n,m = size(A)
    ids_keep = setdiff(collect(1:n),ids)
    values = isnothing(values) ? zeros(length(ids)) : values
    rows = Int[]
    scalings = Float64[]
    bout = Float64[]
    for i = 1:m 
        anorm = norm(view(A,ids_keep,i))
        bi = b[i]-A[ids,i]'*values;
        if(anorm > 1e-10)
            push!(rows,i)
            push!(scalings,1/anorm)
            push!(bout,bi/anorm)
        elseif(bi < -1e-10) # infeasible or 0 == 0 
            return zeros(n-length(ids),0),zeros(0) 
        end

    end
    return A[ids_keep,rows]*Diagonal(scalings), bout 
end

## Containment
function iscontained(x,A,b;tol=1e-10)
    n,m = size(A)
    for i = 1:m
        (A[:,i]'*x > b[i]+tol) && return false
    end
    return true 
end

## isempty
function poly_isempty(A::Matrix{<:Real}, b::Vector{<:Real};sense = nothing)
    nth,m = size(A)  
    bl = fill(-1e30,m)
    sense = isnothing(sense) ? zeros(Int32,m) : sense
    x,fval,exitflag,info = DAQP.quadprog(DAQP.QPj(zeros(0,0),zeros(0),A,b,bl,sense;A_rowmaj=true));
    return exitflag == -1 
end

