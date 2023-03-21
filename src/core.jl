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
import Base.isempty
function isempty(A::Matrix{<:Real}, b::Vector{<:Real};sense = nothing)
    nth,m = size(A)  
    bl = fill(-1e30,m)
    sense = isnothing(sense) ? zeros(Int32,m) : sense
    x,fval,exitflag,info = DAQP.quadprog(DAQP.QPj(zeros(0,0),zeros(0),A,b,bl,sense;A_rowmaj=true));
    return exitflag == -1 
end

function isfeasible(p::Ptr{Cvoid}; m=nothing, ms=nothing)::Bool
    # Update A and bupper/blower dimensions 
    !isnothing(m)  && unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,3)),m)
    !isnothing(ms) && unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,4)),ms)

    exitflag =ccall((:daqp_ldp,DAQP.libdaqp), Int32, (Ptr{Cvoid},),p);

    # Reset the workspace 
    ccall((:deactivate_constraints,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),p);
    ccall((:reset_daqp_workspace,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),p);
    return exitflag==1;
end
## Projection
function project(x::Vector{<:Real},A::Matrix{<:Real}, b::Vector{<:Real})
    m = length(b)
    xp,_,_,_ = DAQP.quadprog(DAQP.QPj(zeros(0,0),x,A,b,fill(-1e30,m),zeros(Int32,m);A_rowmaj=true))
    return xp
end
## Minkowski
import Base.sum
function sum(p::Polyhedron, q::Polyhedron)
    np,mp = size(p.A)
    nq,mq = size(q.A)
    Alift = [zeros(nq,mp) q.A;p.A -q.A] 
    blift = [p.b;q.b]
    As,bs = eliminate(Alift,blift,collect(np+1:2np))
    return Polyhedron(As,bs)
end

function sum(p::Polyhedron, v::Vector{<:Real})
    return Polyhedron(p.A,p.b+p.A'*v)
end

import Base: + 
+(p::Polyhedron,q::Polyhedron) = sum(p,q)
⊕(p::Polyhedron,q::Polyhedron) = sum(p,q)
+(p::Polyhedron,v::Vector{<:Real}) = sum(p,v)
+(v::Vector{<:Real},p::Polyhedron) = sum(p,v)

## Linear  transformation
function linear_transform(p::Polyhedron, F::Matrix{<:Real})
    m,n = size(F)
    if(m <= n)
        Q = qr(F')
        R1 = Q.R 
        ny = size(R1,1)
        nz = n-size(R1,1)
        Q1,Q2 = Q.Q[:,1:ny],Q.Q[:,ny+1:end];
        Alift = [inv(R1)*Q1'*p.A; Q2'*p.A] 

        Ar,br = eliminate(Alift,p.b,collect(ny+1:ny+nz))
        return Polyhedron(Ar,br)
    else # m > n
        #TODO (will lead to a lower dimensional...) 
    end
end
import Base: * 
*(F::Matrix{<:Real},p::Polyhedron) = linear_transform(p,F)
*(p::Polyhedron, F::Matrix{<:Real}) = Polyhedron(F'*p.A,p.b) 
## Size 
import Base: size
size(p::Polyhedron) = size(p.A)
