## Containment
import Base.∈
∈(v::Vector{<:Real},p::Polyhedron) = contains(p,v)

function contains(p::Polyhedron,v::Vector{<:Real};tol = 1e-10)
    return contains(p.A,p.b,v;tol)
end
function contains(A,b,x;tol=1e-10)
    n,m = size(A)
    for i = 1:m
        (A[:,i]'*x > b[i]+tol) && return false
    end
    return true 
end

## isempty
import Base.isempty
function isempty(p::Polyhedron;sense=nothing)
    return isempty(p.A,p.b;sense)
end
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
## isfulldim
function isfulldim(p::Polyhedron;tol=1e-12)
    _,r = center(p) 
    return r > tol
end
## Projection
function project(x::Vector{<:Real}, p::Polyhedron)
    return project(x,p.A,p.b)
end
function project(x::Vector{<:Real},A::Matrix{<:Real}, b::Vector{<:Real})
    m = length(b)
    xp,_,_,_ = DAQP.quadprog(DAQP.QPj(zeros(0,0),x,A,b,fill(-1e30,m),zeros(Int32,m);A_rowmaj=true))
    return xp
end
## Size 
import Base: size
size(p::Polyhedron) = size(p.A)
