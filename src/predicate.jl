## Containment
"""
    x ∈ p
Return true if `x` is contained in the polyhedron `p`
"""
Base.:∈(v::Vector{<:Real},p::Polyhedron) = contains(p,v)

function Base.:contains(p::Polyhedron,v::Vector{<:Real};tol = 1e-10)
    return contains(p.A,p.b,v;tol)
end
function Base.:contains(A::Matrix{<:Real},b::Vector{<:Real},x::Vector{<:Real} ;tol=1e-10)
    n,m = size(A)
    for i = 1:m
        (A[:,i]'*x > b[i]+tol) && return false
    end
    return true 
end

## isempty
"""
    isempty(p)
Return true if `p` ``\\neq \\emptyset`` 
"""
function Base.:isempty(p::Polyhedron;sense=nothing,validate=false)
    return isempty(p.A,p.b;sense,validate)
end
function Base.:isempty(A::Matrix{<:Real}, b::Vector{<:Real};sense = nothing,validate=false)
    nth,m = size(A)  
    bl = fill(-1e30,m)
    sense = isnothing(sense) ? zeros(Int32,m) : sense
    x,fval,exitflag,info = DAQP.quadprog(DAQP.QPj(zeros(0,0),zeros(0),A,b,bl,sense;A_rowmaj=true));
    if(validate && exitflag == -1)
        @inbounds begin
            Atlam = zeros(nth)
            btlam = 0
            for i = 1:m
                if(abs(info.λ[i]) > 0)
                    Atlam += info.λ[i]*view(A,:,i)
                    btlam += b[i]*info.λ[i]
                end
            end
            if(norm(Atlam)+btlam >0)
                @warn "Could not validate infeasibility through Farkas lemma"
            end
        end
    end
    return exitflag == -1 
end

function isfeasible(p::Ptr{Cvoid}; m=nothing, ms=nothing,validate=false)::Bool
    # Update A and bupper/blower dimensions 
    !isnothing(m)  && unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,3)),m)
    !isnothing(ms) && unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,4)),ms)

    exitflag =ccall((:daqp_ldp,DAQP.libdaqp), Int32, (Ptr{Cvoid},),p);

    if(validate && exitflag == -1)
        daqp_ws = unsafe_load(Ptr{DAQP.Workspace}(p))
        m,ms,n = daqp_ws.m, daqp_ws.ms, daqp_ws.n
        A = unsafe_wrap(Matrix{Cdouble}, daqp_ws.M, (n,m), own=false)
        b= unsafe_wrap(Vector{Cdouble}, daqp_ws.dupper, m+ms, own=false)
        AS = copy(unsafe_wrap(Vector{Cint}, daqp_ws.WS, daqp_ws.n_active, own=false))
        AS .+= 1 # Offset for C -> Julia
        lam_star = unsafe_wrap(Vector{Cdouble}, daqp_ws.lam_star, daqp_ws.n_active, own=false)
        err  = dot(b[AS],lam_star)+norm(A[:,AS]*lam_star)
        if(err>0)
            @warn "Couldn't validate infeas. with Frakas (err:$(err), fval=$(daqp_ws.fval))"
        end
    end

    # Reset the workspace 
    ccall((:deactivate_constraints,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),p);
    ccall((:reset_daqp_workspace,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),p);
    #if(exitflag != 1 && exitflag != -1)
    #    @warn "exitflag: $exitflag"
    #end
    return exitflag == 1
end
## isfulldim
"""
    isfulldim(p)
Return true if `p` is full dimensional
"""
function isfulldim(p::Polyhedron;tol=1e-12)
    _,r = center(p) 
    return r > tol
end
