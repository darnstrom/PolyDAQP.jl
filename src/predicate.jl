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
function Base.:isempty(p::Polyhedron;sense=nothing,validate=false, daqp_settings=nothing)
    return isempty(p.A,p.b;sense,validate)
end
function Base.:isempty(A::Matrix{<:Real}, b::Vector{<:Real};sense = nothing,
                       validate=false, daqp_settings=nothing)
    nth,m = size(A)  
    bl = fill(-1e30,m)
    sense = isnothing(sense) ? zeros(Int32,m) : sense
    x,fval,exitflag,info = DAQP.quadprog(DAQP.QPj(zeros(0,0),zeros(0),A,b,bl,sense;A_rowmaj=true),
                                        settings=daqp_settings);
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

## isfulldim
"""
    isfulldim(p)
Return true if `p` is full dimensional
"""
function isfulldim(p::Polyhedron;tol=1e-12)
    _,r = center(p) 
    return r > tol
end
