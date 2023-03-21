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

## Linear transformation
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
