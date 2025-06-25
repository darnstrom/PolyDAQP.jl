function LinearAlgebra.normalize!(p::Polyhedron;start=1)
    normalize!(p.A,p.b;start) 
end
function LinearAlgebra.normalize!(A::Matrix{<:Real},b::Vector{<:Real};start=1)
    for i = start:length(b)
        norm_factor = norm(view(A,:,i),2);
        rdiv!(view(A,:,i),norm_factor)
        b[i]/=norm_factor
    end
end

