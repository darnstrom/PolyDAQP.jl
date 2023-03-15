using PolyDAQP
using LinearAlgebra
using Test

δ = 1e-8
@testset "Chebyshev center" begin
    n,m = 5,20
    cref,rref  = randn(n), 1e-2*rand()
    for p in [1,2,Inf]
        A = randn(n,m)
        Anorms = [norm(c,p) for c in eachcol(A)] 
        b = A'*cref+Anorms*rref
        c,r = center(A,b;p) 
        @test abs(r-rref) < δ 
    end
end
