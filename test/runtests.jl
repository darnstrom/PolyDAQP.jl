using PolyDAQP
using LinearAlgebra
using Test

δ = 1e-8

@testset "Center" begin
    n,m = 10,100
    cref,rref  = randn(n), 1e-2*rand()
    for p in [1,2,Inf]
        A = randn(n,m)
        Anorms = [norm(c,p) for c in eachcol(A)] 
        b = A'*cref+Anorms*rref
        c,r = center(A,b;p) 
        @test abs(r-rref) < δ 
    end
end

@testset "Bounding box " begin
    n,m = 10,100
    A = randn(n,m)
    b = rand(m)
    ub,lb = bounding_box(A,b) 
    @test all(ub.>=lb)


    # Simple 2D example to get vertices 
    A,b = randn(2,100), rand(100)
    ub,lb = bounding_box(A,b)
    Am,bm = minrep(A,b)
    vs = PolyDAQP.vrep_2d(Am,bm) 
    ubref = maximum(reduce(hcat,vs),dims=2)
    lbref = minimum(reduce(hcat,vs),dims=2)
    @test all(abs.(ub-ubref) .< δ)
    @test all(abs.(lb-lbref) .< δ)
end

@testset "Minrep" begin
    n,m = 10,100
    A = randn(n,m)
    b = rand(m)
    Am,bm = minrep(A,b)
    @test length(bm) < m
    # Make sure the bounding box is the same 
    ub,lb= bounding_box(A,b) 
    ubm,lbm = bounding_box(A,b) 
    @test all(abs.(ub-ubm) .< δ)
    @test all(abs.(lb-lbm) .< δ)
end

