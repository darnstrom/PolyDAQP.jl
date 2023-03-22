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

@testset "Elimination" begin
    n,m = 6,40
    A = randn(n,m)
    b = rand(m)
    ub,lb= bounding_box(A,b) 

    An,bn = eliminate(A,b,[4,5,6])
    ubn,lbn = bounding_box(An,bn) 
    @test all(abs.(ub[1:3]-ubn) .< δ)
    @test all(abs.(lb[1:3]-lbn) .< δ)
end

@testset "Affine maps" begin
    n,m = 10,100
    A = randn(n,m)
    b = rand(m)
    ub,lb= bounding_box(A,b) 

    v = randn(n)
    pv = Polyhedron(A,b)+v
    ubv,lbv = bounding_box(pv.A,pv.b)
    @test all(abs.(ub+v-ubv) .< δ)
    @test all(abs.(lb+v-lbv) .< δ)


    # Simple 2D example to get vertices 
    n,m = 2,50;
    A,b = randn(2,100), rand(100)
    Am,bm = minrep(A,b)
    vs = PolyDAQP.vrep_2d(Am,bm) 
    F = randn(2,2)
    Fvsr = [F*v for v in vs]
    p = Polyhedron(Am,bm) 
    Fp = F*p
    Fvs= PolyDAQP.vrep_2d(minrep(Fp.A,Fp.b)...)

    errors = Float64[] 
    for v in Fvs 
        push!(errors, minimum(norm(v-vr) for vr in Fvsr))
    end
    @test all(errors .< δ)

end
