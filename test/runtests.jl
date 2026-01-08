using PolyDAQP
using LinearAlgebra
using Test
using RecipesBase

δ = 1e-8

@testset "Center" begin
    n,m = 10,100
    cref,rref  = randn(n), 1e-2*rand()
    for p in [1,2,Inf]
        A = randn(n,m)
        Anorms = [norm(c,p) for c in eachcol(A)] 
        b = A'*cref+Anorms*rref
        q = Polyhedron(A,b)
        c,r = center(q;p) 
        @test abs(r-rref) < δ 
    end
end

@testset "Bounding box " begin
    n,m = 10,100
    A = randn(n,m)
    b = rand(m)
    p = Polyhedron(A,b)
    ub,lb = bounding_box(p) 
    @test all(ub.>=lb)


    # Simple 2D example to get vertices 
    A,b = randn(2,100), rand(100)
    ub,lb = bounding_box(A,b)
    Am,bm = minrep(A,b)
    pm = Polyhedron(Am,bm)
    vs = PolyDAQP.vrep_2d(pm) 
    ubref = maximum(reduce(hcat,vs),dims=2)
    lbref = minimum(reduce(hcat,vs),dims=2)
    @test all(abs.(ub-ubref) .< δ)
    @test all(abs.(lb-lbref) .< δ)
end

@testset "Minrep" begin
    n,m = 10,100
    A,b = randn(n,m), rand(m)
    p = Polyhedron(A,b)
    pm = minrep(p)
    @test length(pm.b) < m
    # Make sure the bounding box is the same 
    ub,lb= bounding_box(A,b) 
    ubm,lbm = bounding_box(pm.A,pm.b) 
    @test all(abs.(ub-ubm) .< δ)
    @test all(abs.(lb-lbm) .< δ)
end

@testset "Elimination" begin
    n,m = 6,40
    A,b = randn(n,m), rand(m)
    p = Polyhedron(A,b)
    ub,lb= bounding_box(p) 

    pn = eliminate(p,[4,5,6])
    ubn,lbn = bounding_box(pn) 
    @test all(abs.(ub[1:3]-ubn) .< δ)
    @test all(abs.(lb[1:3]-lbn) .< δ)
end

@testset "Affine maps" begin
    n,m = 10,100
    A,b = randn(n,m), rand(m)
    ub,lb= bounding_box(A,b) 

    v = randn(n)
    pv = Polyhedron(A,b)+v
    ubv,lbv = bounding_box(pv.A,pv.b)
    @test all(abs.(ub+v-ubv) .< δ)
    @test all(abs.(lb+v-lbv) .< δ)


    # Simple 2D example to get vertices 
    n,m = 2,100;
    A,b = randn(n,m), rand(m)
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


    n,m = 3,25
    p = Polyhedron(randn(n,m),rand(m))
    q = Polyhedron(randn(n,m),rand(m))
    ppq = minrep(p+q)
    ppq2 = minrep(p⊕q)
    @test size(ppq) == size(ppq2)

    c,r = center(p)
    F = randn(n,n)
    
    @test F\c ∈ (p*F)

end

@testset "Predicates" begin
    n,m = 10,100
    c = randn(n) 
    A = randn(n,m)
    b = A'*c+rand(m)
    p = Polyhedron(A,b)
    @test c ∈ p 
    @test !isempty(p)
    @test isfulldim(p)

    x = 100*randn(5)
end

@testset "Projection" begin
    n,m = 10,100
    x = 100*randn(n) 
    A,b = randn(n,m), rand(m)
    p = Polyhedron(A,b)
    @test proj(p,x) ∈ p 
    @test proj⊥(p,x) ∉ p 
    @test norm(x-(proj(p,x)+proj⊥(p,x))) < δ
end

@testset "Intersection" begin
    n,m = 5,5
    x = randn(n)
    p = Polyhedron(randn(n,m), rand(m))
    q = Polyhedron(randn(n,m), rand(m))
    pq = p ∩ q
    xpq = proj(pq,x)

    @test xpq ∈ p 
    @test xpq ∈ q 
end

@testset "Plotting" begin
    n,m = 10,100
    A,b = randn(n,m), rand(m)
    p = slice(Polyhedron(A,b),collect(3:n))
    A,b = randn(n,m), rand(m)
    q = slice(Polyhedron(A,b),collect(3:n))
    pl = pplot([p, q]);

    d = RecipesBase.apply_recipe(Dict{Symbol,Any}(),p)
    d = RecipesBase.apply_recipe(Dict{Symbol,Any}(),[(p,(v->ones(2)'v+1))])
    d = RecipesBase.apply_recipe(Dict{Symbol,Any}(),[p,q])
end
