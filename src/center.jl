## Chebyshev center
function center(q::Polyhedron;sense=nothing, isnormalized = false, p =2)
    center(q.A,q.b;sense,isnormalized,p)
end

function center(A::Matrix{<:Real} ,b::Vector{<:Real};sense=nothing, isnormalized = false, p =2)
    n,m = size(A)
    bl = fill(-1e30,m)
    sense = isnothing(sense) ? zeros(Int32,m+1) : sense
    rA = isnormalized ? ones(m) : [norm(a,p) for a in eachcol(A)]

    d = DAQP.Model();
    DAQP.settings(d,Dict(:primal_tol=>1e-8))
    DAQP.setup(d,zeros(0,0),[-1; zeros(n);],[rA';A],[1e30;b],[0;bl],sense;A_rowmaj=true);
    x,fval,exitflag,info = DAQP.solve(d);
    return exitflag ==1 ? (x[2:end], x[1]) :  (NaN, -Inf)
end
