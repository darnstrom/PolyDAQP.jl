## Bounding box
function bounding_box(A,b;sense=nothing)
    n,m = size(A)
    bl = fill(-1e30,m)
    sense = isnothing(sense) ? zeros(Int32,m) : sense
    ub,lb = zeros(n),zeros(n)
    for i = 1:n 
        f = zeros(n)
        for (fi,bound) in [(1,lb),(-1,ub)]
            f[i] = fi 
            d = DAQP.Model();
            DAQP.setup(d,zeros(0,0),f,A,b,bl,sense;A_rowmaj=true);
            x,fval,exitflag,info = DAQP.solve(d);
            bound[i] = (exitflag==1) ? x[i] : -fi*Inf 
        end
    end
    return ub,lb 
end
