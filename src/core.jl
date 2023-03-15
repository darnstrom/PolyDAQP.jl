## Chebyshev center
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
## Minimal representation 
# Ar,br = minrep(A,b)
# remove constraints for the polyhedron P = {x : A' x ≤ b} 
# such that {x : Ar' x ≤ br}  = {x : A' x ≤ b} 
function minrep(A::Matrix{<: Real},b::Vector{<: Real};sense=[],max_radius=1e30, check_unique=true, tol_weak=0)

    # Setup DAQP workspace 
    nth,m = size(A)  
    ms = length(b)-m;
    p=DAQP.setup_c_workspace(nth);
    blower = fill(-1e30,m);
    sense = isempty(sense) ? zeros(Cint,m) : sense

    DAQP.init_c_workspace_ldp(p,A,b,blower,sense;max_radius)
    unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,3)),m) # set m 
    unsafe_store!(Ptr{Cint}(p+fieldoffset(DAQP.Workspace,4)),ms) # set ms 

    # Produce minimal representation
    A,b = minrep(p;check_unique,tol_weak)

    # Free DAQP workspace
    DAQP.free_c_workspace(p);
    return A,b
end

# Minrep interal DAQP 
function minrep(p::Ptr{Cvoid}; check_unique=true, tol_weak=0)


    daqp_ws = unsafe_load(Ptr{DAQP.Workspace}(p))
    m,n = daqp_ws.m, daqp_ws.n
    sense = unsafe_wrap(Vector{Cint}, daqp_ws.sense, m, own=false)
    A = unsafe_wrap(Matrix{Cdouble}, daqp_ws.M, (n,m), own=false)
    b= unsafe_wrap(Vector{Cdouble}, daqp_ws.dupper, m, own=false)

    # Start finding redundant constraints
    is_redundant = -ones(Cint,m); 
    for i = 1:m
        (is_redundant[i]!= -1 || sense[i]&4 != 0) && continue; # Decided from previous iteration 

        ccall((:reset_daqp_workspace,DAQP.libdaqp),Cvoid,(Ptr{Cvoid},),p);

        sense[i] = 5; # Force ith constraint to equality
        b[i]+=tol_weak; # displace  
        test =ccall((:add_constraint,DAQP.libdaqp),Cint,(Ptr{Cvoid},Cint,Float64),p,i-1,1.0)

        # Check if system is infeasible (infeasible => reudandant) 
        exitflag =ccall((:daqp_ldp,DAQP.libdaqp), Int32, (Ptr{Cvoid},),p);
        #ws.nLPs +=1
        daqp_ws = unsafe_load(Ptr{DAQP.Workspace}(p))
        AS = unsafe_wrap(Vector{Cint}, daqp_ws.WS, daqp_ws.n_active, own=false)
        if(exitflag==-1)
            is_redundant[i]=1
        else
            is_redundant[i]=0
            sense[i] &=~4; # Should be modifiable in later iterations 
            if(exitflag==1) # All activate constraints must also be nonredundant 
                is_redundant[AS.+1] .=0; 
            end
        end
        b[i] -= tol_weak;# restore 
        sense[AS.+1].&=~1; # Deactivate TODO: make sure equality constraint are not deactivated 
    end

    # Check for unique
    if(check_unique)
        tol = 1e-10
        for i = 1:m
            is_redundant[i] == 1 && continue
            for j = i+1:m
                is_redundant[j] == 1 && continue
                abs(b[i]-b[j]) > tol && continue 
                abs(A[1,i]-A[1,j]) > tol && continue 
                if(norm(view(A,:,i)-view(A,:,j)) < tol) 
                    is_redundant[j] = 1
                end
            end
        end
    end
    sense[1:m] .= 0 # reset sense
    nonred_ids = findall(is_redundant.==0)
    return A[:,nonred_ids],b[nonred_ids] 
end

## Slice
function slice(A,b,ids;values = nothing)
    n,m = size(A)
    ids_keep = setdiff(collect(1:n),ids)
    values = isnothing(values) ? zeros(length(ids)) : values
    rows = Int[]
    for i = 1:m 
        if(norm(A[ids_keep,i])>1e-12)
            push!(rows,i)
        elseif(b[i] < 0)
            @error "Polyhedron is empty"
        end
    end
    return A[ids_keep,rows], b[rows]-A[ids,rows]'*values
end

## Containment
function iscontained(x,A,b;tol=1e-6)
    n,m = size(A)
    for i = 1:m
        (A[:,i]'*x > b[i]+tol) && return false
    end
    return true 
end

## Naive vertex enumeration for 2D
function vrep_2d(A,b)
    n,m = size(A)
    vs = Vector{Float64}[]
    visited,i = falses(m),1
    visited[i] = true
    while(!all(visited))
        for j = 1:m 
            visited[j] && continue
            L = lu!([A[:,i]'; A[:,j]'],check=false)
            r = [b[i];b[j]]
            !issuccess(L) && continue
            v = ldiv!(L,r)
            if(iscontained(v,A,b))
                push!(vs,v)
                i, visited[j] = j, true
                break
            end
        end
    end
    # The last point has to connect to the first facet 
    L = lu!([A[:,i]'; A[:,1]'],check=false)
    r = [b[i];b[1]]
    v = ldiv!(L,r)
    if(iscontained(v,A,b))
        push!(vs,v)
    else
        @error "Could not form a cycle"
    end
    return vs
end
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
