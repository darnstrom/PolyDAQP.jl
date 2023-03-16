## Minimal representation 
# Ar,br = minrep(A,b)
# remove constraints for the polyhedron P = {x : A' x ≤ b} 
# such that {x : Ar' x ≤ br}  = {x : A' x ≤ b} 
function minrep(A::Matrix{<: Real},b::Vector{<: Real};sense=[],max_radius=1e30, check_unique=true, tol_weak=0)

    # Setup DAQP workspace 
    nth,m = size(A)  
    poly_isempty(A,b) && return zeros(nth,0), zeros(0)
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

