function eliminate(A,b,ids;tol_weak=0)
    ids = sort(ids,rev = true) #TODO: allow for other ordering
    n0,m0 = size(A)
    level=1
    H = falses(m0,m0)  
    for i in 1:m0 H[i,i] = true end

    for id in ids
        I0,Ip,In = Int64[],Int64[],Int64[]
        for i in 1:length(b)
            if(A[id,i] > 0)
                push!(Ip,i)
            elseif(A[id,i] < 0)
                push!(In,i)
            else
                push!(I0,i)
            end
        end

        mask = trues(n0-level+1); mask[id] = false;
        An,bn = A[mask,I0], b[I0]
        Hn = H[:,I0] 
        for jp in Ip
            for jn in In
                Hcand = H[:,jp] .| H[:,jn]
                sum(Hcand) > level+1 && continue
                An = hcat(An, A[id,jp]*A[mask,jn] - A[id,jn]*A[mask,jp])
                bn = vcat(bn, A[id,jp]*b[jn] - A[id,jn]*b[jp])
                #normalize
                nrm = norm(view(An,:,size(An,2)))
                An[:,end] ./= nrm; bn[end]/=nrm

                Hn = hcat(Hn,Hcand)
            end
        end
        A,b,nonred_ids = minrep(An,bn;tol_weak,return_ids=true)
        #A,b = minrep(An,bn;tol_weak=0)
        H = Hn[:,nonred_ids]
        level+=1
    end
    return A,b
end
