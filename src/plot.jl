using RecipesBase
@recipe function f(p::Polyhedron)
    seriestype --> :shape
    legend --> false
    #A,b = minrep(p.A,p.b)
    size(p.A,1) != 2 && error("Only plotting for two-dimensional Polyhedron is supported")
    vs=vrep_2d(minrep(p.A,p.b)...)
    isempty(vs) && return nothing
    return [first.(vs);first(vs[1])], [last.(vs);last(vs[1])]
end

@recipe function f(ps::Vector{Polyhedron};z=nothing)
    if (!isnothing(z) && length(z) != length(ps)) 
        error("z has incorrect length: $(length(z)) ≂̸ $(length(ps))")
    end
    for (k,p) in enumerate(ps)
        size(p.A,1) != 2 && error("Only plotting for two-dimensional Polyhedron is supported")
        vs=vrep_2d(minrep(p.A,p.b)...)
        isempty(vs) && continue
        @series begin
            seriestype --> :shape
            legend --> false
            if(!isnothing(z))
                fill_z --> z[k]
            end
            [first.(vs);first(vs[1])], [last.(vs);last(vs[1])]
        end
    end
end
