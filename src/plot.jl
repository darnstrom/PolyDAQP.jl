## Plots.jl
using RecipesBase
@recipe function f(p::Polyhedron)
    seriestype --> :shape
    legend --> false
    size(p.A,1) != 2 && error("Only plotting for two-dimensional Polyhedron is supported")
    vs=vrep_2d(minrep(p))
    isempty(vs) && return nothing
    return [first.(vs);first(vs[1])], [last.(vs);last(vs[1])]
end

@recipe function f(ps::Vector{Polyhedron})
    fill_z = get(plotattributes, :fill_z, nothing)
    if(!isnothing(fill_z))
        clims --> extrema(fill_z).+ (-0.5,0.5)
        fill_z  = (fill_z isa AbstractVector) ? fill_z : [fill_z for _ in 1:length(ps)]
    end
    for (k,p) in enumerate(ps)
        size(p.A,1) != 2 && error("Only plotting for two-dimensional Polyhedron is supported")
        vs=vrep_2d(minrep(p))
        isempty(vs) && continue
        @series begin
            seriestype --> :shape
            legend --> false
            if !isnothing(fill_z)
                fill_z := fill_z[k]
            end
            [first.(vs);first(vs[1])], [last.(vs);last(vs[1])]
        end
    end
end

@recipe function f(ps::Vector{<:Tuple{Polyhedron,Function}})
    for (k,(p,fz)) in enumerate(ps)
        size(p.A,1) != 2 && error("Only plotting for two-dimensional Polyhedron is supported")
        vs=vrep_2d(minrep(p))
        nv = length(vs)
        nv == 0  && continue
        z = fz.(vs)
        @series begin
            st --> :mesh3d
            legend --> false
            connections--> (zeros(Int, nv-2),collect(1:nv-2),collect(2:nv-1))
            (first.(vs), last.(vs),z)
        end
    end
end


## PGFPlots
using PGFPlotsX
function vertex_connections(vss; cs=nothing)
    cs = isnothing(cs) ? collect(1:length(vss)) : cs
    nvmax = maximum(length(vs) for vs in vss)
    offset = Int(0);
    patches = zeros(0,nvmax+1); # vertex ids + color value
    vtable = zeros(0,3); # x,y,z for vertices
    for (vs,c) in zip(vss,cs)
        nv = length(vs)
        nv == 0 && continue
        patch = [fill(offset,nvmax-nv); collect(0:nv-1).+offset; c] # First part is padding
        patches = vcat(patches,patch')

        zcoords = length(vs[1]) > 2 ? getindex.(vs,3) : zeros(nv)
        vtable = [vtable; getindex.(vs,1) getindex.(vs,2) zcoords] 

        offset+=nv
    end
    return vtable, patches
end


function categorical_cbar(vals;clabel=nothing)
    max_val = maximum(vals);
    min_val = minimum(vals);
    cbar_ticks = collect(min_val:Int64(ceil((max_val-min_val)/6)):max_val)
    cbar_ticks .+= Int64(floor((max_val-cbar_ticks[end])/2));
    @pgf opts = {
                 point_meta_min=min_val-0.5,
                 point_meta_max=max_val+0.5,
                 colorbar,
                 colorbar_horizontal,
                 colorbar_sampled,
                 colorbar_style =
                 {
                  samples=max_val-min_val+2,
                  xlabel=clabel,
                  xticklabel_style={yshift="13pt"},
                  xtick_style={draw="none"},
                  xtick=cbar_ticks
                 },
                }
end

function pgfplot(vtable::Matrix{<:Real},patches::Matrix{<:Real})
    @pgf  Plot3(
                {
                 patch,
                 opacity=1.0,
                 line_width="0.5pt",
                 faceted_color="black",
                 "patch type" = "polygon",
                 "vertex count" = size(patches,2)-1,
                 "table/row sep" = "\\\\",
                 patch_table_with_point_meta = TableData(patches)
                },
                Table(:x => vtable[:,1], :y => vtable[:,2],:z => vtable[:,3])
               )
end

function pplot(ps::Vector{Polyhedron};cs=nothing, fs=nothing, opts=Dict{Symbol,Any}())
    vss = [vrep_2d(minrep(p.A,p.b)...) for p in ps]
    if(!isnothing(fs) && !isempty(fs))
        for (vs,f) in zip(vss,fs)
            for v in vs push!(v,f(v)) end
        end
    end
    vtable,patches = vertex_connections(vss;cs)

    if(raw"\usepgfplotslibrary{patchplots}" ∉ PGFPlotsX.CUSTOM_PREAMBLE)
        push!(PGFPlotsX.CUSTOM_PREAMBLE,raw"\usepgfplotslibrary{patchplots}")
    end
    pl = pgfplot(vtable,patches)
    @pgf Aopts ={
                }
    if isnothing(fs) || isempty(fs)
        push!(Aopts,:view=>(0,90))
    else
        push!(Aopts,
              :xmajorgrids=>"true",
              :ymajorgrids=>"true",
              :zmajorgrids=>"true")
    end
    if eltype(cs) <: Integer
        merge!(Aopts, categorical_cbar(cs))
    end

    for opt in opts
        push!(Aopts,opt)
    end
    Axis(Aopts,pl)
end
