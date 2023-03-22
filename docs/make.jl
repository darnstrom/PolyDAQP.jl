using Documenter, PolyDAQP
push!(LOAD_PATH,"../src/")
makedocs(sitename="PolyDAQP.jl",
         pages = [
                  "Home" => "index.md",
                  "API" => "api.md",
                 ]
        )
