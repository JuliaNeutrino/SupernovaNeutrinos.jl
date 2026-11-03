using SupernovaNeutrinos
using Test
using LinearAlgebra

@testset "SupernovaNeutrinos" begin
    include("flavors.jl")
    include("flavor_matrix.jl")
    include("mixing.jl")
    include("container.jl")
    include("transformations.jl")
    include("pinched.jl")
    include("extended.jl")
    include("aqua.jl")
end
