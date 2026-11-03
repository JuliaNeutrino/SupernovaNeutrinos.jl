using SupernovaNeutrinos
using Test
using LinearAlgebra

# Tests that download model files and the SNOwGLoBES data (with SNOwGLoBES.jl,
# which needs Julia 1.11); turn them off with
# SUPERNOVANEUTRINOS_TEST_NETWORK=false
const NETWORK = get(ENV, "SUPERNOVANEUTRINOS_TEST_NETWORK", "true") != "false"

@testset "SupernovaNeutrinos" begin
    include("flavors.jl")
    include("flavor_matrix.jl")
    include("mixing.jl")
    include("container.jl")
    include("transformations.jl")
    include("pinched.jl")
    include("extended.jl")
    include("ccsn.jl")
    include("fornax.jl")
    include("presn.jl")
    NETWORK && include("snowglobes.jl")
    include("aqua.jl")
end
