@testset "Fornax models" begin
    @test length(param_combinations(Fornax_2022)) == 100
    @test only(data_files(Fornax_2022, (progenitor_mass = 12.0,))) == "lum_spec_12.00.bh_dat.h5"
    @test only(data_files(Fornax_2024, (progenitor_mass = 9.6,))) == "lum_spec_z9.6_dat.h5"
end

@testset "spherical harmonics" begin
    # orthonormality on the sphere
    nθ, nϕ = 200, 400
    θs = (0.5:nθ) .* (π / nθ)
    ϕs = (0.5:nϕ) .* (2π / nϕ)
    lm = [(l, m) for l in 0:2 for m in -l:l]
    G = [sum(SupernovaNeutrinos.real_spherical_harmonic(a..., θ, ϕ) * SupernovaNeutrinos.real_spherical_harmonic(b..., θ, ϕ) * sin(θ)
             for θ in θs, ϕ in ϕs) * (π / nθ) * (2π / nϕ) for a in lm, b in lm]
    @test G ≈ I atol = 1e-3
end
