@testset "pinched model" begin
    # spectrum normalization: ∫ dN/dE dE = L/⟨E⟩ and ∫ E dN/dE dE = L
    L, Emean, α = 1e52, 12.0, 2.5
    Es = range(0, 200, length = 40001)
    S = SupernovaNeutrinos.pinched_spectrum.(L, Emean, α, Es)
    trapz(y) = sum(y[1:end-1] .+ y[2:end]) / 2 * step(Es)
    @test trapz(S) ≈ L / SupernovaNeutrinos.MEV_IN_ERG / Emean rtol = 1e-6
    @test trapz(S .* Es) ≈ L / SupernovaNeutrinos.MEV_IN_ERG rtol = 1e-6
    @test trapz(S .* Es .^ 2) / trapz(S) ≈ Emean^2 * (2 + α) / (1 + α) rtol = 1e-6
    @test SupernovaNeutrinos.pinched_spectrum(L, Emean, α, 0.0) == 0 && SupernovaNeutrinos.pinched_spectrum(0.0, Emean, α, 10.0) == 0

    time = [0.0, 1.0, 2.0]
    m = PinchedModel(time; luminosity = Dict("NU_E" => [1e52, 2e52, 1e52], "NU_E_BAR" => fill(1e52, 3), "NU_X" => fill(5e51, 3)),
                     mean_energy = Dict("NU_E" => fill(10.0, 3), "NU_E_BAR" => fill(14.0, 3), "NU_X" => fill(16.0, 3)),
                     pinch = Dict("NU_E" => fill(2.0, 3), "NU_E_BAR" => fill(2.0, 3), "NU_X" => fill(2.0, 3)))
    s = initial_spectra(m, [0.5, 3.0], [10.0])
    @test s.array[1, 1, 1] ≈ SupernovaNeutrinos.pinched_spectrum(1.5e52, 10, 2, 10)
    @test s.array[4, 1, 1] == s.array[6, 1, 1] ≈ SupernovaNeutrinos.pinched_spectrum(5e51, 16, 2, 10)
    @test all(iszero, s.array[:, 2, :])                 # after the end of the simulation
    @test times(m) == time

    # flux at 10 kpc
    F = flux(m, time, [10.0]; distance = 10)
    @test F.array[1, 2, 1] ≈ SupernovaNeutrinos.pinched_spectrum(2e52, 10, 2, 10) / (4π * (10 * SupernovaNeutrinos.KPC_IN_CM)^2)
    @test F.unit == SupernovaNeutrinos.FLUX_UNIT
    Fx = flux(m, time, [10.0]; transformation = CompleteExchange())
    @test Fx[:e].array ≈ (F[:mu].array .+ F[:tau].array) ./ 2

    # a restart of the simulation: the times step back and the later rows win
    m2 = PinchedModel(Dict("TIME" => [0.0, 1.0, 2.0, 1.5, 3.0], "L_NU_E" => [1, 1, 9, 2, 2] .* 1e52, "E_NU_E" => fill(10.0, 5),
        "ALPHA_NU_E" => fill(2.0, 5), "L_NU_E_BAR" => fill(1e52, 5), "E_NU_E_BAR" => fill(10.0, 5), "ALPHA_NU_E_BAR" => fill(2.0, 5),
        "L_NU_X" => fill(1e52, 5), "E_NU_X" => fill(10.0, 5), "ALPHA_NU_X" => fill(2.0, 5)))
    @test times(m2) == [0.0, 1.0, 1.5, 3.0]
end
