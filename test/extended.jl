@testset "extended model" begin
    m = PinchedModel([0.0, 1.0, 2.0]; luminosity = Dict("NU_E" => [1e52, 2e52, 1e52], "NU_E_BAR" => fill(1e52, 3), "NU_X" => fill(5e51, 3)),
                     mean_energy = Dict("NU_E" => fill(10.0, 3), "NU_E_BAR" => fill(14.0, 3), "NU_X" => fill(16.0, 3)),
                     pinch = Dict("NU_E" => fill(2.0, 3), "NU_E_BAR" => fill(2.0, 3), "NU_X" => fill(2.0, 3)))

    # continuous at the end of the simulation, then a cooling tail
    x = ExtendedModel(m)
    s = initial_spectra(x, [2.0, 3.0, 10.0], [10.0])
    @test s.array[:, 1, 1] ≈ initial_spectra(m, [2.0], [10.0]).array[:, 1, 1]
    @test s.array[1, 2, 1] / s.array[1, 1, 1] ≈ only(extended_time_dependence(x, 3.0))
    @test is_extended_tail(x, [1.0, 3.0]) == [false, true]
end
