@testset "model files" begin
    dir = mktempdir()
    file = joinpath(dir, "analytic.dat")
    write(file, """
    TIME L_NU_E L_NU_E_BAR L_NU_X E_NU_E E_NU_E_BAR E_NU_X ALPHA_NU_E ALPHA_NU_E_BAR ALPHA_NU_X
    0.0 1e52 1e52 1e52 10 12 14 2 2 2
    1.0 2e52 2e52 2e52 10 12 14 2 2 2
    """)
    m = Analytic3Species(file)
    @test times(m) == [0.0, 1.0]
    @test initial_spectra(m, [1.0], [12.0]).array[2, 1, 1] ≈ SupernovaNeutrinos.pinched_spectrum(2e52, 12, 2, 12)

    T = read_named_table(IOBuffer("# comment line\n# a b c\n1 2 3\n4 5 6\n"))
    @test collect(keys(T)) == ["a", "b", "c"] && T["b"] == [2.0, 5.0]
    @test read_table(IOBuffer("x y\nunits\n1 2.5D0\n3 4\n"); skip = 2) == [1 2.5; 3 4]

    # the Garching format: one file per flavor, with ⟨E²⟩ in the last column
    for (sfx, L) in (("nue", 1.0), ("nuebar", 2.0), ("nux", 3.0))
        open(joinpath(dir, "s1.0c_LS220_$sfx"), "w") do io
            println(io, "0.1 $L 10.0 120.0")
            sfx == "nux" || println(io, "0.2 $L 11.0 150.0")
        end
    end
    g = SupernovaNeutrinos.read_garching(joinpath(dir, "s1.0c"))
    @test times(g) == [0.1, 0.2]
    @test g.luminosity[:, 1] == [1, 2, 3, 3, 3, 3] .* 1e51
    @test g.pinch[1, 1] ≈ (2 * 100 - 120) / (120 - 100)
    @test g.luminosity[3, 2] == 0                          # missing νx row
end

@testset "model registry and core-collapse models" begin
    @test keys(parameters(Nakazato_2013)) == (:progenitor_mass, :revival_time, :metallicity, :eos)
    combos = param_combinations(Nakazato_2013)
    @test length(combos) == 24
    @test all(p -> p.eos == "shen" || p.revival_time == 0, combos)
    @test length(param_combinations(Warren_2020)) == 600
    @test length(param_combinations(Bugli_2021)) == 20
    @test param_combinations(Fischer_2020) == [(;)]
    for M in all_models(), p in param_combinations(M)
        nameof(M) == :Kato_2017 && continue               # its file list needs the light curve
        @test !isempty(data_files(M, p))
    end
    @test only(data_files(Nakazato_2013, (progenitor_mass = 13.0, revival_time = 100.0, metallicity = 0.02, eos = "shen"))) ==
          "nakazato-shen-z0.02-t_rev100ms-s13.0.fits"
    @test only(data_files(Warren_2020, (progenitor_mass = 10.0, turbmixing_param = 1.25))) ==
          "stir_a1.25/stir_multimessenger_a1.25_m10.0.h5"
    @test only(data_files(Warren_2020, (progenitor_mass = 10.25, turbmixing_param = 1.23))) ==
          "stir_a1.23/stir_multimessenger_a1.23_m10.25.h5"
    @test only(data_files(Takata_2025, (progenitor_mass = 11.2, axion_mass = 40, axion_coupling = 4))) == "11_040_04.dat"
    @test only(data_files(Kuroda_2020, (rotational_velocity = 1, magnetic_field_exponent = 12, progenitor_mass = 20.0))) == "LnuR10B12.dat"

    err = SupernovaNeutrinos.SupernovaNeutrinosError
    @test_throws err Bollig_2016(progenitor_mass = 12)
    @test_throws err Bollig_2016()
    @test_throws err Bollig_2016(progenitor_mass = 27, mass = 3)
    @test_throws err Nakazato_2013(progenitor_mass = 13, revival_time = 0, metallicity = 0.02, eos = "shen")
    @test SupernovaNeutrinos.validate(Bollig_2016; progenitor_mass = 11.2000001) == (progenitor_mass = 11.2, eos = "LS220")
end
