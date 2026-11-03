@testset "pre-supernova models" begin
    @test length(all_models()) == 24
    @test only(data_files(Odrzywolek_2010, (progenitor_mass = 15.0,))) == "s15_nuebar_data.txt"
    for M in (Odrzywolek_2010, Patton_2017, Yoshida_2016, Myers_2026)
        @test !isempty(param_combinations(M))
    end
end
