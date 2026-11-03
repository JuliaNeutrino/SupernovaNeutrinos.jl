@testset "flavor schemes" begin
    @test length(ThreeFlavor) == 6 && length(TwoFlavor) == 4 && length(FourFlavor) == 8
    @test ThreeFlavor.NU_E == ThreeFlavor[:e] == ThreeFlavor["nu_e"] == ThreeFlavor[1]
    @test Symbol(ThreeFlavor[:mu_bar]) == :NU_MU_BAR
    @test ThreeFlavor[:NU] == [ThreeFlavor.NU_E, ThreeFlavor.NU_MU, ThreeFlavor.NU_TAU]
    @test ThreeFlavor["NU_BAR"] == [ThreeFlavor.NU_E_BAR, ThreeFlavor.NU_MU_BAR, ThreeFlavor.NU_TAU_BAR]
    @test ThreeFlavor[[:e, :tau_bar]] == [ThreeFlavor.NU_E, ThreeFlavor.NU_TAU_BAR]
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError ThreeFlavor[:s]
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError ThreeFlavor[TwoFlavor.NU_X]
    f = ThreeFlavor.NU_MU_BAR
    @test is_antineutrino(f) && !is_neutrino(f) && is_muon(f) && lepton(f) == "MU"
    @test is_electron(ThreeFlavor.NU_E) && is_tauon(ThreeFlavor.NU_TAU) && is_sterile(FourFlavor.NU_S)
    @test to_tex(f) == "\$\\overline{\\nu}_{\\mu}\$"
    @test ThreeFlavor.NU_E != ThreeFlavorFlavorBasis.NU_E
    @test [Symbol(s) for s in ThreeFlavorMassBasis] == [:NU_1, :NU_2, :NU_3, :NU_1_BAR, :NU_2_BAR, :NU_3_BAR]
end
