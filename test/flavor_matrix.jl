@testset "flavor matrices" begin
    M = FlavorMatrix((β, α) -> 10 * β.index + α.index, ThreeFlavor)
    @test M[:e, :mu] == 13
    @test M[:NU, :NU] == [11 13 15; 31 33 35; 51 53 55]
    @test transpose(M)[:e, :mu] == 31
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError FlavorMatrix(zeros(3, 3), ThreeFlavor)
    Id = FlavorMatrix(I, ThreeFlavor)
    @test Id * M == M
    # broadcasting over energy
    A = FlavorMatrix(I, ThreeFlavor; ne = 4)
    @test size(A * M) == (6, 6, 1, 4)
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError FlavorMatrix(I, TwoFlavor) * M

    # conversion matrices between schemes
    C = conversion_matrix(TwoFlavor, ThreeFlavor)
    @test C[:x, :mu] == 0.5 && C[:x, :tau] == 0.5 && C[:e, :e] == 1 && C[:x, :e] == 0
    @test conversion_matrix(ThreeFlavor, TwoFlavor)[:mu, :x] == 1
    @test conversion_matrix(ThreeFlavor, TwoFlavor)[:mu, :x_bar] == 0
    P = conversion_matrix(ThreeFlavorFlavorBasis, ThreeFlavor)
    @test P * conversion_matrix(ThreeFlavor, ThreeFlavorFlavorBasis) == FlavorMatrix(I, ThreeFlavorFlavorBasis)
    @test (TwoFlavor * Id * TwoFlavor) ≈ FlavorMatrix(I, TwoFlavor)
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError conversion_matrix(ThreeFlavorMassBasis, ThreeFlavor)
end
