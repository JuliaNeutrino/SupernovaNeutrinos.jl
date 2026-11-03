@testset "mixing parameters" begin
    p = MixingParameters()
    @test p.mass_order == NormalOrdering && p.theta12 == 33.68 && p.dm31_2 == 2.513e-3
    @test p.dm32_2 ≈ p.dm31_2 - p.dm21_2
    q = MixingParameters(:inverted; version = "NuFIT5.0")
    @test q.mass_order == InvertedOrdering && q.dm32_2 == -2.498e-3 && q.dm31_2 ≈ q.dm32_2 + q.dm21_2
    @test MixingParameters("INVERTED") == MixingParameters(InvertedOrdering)
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError MixingParameters(; version = "NuFIT1.0")
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError ThreeFlavorMixingParameters(; theta12 = 33, theta13 = 8,
        theta23 = 45, deltaCP = 0, dm21_2 = 7e-5, dm31_2 = 2e-3, mass_order = :inverted)
    @test ThreeFlavorMixingParameters(p; theta23 = 45).theta23 == 45
    flipped = ThreeFlavorMixingParameters(p; dm32_2 = -2.5e-3)
    @test flipped.mass_order == InvertedOrdering && flipped.dm31_2 ≈ -2.5e-3 + p.dm21_2
    @test mixing_angles(p) == (33.68, 8.56, 43.3)
    for v in ("NuFIT5.0", "NuFIT5.2", "NuFIT6.0", "PDG2022", "PDG2024"), o in (:normal, :inverted)
        @test MixingParameters(o; version = v).mass_order == SupernovaNeutrinos.MassOrdering(o)
    end

    p4 = FourFlavorMixingParameters(p; theta14 = 1, dm41_2 = 1.0)
    @test p4.theta12 == p.theta12 && p4.mass_order == NormalOrdering
    @test p4.dm42_2 ≈ 1.0 - p.dm21_2 && p4.dm43_2 ≈ 1.0 - p.dm31_2
    @test length(mixing_angles(p4)) == 6

    for q in (p, q, p4)
        U = vacuum_mixing_matrix(q).array[:, :, 1, 1]
        @test U * U' ≈ I
    end
    # PDG parametrization: Ue3 = sin θ13 e^{-iδ}
    U = vacuum_mixing_matrix(p)
    @test U[:e, :NU_3] ≈ sind(p.theta13) * cis(-deg2rad(p.deltaCP))
    @test U[:e_bar, :NU_3_BAR] ≈ conj(U[:e, :NU_3])
end

@testset "matter states at high density" begin
    for p in (MixingParameters(:normal), MixingParameters(:inverted), MixingParameters(:normal; version = "PDG2024"))
        P = pmf_high_density_limit(p)
        @test vec(sum(P.array; dims = 1)) ≈ ones(6)
        @test vec(sum(P.array; dims = 2)) ≈ ones(6)
        # the μ–τ content of the matter states is that of the eigenvectors of the μ–τ block
        H = SupernovaNeutrinos._vacuum_hamiltonian(p)
        for bar in (false, true)
            k = bar ? (5:6) : (2:3)
            v = abs2.(eigen(Hermitian(H.array[k, k, 1, 1])).vectors)
            μ, τ = bar ? (:mu_bar, :tau_bar) : (:mu, :tau)
            states = if p.mass_order == NormalOrdering
                bar ? (:NU_2_BAR, :NU_3_BAR) : (:NU_1, :NU_2)
            else
                bar ? (:NU_1_BAR, :NU_2_BAR) : (:NU_3, :NU_1)
            end
            @test [P[states[1], μ] P[states[2], μ]; P[states[1], τ] P[states[2], τ]] ≈ v
        end
        e = p.mass_order == NormalOrdering ? (:NU_3, :NU_1_BAR) : (:NU_2, :NU_3_BAR)
        @test P[e[1], :e] == 1 && P[e[2], :e_bar] == 1
    end
end
