const t = collect(0.0:9.0)
const E = collect(range(1, 100, length = 21))

two(P) = TwoFlavor * P * TwoFlavor
chain(sn, p; kwargs...) = TransformationChain(sn; mixing_params = p, kwargs...)

@testset "simple transformations" begin
    P = P_ff(NoTransformation(), t, E)
    @test P.flavor_in == ThreeFlavor && P == FlavorMatrix(I, ThreeFlavor)
    P = P_ff(CompleteExchange(), t, E)
    for β in ThreeFlavor, α in ThreeFlavor
        @test P[β, α] ≈ (is_neutrino(β) == is_neutrino(α)) * (β != α) * 0.5
    end
    P = ThreeFlavor * P_ff(ThreeFlavorDecoherence(), t, E) * ThreeFlavor
    for β in ThreeFlavor, α in ThreeFlavor
        @test P[β, α] ≈ (is_neutrino(β) == is_neutrino(α)) / 3
    end
end

@testset "MSW in the supernova" begin
    for (p, adiabatic) in ((MixingParameters(:normal), true), (MixingParameters(:inverted), true),
                           (MixingParameters(:normal), false), (MixingParameters(:inverted), false))
        th12, th13, _ = mixing_angles(p)
        P = two(P_ff(chain(adiabatic ? AdiabaticMSW() : NonAdiabaticMSWH(), p), t, E))
        # survival probabilities: SNEWPY test_04_xforms
        pee = (adiabatic && p.mass_order == NormalOrdering) ? sind(th13)^2 : (sind(th12) * cosd(th13))^2
        pbar = (adiabatic && p.mass_order == InvertedOrdering) ? sind(th13)^2 : (cosd(th12) * cosd(th13))^2
        @test P[:e, :e] ≈ pee
        @test P[:e, :x] ≈ 1 - pee
        @test P[:x, :x] ≈ (1 + pee) / 2
        @test P[:x, :e] ≈ (1 - pee) / 2
        @test P[:e_bar, :e_bar] ≈ pbar
        @test P[:e_bar, :x_bar] ≈ 1 - pbar
        @test P[:x_bar, :x_bar] ≈ (1 + pbar) / 2
        @test P[:x_bar, :e_bar] ≈ (1 - pbar) / 2
    end
end

@testset "two-flavor decoherence" begin
    for p in (MixingParameters(:normal), MixingParameters(:inverted))
        th12, th13, _ = mixing_angles(p)
        De1, De2, De3 = (cosd(th12) * cosd(th13))^2, (sind(th12) * cosd(th13))^2, sind(th13)^2
        P = two(P_ff(chain(TwoFlavorDecoherence(), p), t, E))
        if p.mass_order == NormalOrdering
            @test P[:e, :e] ≈ (De2 + De3) / 2
            @test P[:e_bar, :e_bar] ≈ De1
        else
            @test P[:e, :e] ≈ De2
            @test P[:e_bar, :e_bar] ≈ (De1 + De3) / 2
        end
    end
end

@testset "sterile neutrino" begin
    for o in (:normal, :inverted), adiabatic in (true, false)
        p = FourFlavorMixingParameters(MixingParameters(o); theta14 = 1)
        th12, th13, _, th14 = mixing_angles(p)
        De = ((cosd(th12) * cosd(th13) * cosd(th14))^2, (sind(th12) * cosd(th13) * cosd(th14))^2,
              (sind(th13) * cosd(th14))^2, sind(th14)^2)
        Ds = ((cosd(th12) * cosd(th13) * sind(th14))^2, (sind(th12) * cosd(th13) * sind(th14))^2,
              (sind(th13) * sind(th14))^2, cosd(th14)^2)
        P = two(P_ff(chain(adiabatic ? AdiabaticMSWes() : NonAdiabaticMSWes(), p), t, E))
        # electron flavor: matter state of νe and of ν̄e, and the two states of νx (SNEWPY test_04_xforms)
        (ie, ix, ibar, ixbar) = if adiabatic
            o == :normal ? (4, (1, 2), 1, (3, 4)) : (4, (1, 3), 3, (2, 4))
        else
            o == :normal ? (3, (1, 2), 1, (2, 3)) : (2, (1, 3), 3, (1, 2))
        end
        @test P[:e, :e] ≈ De[ie]
        @test P[:e, :x] ≈ De[ix[1]] + De[ix[2]]
        @test P[:x, :x] ≈ (2 - De[ix[1]] - De[ix[2]] - Ds[ix[1]] - Ds[ix[2]]) / 2
        @test P[:x, :e] ≈ (1 - De[ie] - Ds[ie]) / 2 atol = 1e-12
        @test P[:e_bar, :e_bar] ≈ De[ibar]
        @test P[:e_bar, :x_bar] ≈ De[ixbar[1]] + De[ixbar[2]]
        @test P[:x_bar, :e_bar] ≈ (1 - De[ibar] - Ds[ibar]) / 2 atol = 1e-12
    end
    P = P_ff(FourFlavorDecoherence(), t, E)
    @test P.flavor_in == FourFlavorFlavorBasis && size(P, 1) == 8
    @test_throws MethodError P_ff(chain(AdiabaticMSWes(), MixingParameters()), t, E)
end

@testset "neutrino decay" begin
    x = NeutrinoDecay(mass = 1.0, lifetime = 86400.0, distance = 10.0)
    # exp(-m c d / (E τ)) with m in eV, E in MeV
    d = @. exp(-(1.0 / (E * 1e6)) * 10 * 3.0856775814913673e21 / (2.99792458e10 * 86400))
    @test SupernovaNeutrinos.decay_factor(x, E) ≈ d
    for p in (MixingParameters(:normal), MixingParameters(:inverted))
        th12, th13, _ = mixing_angles(p)
        De1, De2, De3 = (cosd(th12) * cosd(th13))^2, (sind(th12) * cosd(th13))^2, sind(th13)^2
        P = two(P_ff(TransformationChain(AdiabaticMSW(), x; mixing_params = p), t, E))
        if p.mass_order == NormalOrdering
            pee = @. De1 * (1 - d) + De3 * d
            pexbar = @. De1 * (1 - d) + De2 + De3 * d
        else
            pee = @. De2 * d + De3 * (1 - d)
            pexbar = @. De1 + De2 * d + De3 * (1 - d)
            @test vec(P[:e_bar, :e_bar]) ≈ fill(De3, length(E))
        end
        @test vec(P[:e, :e]) ≈ pee
        @test vec(P[:x, :e]) ≈ (1 .- pee) ./ 2
        @test vec(P[:e_bar, :x_bar]) ≈ pexbar
        @test vec(P[:x_bar, :x_bar]) ≈ 1 .- pexbar ./ 2
    end
end

@testset "quantum decoherence" begin
    x = QuantumDecoherence(Gamma3 = 1e-27, Gamma8 = 2e-27, distance = 10.0, n = 1, E0 = 10.0)
    p = MixingParameters()
    Pmm = P_mm(x, p, t, E)
    g(Γ) = @. Γ * (E / 10) * 10 * 3.0856775814913673e21 / 1.973269804e-5
    @test vec(Pmm[:NU_1, :NU_1]) ≈ @. 1 / 3 + exp(-$g(1e-27) - $g(2e-27) / 3) / 2 + exp(-$g(2e-27)) / 6
    @test vec(Pmm[:NU_3, :NU_1]) ≈ @. 1 / 3 - exp(-$g(2e-27)) / 3
    @test vec(Pmm[:NU_3, :NU_3]) ≈ @. 1 / 3 + 2exp(-$g(2e-27)) / 3
    @test Pmm[:NU_BAR, :NU_BAR] == Pmm[:NU, :NU]
    @test sum(Pmm.array; dims = 1) ≈ ones(1, 6, 1, length(E))
    th12, th13, _ = mixing_angles(p)
    De = ((cosd(th12) * cosd(th13))^2, (sind(th12) * cosd(th13))^2, sind(th13)^2)
    P = two(P_ff(TransformationChain(AdiabaticMSW(), x; mixing_params = p), t, E))
    pee = sum(vec(Pmm[Symbol("NU_3"), Symbol("NU_$i")]) .* De[i] for i in 1:3)
    @test vec(P[:e, :e]) ≈ pee
end

@testset "transformation chains" begin
    c = TransformationChain(AdiabaticMSW())
    @test c.mixing_params == MixingParameters() && c.in_vacuum isa NoVacuumTransformation && c.in_earth isa NoEarthMatter
    @test TransformationChain(AdiabaticMSWes()).mixing_params isa FourFlavorMixingParameters
    @test P_ff(AdiabaticMSW(), t, E) == P_ff(c, t, E)
    @test string(c) == "AdiabaticMSW+NoVacuumTransformation+NoEarthMatter_NormalOrdering"
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError TransformationChain(NeutrinoDecay())
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError TransformationChain(AdiabaticMSW(), NeutrinoDecay(), NeutrinoDecay())
    for x in (AdiabaticMSW(), NonAdiabaticMSWH(), TwoFlavorDecoherence(), ThreeFlavorDecoherence(), AdiabaticMSWes(),
              NonAdiabaticMSWes(), FourFlavorDecoherence(), CompleteExchange(), NoTransformation(),
              TransformationChain(AdiabaticMSW(), NeutrinoDecay()), TransformationChain(AdiabaticMSW(), QuantumDecoherence()))
        P = P_ff(x, t, E)
        @test all(≈(1), sum(P.array; dims = 1))
    end
end
