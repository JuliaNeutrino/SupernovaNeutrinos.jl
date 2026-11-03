@testset "containers" begin
    tc = collect(range(0, 2, length = 21))
    Ec = collect(range(0, 10, length = 11))
    # a(f) · t · E² for each flavor
    A = [f * ti * Ei^2 for f in 1:6, ti in tc, Ei in Ec]
    c = Container(A, collect(ThreeFlavor), tc, Ec; unit = SupernovaNeutrinos.FLUX_UNIT)
    @test can_integrate(c, :time) && can_integrate(c, :energy) && can_sum(c, :flavor)
    @test c[:e_bar].array == A[2:2, :, :]
    @test size(c[:NU]) == (3, 21, 11)

    ct = integrate(c, :time)
    @test size(ct) == (6, 1, 11) && ct.time == [0, 2]
    @test ct.unit == SupernovaNeutrinos.ContainerUnit(-1, 0, -2)
    @test ct.array[3, 1, :] ≈ 3 .* 2 .* Ec .^ 2           # ∫₀² t dt = 2
    @test can_sum(ct, :time) && !can_integrate(ct, :time)
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError integrate(ct, :time)
    @test sum(ct, :time).array == ct.array

    # bins, and limits beyond the samples (zero outside)
    cb = integrate(c, :time; limits = [0, 0.5, 1, 3])
    @test cb.array[1, :, 5] ≈ [0.125, 0.375, 1.5] .* 16
    # partial bins are exact for a linear function
    cp = integrate(c, :time; limits = [0.05, 0.37])
    @test cp.array[1, 1, 3] ≈ (0.37^2 - 0.05^2) / 2 * 4

    ce = integrate(ct, :energy)
    @test ce.array[:, 1, 1] ≈ (1:6) .* 2 .* (1000 / 3 + 5 / 3)   # trapezoid rule of E² on 1 MeV steps
    @test ce.unit == SupernovaNeutrinos.ContainerUnit(0, 0, -2)
    total = sum(ce, :flavor)
    @test size(total) == (1, 1, 1) && only(total.array) ≈ sum(ce.array)
    @test integrate_or_sum(ct, :time) == sum(ct, :time)

    # flavor conversion: TwoFlavor x is the average of μ and τ
    c2 = TwoFlavor * c
    @test c2.scheme == TwoFlavor && c2[:x].array ≈ (A[3:3, :, :] .+ A[5:5, :, :]) ./ 2
    @test (2c).array == 2A
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError Container(zeros(6, 3, 3), collect(ThreeFlavor), tc, Ec; unit = SupernovaNeutrinos.FLUX_UNIT)
    @test occursin("flux", sprint(show, MIME"text/plain"(), c))
end
