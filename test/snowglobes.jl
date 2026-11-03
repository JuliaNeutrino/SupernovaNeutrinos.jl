# Downloads a model and, through SNOwGLoBES.jl, the SNOwGLoBES data. Reference
# totals from SNEWPY's crosscheck (Bollig_2016 11.2 M☉, AdiabaticMSW with
# normal ordering, 10 kpc).
using SNOwGLoBES

@testset "event rates with SNOwGLoBES.jl" begin
    m = Bollig_2016(progenitor_mass = 11.2)
    fl = integrate(flux(m, times(m), range(0.5, 100, length = 201); transformation = AdiabaticMSW()), :time)
    f = SNOwGLoBES.Fluence(fl)
    @test f.energies == fl.energy
    @test f.nuebar == vec(fl[:e_bar].array) && f.numu == vec(fl[:mu].array)
    @test SNOwGLoBES.Fluence(TwoFlavor * fl).nutau ≈ (f.numu .+ f.nutau) ./ 2
    @test_throws SupernovaNeutrinos.SupernovaNeutrinosError SNOwGLoBES.Fluence(flux(m, [0.1], [10.0]))

    unsmeared = Dict("icecube" => 7.2352e6, "ar40kt" => 2858.7, "halo2" => 146.4671)
    smeared = Dict("icecube" => 320766.3, "ar40kt" => 2696.40, "halo2" => 53.945)
    for (d, mat) in (("icecube", :water), ("ar40kt", :argon), ("halo2", :lead))
        r = event_rates(SupernovaDetector(Symbol(d), mat), fl)
        @test total_events(r; smeared = false) ≈ unsmeared[d] rtol = 0.01
        @test total_events(r) ≈ smeared[d] rtol = 0.01
    end

    # several time bins give one result per bin, adding up to the whole
    bins = integrate(flux(m, times(m), range(0.5, 100, length = 201); transformation = AdiabaticMSW()), :time;
                     limits = [first(times(m)), 0.5, last(times(m))])
    det = SupernovaDetector(:halo2, :lead)
    rs = event_rates(det, bins)
    @test length(rs) == 2
    @test sum(total_events, rs) ≈ total_events(event_rates(det, fl)) rtol = 1e-6
end
