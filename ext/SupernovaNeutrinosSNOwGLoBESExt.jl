module SupernovaNeutrinosSNOwGLoBESExt

# Event rates with SNOwGLoBES.jl: a fluence Container becomes a
# SNOwGLoBES.Fluence, so that the models and flavor transformations of this
# package feed SNOwGLoBES detectors directly.

using SupernovaNeutrinos
using SupernovaNeutrinos: SupernovaNeutrinosError, ContainerUnit
import SNOwGLoBES

# SNOwGLoBES column order: νe, νμ, ντ, ν̄e, ν̄μ, ν̄τ
const COLUMNS = (:e, :mu, :tau, :e_bar, :mu_bar, :tau_bar)

"""
    SNOwGLoBES.Fluence(fluence::Container; time_bin = 1)

The SNOwGLoBES fluence of a fluence `Container` (1/(MeV cm²), e.g.
`integrate(flux(model, times(model), E), :time)`) in one of its time bins.
"""
function SNOwGLoBES.Fluence(c::Container; time_bin::Integer = 1)
    c.unit == ContainerUnit(-1, 0, -2) ||
        throw(SupernovaNeutrinosError("a SNOwGLoBES fluence needs a fluence (1/(MeV cm²)), not [$(c.unit)]; integrate the flux over time"))
    can_integrate(c, :energy) || throw(SupernovaNeutrinosError("the energies of the fluence must be points, not bins"))
    c3 = c.scheme == ThreeFlavor ? c : ThreeFlavor * c
    spectra = reduce(hcat, [vec(c3[f].array[1, time_bin, :]) for f in COLUMNS])
    return SNOwGLoBES.Fluence(c.energy, spectra)
end

"""
    event_rates(detector::SNOwGLoBES.SupernovaDetector, fluence::Container)

SNOwGLoBES events for a fluence `Container`: a `SupernovaRates`, or a vector of
them, one per time bin, when the fluence has several.
"""
function SNOwGLoBES.event_rates(d::SNOwGLoBES.SupernovaDetector, c::Container)
    n = size(c, 2)
    n == 1 && return SNOwGLoBES.event_rates(d, SNOwGLoBES.Fluence(c))
    return SNOwGLoBES.event_rates(d, [SNOwGLoBES.Fluence(c; time_bin = k) for k in 1:n])
end

end
