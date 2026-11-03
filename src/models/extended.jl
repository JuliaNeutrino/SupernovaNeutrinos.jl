# A model extended beyond the end of its simulation by a cooling tail.

"""
    ExtendedModel(base; k = -1, A = nothing, tau_c = 36, alpha = 2.66)

`base` model followed, after its last time point t_end, by a cooling tail
with the shape of its last spectrum and the time dependence
`A tᵏ exp(-(t/τc)^α)` (t and τc in s), a generalization of Eq. 2 of Li,
Roberts and Beacom, PRD 103, 023016 (2021). By default `A` makes the tail
continuous at t_end.
"""
struct ExtendedModel{M <: SupernovaModel} <: SupernovaModel
    base::M
    k::Float64
    A::Float64
    tau_c::Float64
    alpha::Float64
end

function ExtendedModel(base::SupernovaModel; k::Real = -1.0, A::Union{Real, Nothing} = nothing, tau_c::Real = 36.0, alpha::Real = 2.66)
    t_end = last(times(base))
    A = isnothing(A) ? 1 / (t_end^k * exp(-(t_end / tau_c)^alpha)) : A
    return ExtendedModel(base, float(k), float(A), float(tau_c), float(alpha))
end

times(m::ExtendedModel) = times(m.base)
metadata(m::ExtendedModel) = metadata(m.base)

"Time dependence of the cooling tail at times `t` (s)."
function extended_time_dependence(m::ExtendedModel, t)
    t = _vector(t)
    any(<(0.5), t) && @warn "the cooling tail does not apply to times before 0.5 s"
    return @. m.A * t^m.k * exp(-(t / m.tau_c)^m.alpha)
end

"`true` for the times (s) in the cooling tail, after the end of the base model."
is_extended_tail(m::ExtendedModel, t) = _vector(t) .> last(times(m.base))

function spectra(m::ExtendedModel, t::AbstractVector, E::AbstractVector)
    t_end = last(times(m.base))
    tail = is_extended_tail(m, t)
    S = zeros(length(ThreeFlavor), length(t), length(E))
    any(.!tail) && (S[:, .!tail, :] = spectra(m.base, t[.!tail], E))
    if any(tail)
        # the shape of the spectrum at the end of the simulation
        last_spectrum = spectra(m.base, [t_end], E)
        f = extended_time_dependence(m, t[tail])
        S[:, tail, :] = last_spectrum .* reshape(f, 1, :, 1)
    end
    return S
end
