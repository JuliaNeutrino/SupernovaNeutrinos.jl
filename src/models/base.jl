# Supernova models: neutrino spectra at the source as functions of time and
# energy, and the flux they give at a detector.

"""
    SupernovaModel

A model of the neutrino emission of a supernova. A model implements

- `spectra(model, t, E)`: the `(flavor, time, energy)` array of the emitted
  spectrum d²N/dE dt in 1/(MeV s) for the six flavors of
  [`ThreeFlavor`](@ref), at times `t` (s) and energies `E` (MeV);
- `times(model)`: the time points of the simulation (s);

and optionally `metadata(model)`.
"""
abstract type SupernovaModel end

"Time points (s) of the simulation."
times(m::SupernovaModel) = m.time

"Dictionary of the physics parameters of the model."
metadata(::SupernovaModel) = OrderedDict{String, Any}()

function spectra end

_vector(x::Real) = [float(x)]
_vector(x::AbstractVector) = collect(float.(x))

"""
    initial_spectra(model, t, E) -> Container

Spectrum d²N/dE dt (1/(MeV s)) emitted by the supernova, before any flavor
transformation, at times `t` (s) and energies `E` (MeV), scalars or vectors.
"""
function initial_spectra(m::SupernovaModel, t, E)
    t, E = _vector(t), _vector(E)
    return Container(spectra(m, t, E), collect(ThreeFlavor), t, E; unit = SPECTRUM_UNIT)
end

"""
    transformed_spectra(model, t, E, transformation) -> Container

[`initial_spectra`](@ref) after the flavor `transformation`.
"""
transformed_spectra(m::SupernovaModel, t, E, transformation) = apply(transformation, initial_spectra(m, t, E))

"""
    flux(model, t, E; distance = 10, transformation = NoTransformation()) -> Container

Neutrino flux d²Φ/dE dt in 1/(MeV s cm²) at `distance` (kpc) from the
supernova, after the flavor `transformation`, at times `t` (s) and energies
`E` (MeV).

```julia
model = Bollig_2016(progenitor_mass = 27)
F = flux(model, times(model), 0:0.5:100; transformation = AdiabaticMSW())
```
"""
function flux(m::SupernovaModel, t, E; distance::Real = 10.0, transformation = NoTransformation())
    c = transformed_spectra(m, t, E, transformation)
    area = 4π * (distance * KPC_IN_CM)^2
    return _with(c; array = c.array ./ area, unit = c.unit * ContainerUnit(0, 0, -2))
end

# Linear interpolation of (x, y) at xq, `outside` beyond the ends of x
function _interp(x::AbstractVector, y::AbstractVector, xq::Real; outside = 0.0)
    (xq < first(x) || xq > last(x)) && return outside
    k = searchsortedlast(x, xq)
    k >= length(x) && return y[end]
    x0, x1 = x[k], x[k+1]
    return x1 == x0 ? y[k] : y[k] + (y[k+1] - y[k]) * (xq - x0) / (x1 - x0)
end

# Linear interpolation clamped to the end values
_interp_clamped(x, y, xq) = xq <= first(x) ? y[1] : xq >= last(x) ? y[end] : _interp(x, y, xq)

"""
    pinched_spectrum(L, Emean, α, E)

Spectrum dN/dE dt (1/(MeV s)) of the pinched (Keil–Raffelt–Janka)
distribution with luminosity `L` (erg/s), mean energy `Emean` (MeV) and pinch
parameter `α` at energy `E` (MeV):

    L/⟨E⟩² · (1+α)^(1+α)/Γ(1+α) · (E/⟨E⟩)^α · exp(-(1+α) E/⟨E⟩)
"""
function pinched_spectrum(L::Real, Emean::Real, α::Real, E::Real)
    (E > 0 && L > 0 && Emean > 0 && α > -1) || return 0.0
    # in logarithms for numerical stability at large α
    f = exp(log(L / MEV_IN_ERG) - (2 + α) * log(Emean) + (1 + α) * log1p(α) - loggamma(1 + α) +
            α * log(E) - (1 + α) * E / Emean)
    return isfinite(f) ? f : 0.0
end

"""
    PinchedModel(time; luminosity, mean_energy, pinch)
    PinchedModel(columns::AbstractDict)

A model whose spectra follow the pinched distribution
([`pinched_spectrum`](@ref)), from the luminosity (erg/s), mean energy (MeV)
and pinch parameter of each flavor at the times `time` (s). The three
keywords are dictionaries `flavor => vector`, with flavors of
[`ThreeFlavor`](@ref) or the names `"NU_E"`, `"NU_E_BAR"`, `"NU_X"`
(standing for μ and τ) and optionally `"NU_X_BAR"` (`"NU_X"` otherwise).

The second form takes the columns `TIME`, `L_<flavor>`, `E_<flavor>` and
`ALPHA_<flavor>`. Outside the simulated times the spectra are zero.
"""
struct PinchedModel <: SupernovaModel
    time::Vector{Float64}
    luminosity::Matrix{Float64}
    mean_energy::Matrix{Float64}
    pinch::Matrix{Float64}
end

function PinchedModel(time::AbstractVector; luminosity::AbstractDict, mean_energy::AbstractDict, pinch::AbstractDict)
    keep = _after_restarts(time)
    table(d) = permutedims(reduce(hcat, [_column(d, f)[keep] for f in ThreeFlavor]))
    return PinchedModel(collect(float.(time[keep])), table(luminosity), table(mean_energy), table(pinch))
end

# Rows of strictly increasing times. When a simulation restarts from an
# earlier time (the time column steps back), the rows it recomputes are taken
# from the restart.
function _after_restarts(time::AbstractVector)
    keep = Int[]
    earliest = Inf
    for i in reverse(eachindex(time))
        time[i] < earliest || continue
        push!(keep, i)
        earliest = time[i]
    end
    return reverse!(keep)
end

function PinchedModel(columns::AbstractDict)
    prefixed(p) = Dict(k[length(p)+1:end] => v for (k, v) in columns if startswith(k, p))
    return PinchedModel(columns["TIME"]; luminosity = prefixed("L_"), mean_energy = prefixed("E_"), pinch = prefixed("ALPHA_"))
end

# The column of flavor f in a dictionary keyed by flavors or by names, where
# NU_X stands for μ and τ and NU_X_BAR defaults to NU_X
function _column(d::AbstractDict, f::Flavor)
    for key in (f, string(f), Symbol(f))
        haskey(d, key) && return collect(float.(d[key]))
    end
    if is_muon(f) || is_tauon(f)
        is_antineutrino(f) && haskey(d, "NU_X_BAR") && return collect(float.(d["NU_X_BAR"]))
        haskey(d, "NU_X") && return collect(float.(d["NU_X"]))
    end
    throw(SupernovaNeutrinosError("no data for $(string(f)) in the model columns $(join(keys(d), ", "))"))
end

function spectra(m::PinchedModel, t::AbstractVector, E::AbstractVector)
    S = zeros(length(ThreeFlavor), length(t), length(E))
    for f in 1:length(ThreeFlavor), (i, ti) in enumerate(t)
        (ti < first(m.time) || ti > last(m.time)) && continue
        L = _interp(m.time, view(m.luminosity, f, :), ti)
        Ea = _interp(m.time, view(m.mean_energy, f, :), ti)
        α = _interp(m.time, view(m.pinch, f, :), ti)
        for (j, Ej) in enumerate(E)
            S[f, i, j] = pinched_spectrum(L, Ea, α, Ej)
        end
    end
    return S
end

"Pinch parameter from the mean energy and the mean squared energy ⟨E²⟩."
alpha_from_e2(E, E2) = @. (2E^2 - E2) / (E2 - E^2)
"Pinch parameter from the mean energy and the rms energy √⟨E²⟩."
alpha_from_rms(E, rms) = alpha_from_e2(E, rms .^ 2)

function Base.show(io::IO, m::SupernovaModel)
    print(io, nameof(typeof(m)))
    md = metadata(m)
    isempty(md) || print(io, "(", join(["$k = $(_show_value(v))" for (k, v) in md], ", "), ")")
end

function Base.show(io::IO, ::MIME"text/plain", m::SupernovaModel)
    print(io, nameof(typeof(m)), " model")
    file = _filename(m)
    isempty(file) || print(io, ": ", file)
    t = times(m)
    isempty(t) || @printf(io, "\n  %d time points, %g … %g s", length(t), first(t), last(t))
    md = metadata(m)
    width = isempty(md) ? 0 : max(18, maximum(length ∘ string, keys(md)) + 2)   # at least two spaces after the longest label
    for (k, v) in md
        print(io, "\n  ", rpad(k, width), _show_value(v))
    end
end
_show_value(v) = v isa AbstractFloat ? @sprintf("%g", v) : isnothing(v) ? "none" : string(v)
_filename(m::SupernovaModel) = hasfield(typeof(m), :filename) ? getfield(m, :filename) : ""
