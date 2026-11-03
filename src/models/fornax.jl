# Models of the Fornax code (Burrows et al.): luminosity spectra dL/dE in 12
# energy groups at each time snapshot, in HDF5 files. Need `using HDF5`.

"""
    BinnedSpectraModel(time, energy, dLdE)

A model given as luminosity spectra in energy groups: `energy[f, k, i]` is the
center (MeV) of group `k` of flavor `f` at time `time[i]` (s), and
`dLdE[f, k, i]` the luminosity spectrum in erg/(s MeV). The number spectrum is
dL/dE / E, interpolated linearly in log E between the group centers (and to
zero one group beyond the last and at zero energy), and linearly in time
between snapshots. The Fornax models are of this kind.
"""
struct BinnedSpectraModel <: SupernovaModel
    time::Vector{Float64}
    energy::Array{Float64, 3}
    dLdE::Array{Float64, 3}
end

# Number spectrum (1/(MeV s)) at E of one snapshot of one flavor
function _snapshot_spectrum(Ek::AbstractVector, dLk::AbstractVector, E::Real)
    E > 0 || return 0.0
    logEk = log10.(Ek)
    # pad with zero flux at ε and one group beyond the last
    x = [log10(eps()); logEk; logEk[end] + (logEk[end] - logEk[end-1])]
    y = [0.0; dLk; 0.0]
    return _interp_clamped(x, y, log10(E)) / MEV_IN_ERG / E
end

function spectra(m::BinnedSpectraModel, t::AbstractVector, E::AbstractVector)
    S = zeros(length(ThreeFlavor), length(t), length(E))
    nt = length(m.time)
    for (i, ti) in enumerate(t)
        (ti < first(m.time) || ti > last(m.time)) && continue
        k = min(searchsortedlast(m.time, ti), nt - 1)
        w = nt == 1 ? 0.0 : (ti - m.time[k]) / (m.time[k+1] - m.time[k])
        for f in 1:length(ThreeFlavor), (j, Ej) in enumerate(E)
            s0 = _snapshot_spectrum(view(m.energy, f, :, k), view(m.dLdE, f, :, k), Ej)
            s1 = w == 0 ? s0 : _snapshot_spectrum(view(m.energy, f, :, k + 1), view(m.dLdE, f, :, k + 1), Ej)
            S[f, i, j] = (1 - w) * s0 + w * s1
        end
    end
    return S
end

# HDF5 key of each flavor; νx stands for the four heavy flavors together
_fornax_key(f::Flavor) = is_electron(f) ? (is_neutrino(f) ? "nu0" : "nu1") : "nu2"
_fornax_factor(f::Flavor) = is_electron(f) ? 1.0 : 0.25

# Fornax 2021, 2022, 2024: groups "nu0", "nu1", "nu2" with "egroup" (MeV) and
# "g0"…"g11", dL/dE of each group in 1e50 erg/(s MeV)
function _read_fornax(h)
    time = vec(h["nu0"]["@attrs"]["time"])
    ng = size(h["nu0"]["egroup"], 1)
    energy = zeros(length(ThreeFlavor), ng, length(time))
    dLdE = zeros(length(ThreeFlavor), ng, length(time))
    for f in ThreeFlavor
        g = h[_fornax_key(f)]
        energy[f.index, :, :] = g["egroup"]
        for k in 1:ng
            dLdE[f.index, k, :] = g["g$(k-1)"] .* (1e50 * _fornax_factor(f))
        end
    end
    return BinnedSpectraModel(time, energy, dLdE)
end

"Model from an HDF5 file of Fornax 2021, 2022 or 2024. Needs `using HDF5`."
read_fornax(path::AbstractString) = _read_fornax(_hdf5(path, "This Fornax model"))

"""
    Fornax_2021(; progenitor_mass)

3D models of Burrows et al. (Fornax), https://www.astro.princeton.edu/~burrows/nu-emissions.3d/:
`progenitor_mass` 12–23, 25, 26, 26.99 M☉. Needs `using HDF5`.
"""
@registry_model Fornax_2021
parameters(::Type{Fornax_2021}) = (progenitor_mass = Parameter([12.0:23.0; 25.0; 26.0; 26.99]; unit = MSUN),)
data_files(::Type{Fornax_2021}, p::NamedTuple) =
    [isinteger(p.progenitor_mass) ? @sprintf("lum_spec_%2dM_r10000_dat.h5", round(Int, p.progenitor_mass)) :
                                    @sprintf("lum_spec_%.2fM_r10000_dat.h5", p.progenitor_mass)]
function Fornax_2021(; kwargs...)
    p = validate(Fornax_2021; kwargs...)
    file = only(data_files(Fornax_2021, p))
    return _registry(Fornax_2021, p, _read_fornax(_hdf5(model_file(Fornax_2021, file), "Fornax_2021")), file)
end

const _FORNAX_2022_PROGENITORS = [
    "9.0", "9.25", "9.5", "9.75", "10.0", "10.25", "10.5", "10.75", "11.0", "11.25", "11.5", "11.75",
    "12.00.bh", "12.03.bh", "12.07.bh", "12.1.bh", "12.13", "12.15", "12.18.bh", "12.20.bh", "12.25", "12.33.bh",
    "12.40.bh", "12.45.bh", "12.50.bh", "12.54.bh", "12.60.bh", "12.63", "12.70", "12.72.bh", "12.75", "12.80.bh",
    "12.85.bh", "12.90.bh", "12.93", "12.97.bh", "13.00.bh", "13.05.bh", "13.11", "13.25.bh", "13.27.bh", "13.32.bh",
    "13.40.bh", "13.45", "13.50.bh", "13.60.bh", "13.75", "13.82.bh", "13.90.bh", "13.96", "14.01", "14.13.bh",
    "14.25.bh", "14.40.bh", "14.41.bh", "14.43", "14.44.bh", "14.70.bh", "14.87.bh", "15.00.bh", "15.01", "15.04.bh",
    "15.05", "15.38.bh", "16.43", "16.65", "16.99", "17.00", "17.07", "17.10", "17.40", "17.48", "17.50", "17.51",
    "17.83", "18.04", "18.05", "18.09", "18.10", "18.50", "19.02", "19.56", "19.83", "19.99", "20.08", "20.09",
    "20.18", "20.37", "21.00", "21.68", "22.00", "22.30", "22.82", "23.00", "23.04", "23.43", "24.00", "25.00",
    "26.00", "26.99"]
_fornax_2022_mass(p) = parse(Float64, replace(p, ".bh" => ""))

"""
    Fornax_2022(; progenitor_mass)

2D models of 100 progenitors of Wang, Vartanyan, Burrows & Coleman, MNRAS
517, 543 (2022), https://www.astro.princeton.edu/~burrows/nu-emissions.2d.large/:
`progenitor_mass` 9–26.99 M☉ (see `parameters(Fornax_2022)`). Needs `using HDF5`.
"""
@registry_model Fornax_2022
parameters(::Type{Fornax_2022}) = (progenitor_mass = Parameter(_fornax_2022_mass.(_FORNAX_2022_PROGENITORS); unit = MSUN),)
_fornax_2022_progenitor(m) = _FORNAX_2022_PROGENITORS[findfirst(p -> _same(_fornax_2022_mass(p), m), _FORNAX_2022_PROGENITORS)]
data_files(::Type{Fornax_2022}, p::NamedTuple) = ["lum_spec_$(_fornax_2022_progenitor(p.progenitor_mass))_dat.h5"]

const _FORNAX_2024_PROGENITORS = ["u8.1", "9b", "9.25", "9.5", "z9.6", "11", "12.25", "14", "15.01", "16.5", "16",
    "17", "18", "18.5", "19", "19.56", "20", "21.68", "23", "24", "25", "40", "60", "100"]
const _FORNAX_2024_BLACK_HOLES = ("12.25", "14", "19.56", "40")
_fornax_2024_mass(p) = parse(Float64, replace(p, r"[A-Za-z]" => ""))

"""
    Fornax_2024(; progenitor_mass)

3D models of 24 progenitors with the Fornax code,
https://www.astro.princeton.edu/~burrows/nu-emissions.3d.update/: `progenitor_mass`
8.1–100 M☉ (see `parameters(Fornax_2024)`). Needs `using HDF5`.
"""
@registry_model Fornax_2024
parameters(::Type{Fornax_2024}) = (progenitor_mass = Parameter(_fornax_2024_mass.(_FORNAX_2024_PROGENITORS); unit = MSUN),)
_fornax_2024_progenitor(m) = _FORNAX_2024_PROGENITORS[findfirst(p -> _same(_fornax_2024_mass(p), m), _FORNAX_2024_PROGENITORS)]
data_files(::Type{Fornax_2024}, p::NamedTuple) = ["lum_spec_$(_fornax_2024_progenitor(p.progenitor_mass))_dat.h5"]

for (M, progenitor, black_hole) in ((:Fornax_2022, :_fornax_2022_progenitor, p -> endswith(p, ".bh")),
                                    (:Fornax_2024, :_fornax_2024_progenitor, p -> p in _FORNAX_2024_BLACK_HOLES))
    @eval function $M(; kwargs...)
        p = validate($M; kwargs...)
        file = only(data_files($M, p))
        h = _hdf5(model_file($M, file), $(string(M)))
        return _registry($M, p, _read_fornax(h), file;
                         var"Black hole" = $black_hole($progenitor(p.progenitor_mass)),
                         var"PNS mass" = only(h["@attrs"]["Mpns"]))
    end
end

# ----------------------------------------------------------------------------
# Fornax 2019: 3D models with the angular dependence of the emission as real
# spherical harmonics up to l = 2

# Associated Legendre functions P_l^m(x) with the Condon–Shortley phase, m ≥ 0, l ≤ 2
function _legendre(l::Int, m::Int, x::Real)
    s = sqrt(max(0.0, 1 - x^2))
    (l, m) == (0, 0) && return 1.0
    (l, m) == (1, 0) && return x
    (l, m) == (1, 1) && return -s
    (l, m) == (2, 0) && return (3x^2 - 1) / 2
    (l, m) == (2, 1) && return -3x * s
    (l, m) == (2, 2) && return 3 * (1 - x^2)
    throw(ArgumentError("P_$l^$m is not implemented"))
end

"Orthonormal real spherical harmonic Y_lm at zenith angle θ and azimuth φ (radians)."
function real_spherical_harmonic(l::Int, m::Int, θ::Real, ϕ::Real)
    m == 0 && return sqrt((2l + 1) / 4π) * _legendre(l, 0, cos(θ))
    norm = sqrt((2l + 1) / 2π * factorial(l - abs(m)) / factorial(l + abs(m)))
    return m > 0 ? norm * _legendre(l, m, cos(θ)) * cos(m * ϕ) : norm * _legendre(l, -m, cos(θ)) * sin(-m * ϕ)
end

"""
    read_fornax_2019(path; theta = 0, phi = 0) -> BinnedSpectraModel

Model from an HDF5 file of Fornax 2019, seen from the direction (`theta`,
`phi`) in radians. Needs `using HDF5`.
"""
read_fornax_2019(path::AbstractString; theta::Real = 0.0, phi::Real = 0.0) =
    _read_fornax_2019(_hdf5(path, "Fornax_2019"), theta, phi)

function _read_fornax_2019(h, θ, ϕ)
    time = vec(h["nu0"]["g0"]["@attrs"]["time"])
    ng = size(h["nu0"]["egroup"], 1)
    energy = zeros(length(ThreeFlavor), ng, length(time))
    dLdE = zeros(length(ThreeFlavor), ng, length(time))
    Y = Dict((l, m) => real_spherical_harmonic(l, m, θ, ϕ) for l in 0:2 for m in -l:l)
    for f in ThreeFlavor
        g = h[_fornax_key(f)]
        energy[f.index, :, :] = g["egroup"]
        for k in 1:ng
            moments = g["g$(k-1)"]
            L = sum(moments["l=$l m=$m"] .* Y[(l, m)] for l in 0:2 for m in -l:l)
            dLdE[f.index, k, :] = abs.(L) .* (1e50 * _fornax_factor(f))
        end
    end
    return BinnedSpectraModel(time, energy, dLdE)
end

"""
    Fornax_2019(; progenitor_mass, theta = 0, phi = 0)

3D models of Vartanyan, Burrows, Radice, Skinner & Dolence, MNRAS 482, 351
(2019), https://www.astro.princeton.edu/~burrows/nu-emissions.3d/:
`progenitor_mass` 9, 10, 12, 13, 14, 15, 19, 25 or 60 M☉, seen from the
direction (`theta`, `phi`) in radians. Needs `using HDF5`.

The emission is evaluated at the exact direction from the spherical-harmonic
moments (SNEWPY uses the nearest pixel of a HEALPix grid with 192 pixels).
"""
@registry_model Fornax_2019
parameters(::Type{Fornax_2019}) = (progenitor_mass = Parameter([9.0, 10.0, 12.0, 13.0, 14.0, 15.0, 19.0, 25.0, 60.0]; unit = MSUN),)
data_files(::Type{Fornax_2019}, p::NamedTuple) = [@sprintf("lum_spec_%dM.h5", round(Int, p.progenitor_mass))]
function Fornax_2019(; theta::Real = 0.0, phi::Real = 0.0, kwargs...)
    p = validate(Fornax_2019; kwargs...)
    file = only(data_files(Fornax_2019, p))
    base = _read_fornax_2019(_hdf5(model_file(Fornax_2019, file), "Fornax_2019"), theta, phi)
    return _registry(Fornax_2019, p, base, file; theta = float(theta), phi = float(phi))
end
