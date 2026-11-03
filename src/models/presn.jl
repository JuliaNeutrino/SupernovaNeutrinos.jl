# Pre-supernova models: the neutrino emission of the star in the hours and
# days before core collapse, from https://github.com/SNEWS2/snewpy-models-presn.
# Times are negative, before collapse.

const _DT, _DV = 1e-3, 1e-10

# Interpolation in log(|t| + dt) of log(v + dv) ("dilog"), zero outside the
# tabulated times; |t| in the native time unit of the model
function _dilog_interp(tgrid::AbstractVector, v::AbstractVector, t::Real)
    (t < first(tgrid) || t > last(tgrid)) && return 0.0
    lv = _interp(log.(tgrid .+ _DT), log.(v .+ _DV), log(t + _DT))
    return max(exp(lv) - _DV, 0.0)
end

"""
    TabulatedModel(time, energy, array; time_unit = 1.0)

Model of tabulated spectra d²N/dE dt (1/(MeV s)): `array[f, i, k]` for the
flavors of [`ThreeFlavor`](@ref) at times `time[i]` (s, negative before
collapse) and energies `energy[k]` (MeV). The spectra are interpolated in
log(|t| + 10⁻³ `time_unit`) and log(value) between times, linearly in
energy, and are zero outside the table. The pre-supernova models are of this
kind.
"""
struct TabulatedModel <: SupernovaModel
    time::Vector{Float64}
    energy::Vector{Float64}
    array::Array{Float64, 3}
    time_unit::Float64
    function TabulatedModel(time::AbstractVector, energy::AbstractVector, array::AbstractArray{<:Real, 3}; time_unit::Real = 1.0)
        p = sortperm(time)
        q = sortperm(energy)
        return new(collect(float.(time[p])) .+ 0.0, collect(float.(energy[q])), Array{Float64}(array[:, p, q]), float(time_unit))
    end
end

function spectra(m::TabulatedModel, t::AbstractVector, E::AbstractVector)
    # interpolate in |t| (native units), sorted increasingly
    order = sortperm(abs.(m.time))
    tgrid = abs.(m.time[order]) ./ m.time_unit
    S = zeros(length(ThreeFlavor), length(t), length(E))
    for f in 1:length(ThreeFlavor), (i, ti) in enumerate(t)
        (ti < first(m.time) || ti > last(m.time)) && continue
        tq = abs(ti) / m.time_unit
        at_t = [_dilog_interp(tgrid, view(m.array, f, order, k), tq) for k in eachindex(m.energy)]
        for (j, Ej) in enumerate(E)
            S[f, i, j] = _interp(m.energy, at_t, Ej)
        end
    end
    return S
end

# ----------------------------------------------------------------------------

"""
    OdrzywolekModel(time, a, alpha, b)

Pre-supernova spectra of Odrzywolek & Heger: dN/dE dt = a E^α exp(-b E) for
ν̄e, with the parameters interpolated in time like [`TabulatedModel`](@ref);
νe has the same spectrum and each heavy flavor 0.19 of it
(arXiv:astro-ph/0311012).
"""
struct OdrzywolekModel <: SupernovaModel
    time::Vector{Float64}
    a::Vector{Float64}
    alpha::Vector{Float64}
    b::Vector{Float64}
end

function spectra(m::OdrzywolekModel, t::AbstractVector, E::AbstractVector)
    order = sortperm(abs.(m.time))
    tgrid = abs.(m.time[order])
    S = zeros(length(ThreeFlavor), length(t), length(E))
    for (i, ti) in enumerate(t)
        (ti < minimum(m.time) || ti > maximum(m.time)) && continue
        a, α, b = (_dilog_interp(tgrid, v[order], abs(ti)) for v in (m.a, m.alpha, m.b))
        for (j, Ej) in enumerate(E), f in ThreeFlavor
            S[f.index, i, j] = a * Ej^α * exp(-b * Ej) * (is_electron(f) ? 1.0 : 0.19)
        end
    end
    return S
end

"""
    Odrzywolek_2010(; progenitor_mass)

Pre-supernova model of Odrzywolek & Heger, Acta Phys. Polon. B 41, 1611
(2010): `progenitor_mass` 15 or 25 M☉.
"""
@registry_model Odrzywolek_2010
parameters(::Type{Odrzywolek_2010}) = (progenitor_mass = Parameter([15.0, 25.0]; unit = MSUN),)
repository(::Type{Odrzywolek_2010}) = PRESN_REPOSITORY
data_files(::Type{Odrzywolek_2010}, p::NamedTuple) = [@sprintf("s%.0f_nuebar_data.txt", p.progenitor_mass)]
function Odrzywolek_2010(; kwargs...)
    p = validate(Odrzywolek_2010; kwargs...)
    file = only(data_files(Odrzywolek_2010, p))
    # columns: step, time before collapse (s), …, a, α, b
    A = read_table(model_file(Odrzywolek_2010, file); skip = 1)
    return _registry(Odrzywolek_2010, p, OdrzywolekModel(0.0 .- A[:, 2], A[:, 7], A[:, 8], A[:, 9]), file)
end

"""
    Patton_2017(; progenitor_mass)

Pre-supernova model of Patton et al., ApJ 851, 6 (2017): `progenitor_mass`
15 or 30 M☉.
"""
@registry_model Patton_2017
parameters(::Type{Patton_2017}) = (progenitor_mass = Parameter([15.0, 30.0]; unit = MSUN),)
repository(::Type{Patton_2017}) = PRESN_REPOSITORY
data_files(::Type{Patton_2017}, p::NamedTuple) = [@sprintf("totalLuminosity_%.0fSolarMass.dat", p.progenitor_mass)]
function Patton_2017(; kwargs...)
    p = validate(Patton_2017; kwargs...)
    file = only(data_files(Patton_2017, p))
    # columns: time before collapse (h), energy (MeV), spectra of νe, ν̄e, νx, ν̄x (1/(MeV s)), per heavy flavor
    A = read_table(model_file(Patton_2017, file))
    tvals, evals = sort(unique(A[:, 1])), sort(unique(A[:, 2]))
    it = Dict(v => i for (i, v) in enumerate(tvals))
    ie = Dict(v => i for (i, v) in enumerate(evals))
    array = zeros(length(ThreeFlavor), length(tvals), length(evals))
    columns = (NU_E = 3, NU_E_BAR = 4, NU_MU = 5, NU_MU_BAR = 6, NU_TAU = 5, NU_TAU_BAR = 6)
    for r in axes(A, 1), f in ThreeFlavor
        array[f.index, it[A[r, 1]], ie[A[r, 2]]] = A[r, columns[Symbol(f)]]
    end
    return _registry(Patton_2017, p, TabulatedModel(-tvals .* HOUR, evals, array; time_unit = HOUR), file)
end

"""
    Kato_2017(; progenitor_mass)

Pre-supernova model of Kato et al., ApJ 848, 48 (2017): `progenitor_mass` 12
or 15 M☉. Downloads one spectrum file per time step.
"""
@registry_model Kato_2017
parameters(::Type{Kato_2017}) = (progenitor_mass = Parameter([12.0, 15.0]; unit = MSUN),)
repository(::Type{Kato_2017}) = PRESN_REPOSITORY

const _KATO_FILES = (NU_E = "total_nue/spe_all", NU_E_BAR = "total_nueb/spe_all", NU_MU = "total_nux/spe_sum_mu_nu",
                     NU_MU_BAR = "total_nux/spe_sum_mu", NU_TAU = "total_nux/spe_sum_mu_nu", NU_TAU_BAR = "total_nux/spe_sum_mu")
_kato_dir(p) = @sprintf("pre_collapse/m%.0f", p.progenitor_mass)
function _kato_steps(p)
    A = read_table(model_file(Kato_2017, "$(_kato_dir(p))/total_nue/lightcurve_nue_all.dat"))
    return A[:, 1], round.(Int, A[:, 4])
end
function data_files(::Type{Kato_2017}, p::NamedTuple)
    _, steps = _kato_steps(p)
    return [["$(_kato_dir(p))/total_nue/lightcurve_nue_all.dat"];
            [@sprintf("%s/%s%05d.dat", _kato_dir(p), base, s) for base in unique(values(_KATO_FILES)) for s in steps]]
end
function Kato_2017(; kwargs...)
    p = validate(Kato_2017; kwargs...)
    tvals, steps = _kato_steps(p)
    energy = Float64[]
    array = zeros(length(ThreeFlavor), length(steps), 0)
    for f in ThreeFlavor, (i, s) in enumerate(steps)
        A = read_table(model_file(Kato_2017, @sprintf("%s/%s%05d.dat", _kato_dir(p), _KATO_FILES[Symbol(f)], s)))
        if isempty(energy)
            energy = A[:, 1]
            array = zeros(length(ThreeFlavor), length(steps), length(energy))
        end
        array[f.index, i, :] = A[:, 2]
    end
    return _registry(Kato_2017, p, TabulatedModel(-tvals, energy, array), _kato_dir(p))
end

"""
    Yoshida_2016(; progenitor_mass)

Pre-supernova model of Yoshida et al., PRD 93, 123012 (2016):
`progenitor_mass` 12, 15 or 20 M☉.
"""
@registry_model Yoshida_2016
parameters(::Type{Yoshida_2016}) = (progenitor_mass = Parameter([12.0, 15.0, 20.0]; unit = MSUN),)
repository(::Type{Yoshida_2016}) = PRESN_REPOSITORY
data_files(::Type{Yoshida_2016}, p::NamedTuple) = [@sprintf("t_spc_m%.0f_1.2.txt", p.progenitor_mass)]
function Yoshida_2016(; kwargs...)
    p = validate(Yoshida_2016; kwargs...)
    file = only(data_files(Yoshida_2016, p))
    # blocks of one header line ("step time …") and, for νe, ν̄e, νx, ν̄x, 100
    # lines of the spectrum at 1500 energies: 0.01–10 MeV (step 0.01) and 10.02–20 MeV (step 0.02)
    energy = [range(0, 10; length = 1001)[2:end]; range(10, 20; length = 501)[2:end]]
    tvals = Float64[]
    blocks = Matrix{Float64}[]
    open(model_file(Yoshida_2016, file)) do io
        while !eof(io)
            header = strip(readline(io))
            isempty(header) && continue
            push!(tvals, parse(Float64, split(header)[2]))
            push!(blocks, reduce(hcat, [reduce(vcat, (_parse_floats(readline(io)) for _ in 1:100)) for _ in 1:4]))
        end
    end
    array = zeros(length(ThreeFlavor), length(tvals), length(energy))
    # νx and ν̄x are the sums over μ and τ
    source = (NU_E = (1, 1.0), NU_E_BAR = (2, 1.0), NU_MU = (3, 0.5), NU_MU_BAR = (4, 0.5), NU_TAU = (3, 0.5), NU_TAU_BAR = (4, 0.5))
    for (i, B) in enumerate(blocks), f in ThreeFlavor
        col, factor = source[Symbol(f)]
        array[f.index, i, :] = B[:, col] .* factor
    end
    return _registry(Yoshida_2016, p, TabulatedModel(-tvals, energy, array), file)
end

"""
    Myers_2026(; progenitor_mass, eta, bounds)

Pre-supernova model of Myers et al. (2026) [arXiv:2604.22605]:
`progenitor_mass` 12, 15, 18 or 20 M☉; `eta` 0.2, 0.4, 0.8 or 1.0; `bounds`
`"all"` or `"core"`. The μ+τ spectra of the pair processes are split equally
between μ and τ.
"""
@registry_model Myers_2026
parameters(::Type{Myers_2026}) = (progenitor_mass = Parameter([12.0, 15.0, 18.0, 20.0]; unit = MSUN),
    eta = Parameter([0.2, 0.4, 0.8, 1.0]), bounds = Parameter(["all", "core"]))
repository(::Type{Myers_2026}) = PRESN_REPOSITORY_MAIN
data_files(::Type{Myers_2026}, p::NamedTuple) = [@sprintf("%.0fM_%.1f_%s.zip", p.progenitor_mass, p.eta, p.bounds)]

function _myers_snapshot(text::AbstractString)
    m = match(r"t_cc\s*=\s*([-+0-9.eE]+)", text)
    isnothing(m) && throw(SupernovaNeutrinosError("no t_cc in a Myers_2026 file"))
    A = read_table(IOBuffer(text))
    A[1, 1] == 0 && (A = A[2:end, :])            # drop the E = 0 row
    return parse(Float64, m[1]), A[:, 1], A[:, 2:end]
end

function Myers_2026(; kwargs...)
    p = validate(Myers_2026; kwargs...)
    file = only(data_files(Myers_2026, p))
    profile(name) = parse(Int, match(r"LuminosityProfile(\d+)\.dat$", name)[1])
    members = read_zip_members(n -> endswith(n, ".dat") && (occursin("/betaLuminosity/", n) || occursin("/pairLuminosity/", n)),
                               model_file(Myers_2026, file))
    beta = Dict(profile(n) => _myers_snapshot(s) for (n, s) in members if occursin("/betaLuminosity/", n))
    pair = Dict(profile(n) => _myers_snapshot(s) for (n, s) in members if occursin("/pairLuminosity/", n))
    # the last 1000 hours before collapse, in the profiles of both processes
    profiles = sort([k for k in intersect(keys(beta), keys(pair)) if beta[k][1] >= -1000]; by = k -> beta[k][1])
    tvals = [beta[k][1] for k in profiles]
    energy = beta[first(profiles)][2]
    array = zeros(length(ThreeFlavor), length(profiles), length(energy))
    for (i, k) in enumerate(profiles)
        b, q = beta[k][3], pair[k][3]
        array[ThreeFlavor.NU_E.index, i, :] = b[:, 1] .+ q[:, 1]
        array[ThreeFlavor.NU_E_BAR.index, i, :] = b[:, 2] .+ q[:, 2]
        for f in (:NU_MU, :NU_TAU)
            array[ThreeFlavor[f].index, i, :] = q[:, 3] ./ 2
        end
        for f in (:NU_MU_BAR, :NU_TAU_BAR)
            array[ThreeFlavor[f].index, i, :] = q[:, 4] ./ 2
        end
    end
    return _registry(Myers_2026, p, TabulatedModel(tvals .* HOUR, energy, array; time_unit = HOUR), file)
end
