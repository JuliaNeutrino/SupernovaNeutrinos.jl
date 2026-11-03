# Core-collapse supernova models of the SNEWS model repository
# (https://github.com/SNEWS2/snewpy-models-ccsn). Each model is created from
# its physics parameters and downloads its data files the first time:
#
#     Nakazato_2013(progenitor_mass = 13, revival_time = 100, metallicity = 0.02, eos = "shen")

const MSUN = "M☉"

# Hooks of the FITSIO and HDF5 extensions, which read the binary model files
function read_fits_table end
function read_hdf5 end

function _need_extension(name::Symbol, pkg::String, model)
    isnothing(Base.get_extension(@__MODULE__, name)) && throw(SupernovaNeutrinosError(
        "$model reads $pkg files: install and load $pkg first (`using $pkg`)"))
end

_fits_table(path, model) = (_need_extension(:SupernovaNeutrinosFITSIOExt, "FITSIO", model); read_fits_table(path))
_hdf5(path, model) = (_need_extension(:SupernovaNeutrinosHDF5Ext, "HDF5", model); read_hdf5(path))

_registry(M, p, base, file; extra...) = M(base, merge!(_metadata(M, p), OrderedDict{String, Any}(String(k) => v for (k, v) in extra)), file)

# ----------------------------------------------------------------------------
# Garching Supernova Archive format: three files per model (nue, nuebar, nux)
# with columns time (s), luminosity (1e51 erg/s), mean energy (MeV) and mean
# squared energy (MeV²)

const _GARCHING_FLAVORS = (("nue", "NU_E"), ("nuebar", "NU_E_BAR"), ("nux", "NU_X"))

"""
    read_garching(prefix; eos = "LS220") -> PinchedModel

Model in the format of the Garching Supernova Archive, from the files
`<prefix>_<eos>_nue`, `<prefix>_<eos>_nuebar` and `<prefix>_<eos>_nux`.
"""
read_garching(prefix::AbstractString; eos::AbstractString = "LS220") = _read_garching(f -> f, prefix, eos)

function _read_garching(file, prefix, eos)
    tables = Dict(sfx => read_table(file("$(prefix)_$(eos)_$(sfx)")) for (sfx, _) in _GARCHING_FLAVORS)
    t = tables["nue"][:, 1]
    cols = Dict{String, Vector{Float64}}("TIME" => t)
    for (sfx, f) in _GARCHING_FLAVORS
        A = tables[sfx]
        # the times of the other flavors are matched to those of νe; missing ones are zero
        row = Dict(A[i, 1] => i for i in axes(A, 1))
        idx = [get(row, ti, 0) for ti in t]
        col(c) = [i == 0 ? 0.0 : A[i, c] for i in idx]
        E, E2 = col(3), col(4)
        cols["L_$f"] = col(2) .* 1e51
        cols["E_$f"] = E
        cols["ALPHA_$f"] = [i == 0 ? 0.0 : a for (i, a) in zip(idx, alpha_from_e2(E, E2))]
    end
    return PinchedModel(cols)
end

"""
    Bollig_2016(; progenitor_mass, eos = "LS220")

Model of Bollig et al. (2016) [arXiv:1508.00785], from the Garching
Supernova Archive. `progenitor_mass`: 11.2 or 27 M☉.
"""
@registry_model Bollig_2016
parameters(::Type{Bollig_2016}) = (progenitor_mass = Parameter([11.2, 27.0]; unit = MSUN), eos = Parameter(["LS220"]))
_garching_prefix(::Type{Bollig_2016}, p) = @sprintf("s%.1fc", p.progenitor_mass)

"""
    Tamborra_2014(; progenitor_mass, direction, eos = "LS220")

3D model of Tamborra et al., PRD 90, 045032 (2014) [arXiv:1406.0006], from the
Garching Supernova Archive. `progenitor_mass` 11.2 and 20 M☉ (direction 1)
or 27 M☉ (directions 1, 2, 3).
"""
@registry_model Tamborra_2014
parameters(::Type{Tamborra_2014}) = (progenitor_mass = Parameter([11.2, 20.0, 27.0]; unit = MSUN),
    direction = Parameter([1, 2, 3]), eos = Parameter(["LS220"]))
valid_combination(::Type{Tamborra_2014}, p::NamedTuple) = p.direction == 1 || p.progenitor_mass == 27
_garching_prefix(::Type{Tamborra_2014}, p) = @sprintf("s%.1fc_3D_dir%d", p.progenitor_mass, p.direction)

"""
    Walk_2018(; rotation, direction, progenitor_mass = 15, eos = "LS220")

SASI-dominated 3D model of Walk et al., PRD 98, 123001 (2018)
[arXiv:1807.02366], from the Garching Supernova Archive. `rotation`:
`"fast"`, `"slow"` or `"non"`; `direction`: 1, 2 or 3.
"""
@registry_model Walk_2018
parameters(::Type{Walk_2018}) = (progenitor_mass = Parameter([15.0]; unit = MSUN),
    rotation = Parameter(["fast", "slow", "non"]), direction = Parameter([1, 2, 3]), eos = Parameter(["LS220"]))
_garching_prefix(::Type{Walk_2018}, p) = @sprintf("s%.1fc_3D_%srot_dir%d", p.progenitor_mass, p.rotation, p.direction)

"""
    Walk_2019(; progenitor_mass, direction, eos = "LS220")

3D black-hole forming model of Walk et al., PRD 101, 123013 (2019)
[arXiv:1910.12971], from the Garching Supernova Archive. `progenitor_mass`
40 M☉ (directions 1, 2, 3) or 75 M☉ (directions 1, 2).
"""
@registry_model Walk_2019
parameters(::Type{Walk_2019}) = (progenitor_mass = Parameter([40.0, 75.0]; unit = MSUN),
    direction = Parameter([1, 2, 3]), eos = Parameter(["LS220"]))
valid_combination(::Type{Walk_2019}, p::NamedTuple) = p.progenitor_mass == 40 || p.direction in (1, 2)
_garching_prefix(::Type{Walk_2019}, p) = @sprintf("s%.1fc_3DBH_dir%d", p.progenitor_mass, p.direction)

const GarchingModel = Union{Bollig_2016, Tamborra_2014, Walk_2018, Walk_2019}

data_files(M::Type{<:GarchingModel}, p::NamedTuple) =
    ["$(_garching_prefix(M, p))_$(p.eos)_$(sfx)" for (sfx, _) in _GARCHING_FLAVORS]

function (M::Type{<:GarchingModel})(; kwargs...)
    p = validate(M; kwargs...)
    prefix = _garching_prefix(M, p)
    return _registry(M, p, _read_garching(f -> model_file(M, f), prefix, p.eos), prefix)
end

# ----------------------------------------------------------------------------

"""
    Nakazato_2013(; progenitor_mass, revival_time, metallicity, eos)

Models of Nakazato et al., ApJS 205, 2 (2013), ApJ 804, 75 (2015), PASJ 73,
639 (2021): `progenitor_mass` 13, 20, 30, 50 M☉; shock `revival_time` 0
(black hole), 100, 200, 300 ms; `metallicity` 0.02 or 0.004; `eos` `"LS220"`,
`"shen"` or `"togashi"`. Needs `using FITSIO`.
"""
@registry_model Nakazato_2013
parameters(::Type{Nakazato_2013}) = (progenitor_mass = Parameter([13.0, 20.0, 30.0, 50.0]; unit = MSUN),
    revival_time = Parameter([0.0, 100.0, 200.0, 300.0]; unit = "ms"), metallicity = Parameter([0.02, 0.004]),
    eos = Parameter(["LS220", "shen", "togashi"]))
function valid_combination(::Type{Nakazato_2013}, p::NamedTuple)
    bh = p.progenitor_mass == 30 && p.metallicity == 0.004
    return p.revival_time == 0 ? bh : (p.eos == "shen" && !bh)
end
function data_files(::Type{Nakazato_2013}, p::NamedTuple)
    m = @sprintf("%3.1f", p.progenitor_mass)
    return [p.revival_time == 0 ? "nakazato-$(p.eos)-BH-z$(p.metallicity)-s$m.fits" :
            "nakazato-$(p.eos)-z$(p.metallicity)-t_rev$(round(Int, p.revival_time))ms-s$m.fits"]
end

"""
    Sukhbold_2015(; progenitor_mass, eos)

Models of Sukhbold et al., ApJ 821, 38 (2016): `progenitor_mass` 9.6 or
27 M☉, `eos` `"LS220"` or `"SFHo"`. Needs `using FITSIO`.
"""
@registry_model Sukhbold_2015
parameters(::Type{Sukhbold_2015}) = (progenitor_mass = Parameter([27.0, 9.6]; unit = MSUN), eos = Parameter(["LS220", "SFHo"]))
data_files(::Type{Sukhbold_2015}, p::NamedTuple) =
    [@sprintf("sukhbold-%s-%s%3.1f.fits", p.eos, p.progenitor_mass == 9.6 ? "z" : "s", p.progenitor_mass)]

"Model from a FITS table with the columns of [`PinchedModel`](@ref) (Nakazato, Sukhbold)."
read_pinched_fits(path::AbstractString) = PinchedModel(_fits_table(path, "This model"))

for M in (:Nakazato_2013, :Sukhbold_2015)
    @eval function $M(; kwargs...)
        p = validate($M; kwargs...)
        file = only(data_files($M, p))
        return _registry($M, p, PinchedModel(_fits_table(model_file($M, file), $(string(M)))), file)
    end
end

# ----------------------------------------------------------------------------

"""
    OConnor_2013(; progenitor_mass, eos)

Black-hole formation models of O'Connor & Ott (2013) [arXiv:1207.1100]:
`progenitor_mass` 12–33, 35–60 (step 5), 70, 80, 100, 120 M☉; `eos`
`"HShen"` or `"LS220"`.
"""
@registry_model OConnor_2013
parameters(::Type{OConnor_2013}) = (
    progenitor_mass = Parameter(Float64[12:33; 35:5:60; 70; 80; 100; 120]; unit = MSUN), eos = Parameter(["HShen", "LS220"]))
data_files(::Type{OConnor_2013}, p::NamedTuple) = ["$(p.eos)_timeseries.tar.gz"]

function OConnor_2013(; kwargs...)
    p = validate(OConnor_2013; kwargs...)
    file = only(data_files(OConnor_2013, p))
    member = "s$(round(Int, p.progenitor_mass))_$(p.eos)_timeseries.dat"
    A = read_table(IOBuffer(read_tar_member(model_file(OConnor_2013, file), member)))
    # the νx luminosity of these files is already per flavor
    return _registry(OConnor_2013, p, PinchedModel(_pinched_from_rms(A)), member)
end

# O'Connor (2015) format: a "#tbounce=…" header line, then columns TIME, L, E,
# RMS energy for νe, ν̄e, νx with L in erg/s and the νx luminosity summed
# over the four heavy flavors
function _read_oconnor_2015(path)
    tbounce = open(path) do io
        line = readline(io)
        m = match(r"tbounce\s*=\s*([-+0-9.eE]+)", line)
        isnothing(m) ? 0.0 : parse(Float64, m[1])
    end
    A = read_table(path)
    A[:, 1] .-= tbounce
    A[:, 4] ./= 4
    # negative luminosities are numerical noise
    for c in 2:4
        A[A[:, c] .< 0, c] .= 1
    end
    return PinchedModel(_pinched_from_rms(A))
end

"""
    OConnor_2015(; progenitor_mass = 40)

Black-hole formation model of O'Connor (2015) [arXiv:1411.7058].
"""
@registry_model OConnor_2015
parameters(::Type{OConnor_2015}) = (progenitor_mass = Parameter([40.0]; unit = MSUN),)
data_files(::Type{OConnor_2015}, p::NamedTuple) = ["M1_neutrinos.dat"]
function OConnor_2015(; kwargs...)
    p = validate(OConnor_2015; kwargs...)
    return _registry(OConnor_2015, p, _read_oconnor_2015(model_file(OConnor_2015, "M1_neutrinos.dat")), "M1_neutrinos.dat"; EOS = "LS220")
end

"""
    Zha_2021(; progenitor_mass)

Hadron–quark phase transition models of Zha et al. (2021) [arXiv:2103.02268]:
`progenitor_mass` 16–26, 19.89, 22.39, 30, 33 M☉.
"""
@registry_model Zha_2021
parameters(::Type{Zha_2021}) = (progenitor_mass = Parameter([16.0:26.0; 19.89; 22.39; 30.0; 33.0]; unit = MSUN),)
data_files(::Type{Zha_2021}, p::NamedTuple) = [@sprintf("s%g.dat", p.progenitor_mass)]
function Zha_2021(; kwargs...)
    p = validate(Zha_2021; kwargs...)
    file = only(data_files(Zha_2021, p))
    return _registry(Zha_2021, p, _read_oconnor_2015(model_file(Zha_2021, file)), file; EOS = "STOS_B145")
end

# ----------------------------------------------------------------------------

const _WARREN_MASSES = [range(9.0, 12.75; length = 16); range(13.0, 30.0; length = 171); 31.0:33.0; 35.0:5.0:55.0;
                        60.0:10.0:80.0; 100.0; 120.0]

"""
    Warren_2020(; progenitor_mass, turbmixing_param)

Models of Warren et al., ApJ 898, 139 (2020): `progenitor_mass` 9–13 (step
0.25), 13–30 (step 0.1), 31–33, 35–55 (step 5), 60, 70, 80, 100, 120 M☉;
`turbmixing_param` 1.23, 1.25 or 1.27. Needs `using HDF5`.
"""
@registry_model Warren_2020
parameters(::Type{Warren_2020}) = (progenitor_mass = Parameter(_WARREN_MASSES; unit = MSUN),
    turbmixing_param = Parameter([1.23, 1.25, 1.27]; label = "Turb. mixing param.", description = "Turbulent mixing parameter α_Λ"))
function data_files(::Type{Warren_2020}, p::NamedTuple)
    m, a = p.progenitor_mass, @sprintf("%3.2f", p.turbmixing_param)
    mass = isinteger(m) && m <= 30 ? @sprintf("%.1f", m) : @sprintf("%g", m)
    return ["stir_a$a/stir_multimessenger_a$(a)_m$mass.h5"]
end

function Warren_2020(; kwargs...)
    p = validate(Warren_2020; kwargs...)
    file = only(data_files(Warren_2020, p))
    return _registry(Warren_2020, p, read_warren(model_file(Warren_2020, file)), file; EOS = "SFHo")
end

"Model from an HDF5 file of Warren et al. (2020). Needs `using HDF5`."
function read_warren(path::AbstractString)
    h = _hdf5(path, "Warren_2020")
    # HDF5 arrays arrive transposed: (column, row)
    shock = h["sim_data"]["shock_radius"]
    i = findfirst(>(1e-5), shock[2, :])
    bounce = isnothing(i) ? 0.0 : shock[1, i]
    cols = Dict{String, Vector{Float64}}("TIME" => h["nue_data"]["lum"][1, :] .- bounce)
    for (key, f) in (("nue_data", "NU_E"), ("nuae_data", "NU_E_BAR"), ("nux_data", "NU_X"))
        E, rms = h[key]["avg_energy"][2, :], h[key]["rms_energy"][2, :]
        cols["L_$f"] = h[key]["lum"][2, :] .* 1e51
        cols["E_$f"] = E
        cols["ALPHA_$f"] = alpha_from_rms(E, rms)
    end
    return PinchedModel(cols)
end

# ----------------------------------------------------------------------------

"""
    Kuroda_2020(; rotational_velocity, magnetic_field_exponent, progenitor_mass = 20)

Magnetorotational models of Kuroda et al. (2020) [arXiv:2009.07733]:
`rotational_velocity` 0 (with `magnetic_field_exponent` 0) or 1 rad/s (with
12 or 13). The files have no pinch parameter; α = 2 is used.
"""
@registry_model Kuroda_2020
parameters(::Type{Kuroda_2020}) = (rotational_velocity = Parameter([0, 1]; unit = "rad/s"),
    magnetic_field_exponent = Parameter([0, 12, 13]; label = "B_0 Exponent", description = "Exponent of the magnetic field (Eq. 46)"),
    progenitor_mass = Parameter([20.0]; unit = MSUN))
valid_combination(::Type{Kuroda_2020}, p::NamedTuple) =
    p.rotational_velocity == 1 ? p.magnetic_field_exponent in (12, 13) : p.magnetic_field_exponent == 0
data_files(::Type{Kuroda_2020}, p::NamedTuple) =
    [@sprintf("LnuR%d0B%02d.dat", p.rotational_velocity, p.magnetic_field_exponent)]

function Kuroda_2020(; kwargs...)
    p = validate(Kuroda_2020; kwargs...)
    file = only(data_files(Kuroda_2020, p))
    T = read_named_table(model_file(Kuroda_2020, file))
    cols = Dict{String, Vector{Float64}}("TIME" => T["Tpb[ms]"] ./ 1000)
    for (key, f) in (("nu_e", "NU_E"), ("nu_ae", "NU_E_BAR"), ("nu_x", "NU_X"))
        cols["L_$f"] = T["<L$key>"] .* 1e51
        cols["E_$f"] = T["<E$key>"]
        cols["ALPHA_$f"] = fill(2.0, length(cols["TIME"]))
    end
    return _registry(Kuroda_2020, p, PinchedModel(cols), file; EOS = "LS220")
end

# ----------------------------------------------------------------------------
# Models with axion-like particles: columns 1:t_sim, 2:t_pb (s), 6–8:
# luminosity (erg/s), 9–11: mean energy, 12–14: rms energy (MeV)

function _read_alp_model(path)
    A = read_table(path)
    A = A[vec(all(isfinite, A; dims = 2)) .& (A[:, 1] .> 0), :]
    cols = Dict{String, Vector{Float64}}("TIME" => A[:, 2])
    for (k, f) in enumerate(("NU_E", "NU_E_BAR", "NU_X"))
        cols["L_$f"] = A[:, 5+k]
        cols["E_$f"] = A[:, 8+k]
        cols["ALPHA_$f"] = alpha_from_rms(A[:, 8+k], A[:, 11+k])
    end
    return PinchedModel(cols)
end

const _MORI_PNS_MASS = Dict((0, 0) => 1.78, (100, 2) => 1.77, (100, 4) => 1.76, (100, 10) => 1.77, (100, 12) => 1.77,
    (100, 14) => 1.77, (100, 16) => 1.77, (100, 20) => 1.74, (200, 2) => 1.77, (200, 4) => 1.76, (200, 6) => 1.75,
    (200, 8) => 1.74, (200, 10) => 1.73, (200, 20) => 1.62)

"""
    Mori_2023(; axion_mass, axion_coupling, progenitor_mass = 20)

2D models with axion-like particles of Mori et al., PRD 108, 063027 (2023):
`axion_mass` 0, 100 or 200 MeV and `axion_coupling` in 10⁻¹⁰ GeV⁻¹ (0 for
the standard model; 2, 4, 10, 12, 14, 16, 20 at 100 MeV; 2, 4, 6, 8, 10, 20
at 200 MeV).
"""
@registry_model Mori_2023
parameters(::Type{Mori_2023}) = (axion_mass = Parameter([0, 100, 200]; unit = "MeV", description = "Axion mass"),
    axion_coupling = Parameter([0, 2, 4, 6, 8, 10, 12, 14, 16, 20]; unit = "1e-10/GeV", description = "Axion–photon coupling"),
    progenitor_mass = Parameter([20.0]; unit = MSUN))
valid_combination(::Type{Mori_2023}, p::NamedTuple) = haskey(_MORI_PNS_MASS, (p.axion_mass, p.axion_coupling))
data_files(::Type{Mori_2023}, p::NamedTuple) =
    [p.axion_mass == 0 ? "t-prof_std.dat" : @sprintf("t-prof_%g_%g.dat", p.axion_mass, p.axion_coupling)]
function Mori_2023(; kwargs...)
    p = validate(Mori_2023; kwargs...)
    file = only(data_files(Mori_2023, p))
    return _registry(Mori_2023, p, _read_alp_model(model_file(Mori_2023, file)), file;
                     var"PNS mass" = _MORI_PNS_MASS[(p.axion_mass, p.axion_coupling)])
end

"""
    Takata_2025(; progenitor_mass, axion_mass, axion_coupling)

1D models with axion-like particles of Takata et al., PRD 111, 103028
(2025): `progenitor_mass` 11.2, 20 or 25 M☉; `axion_mass` 0 (with coupling
0) or 40–800 MeV with `axion_coupling` 4, 6, 8, 10 (10⁻¹⁰ GeV⁻¹).
"""
@registry_model Takata_2025
parameters(::Type{Takata_2025}) = (progenitor_mass = Parameter([11.2, 20.0, 25.0]; unit = MSUN),
    axion_mass = Parameter([0, 40, 100, 150, 200, 300, 400, 600, 800]; unit = "MeV", description = "Axion mass"),
    axion_coupling = Parameter([0, 4, 6, 8, 10]; unit = "1e-10/GeV", description = "Axion–photon coupling"))
valid_combination(::Type{Takata_2025}, p::NamedTuple) = (p.axion_mass == 0) == (p.axion_coupling == 0)
function data_files(::Type{Takata_2025}, p::NamedTuple)
    m = @sprintf("%g", trunc(p.progenitor_mass))
    return [p.axion_mass == 0 ? "$(m)_000_00.dat" : @sprintf("%s_%03d_%02d.dat", m, p.axion_mass, p.axion_coupling)]
end
function Takata_2025(; kwargs...)
    p = validate(Takata_2025; kwargs...)
    file = only(data_files(Takata_2025, p))
    return _registry(Takata_2025, p, _read_alp_model(model_file(Takata_2025, file)), file)
end

# ----------------------------------------------------------------------------

"""
    Bugli_2021(; Bfield, direction, grav = nothing, rotation = nothing)

Magnetorotational 3D models of Bugli et al. (2021) [arXiv:2105.00665]:
`Bfield` `"hydro"`, `"L1"` (with `rotation` 0 or 90) or `"L2"` (with `grav`
`"A"` or `"B"`); `direction` `"average"`, `"equator"`, `"north"` or `"south"`.

The files give the mean squared energy ⟨E²⟩, from which the pinch parameter
is computed. (SNEWPY reads that column as the rms energy and squares it again,
which gives α ≈ -1.)
"""
@registry_model Bugli_2021
parameters(::Type{Bugli_2021}) = (Bfield = Parameter(["hydro", "L1", "L2"]),
    direction = Parameter(["average", "equator", "north", "south"]),
    grav = Parameter(["A", "B", nothing]), rotation = Parameter([0, 90, nothing]))
function valid_combination(::Type{Bugli_2021}, p::NamedTuple)
    p.Bfield == "hydro" && return isnothing(p.grav) && isnothing(p.rotation)
    p.Bfield == "L1" && return isnothing(p.grav) && !isnothing(p.rotation)
    return isnothing(p.rotation) && !isnothing(p.grav)
end
function data_files(::Type{Bugli_2021}, p::NamedTuple)
    p.Bfield == "L1" && return ["L1_b12_3d_$(p.rotation)deg_snewpy_$(p.direction).dat"]
    p.Bfield == "L2" && return ["L2_b12_dipdecay_3d_grav$(p.grav)_snewpy_$(p.direction).dat"]
    return ["hydro_3d_snewpy_$(p.direction).dat"]
end
function Bugli_2021(; grav = nothing, rotation = nothing, kwargs...)
    p = validate(Bugli_2021; grav, rotation, kwargs...)
    file = only(data_files(Bugli_2021, p))
    base = PinchedModel(_pinched_from_rms(read_table(model_file(Bugli_2021, file)); rms_is_e2 = true))
    return _registry(Bugli_2021, p, base, file; var"Progenitor mass" = 35.0, EOS = "LS220")
end

# ----------------------------------------------------------------------------

"""
    Fischer_2020()

Model of Fischer et al. (2020) [arXiv:1804.10890] of an 18 M☉ progenitor with
the HS(DD2) equation of state.
"""
@registry_model Fischer_2020
parameters(::Type{Fischer_2020}) = (;)
data_files(::Type{Fischer_2020}, p::NamedTuple) = ["Fischer_2020.tar.gz"]

function Fischer_2020()
    path = model_file(Fischer_2020, "Fischer_2020.tar.gz")
    # luminosity (erg/s), mean and rms energy (MeV) of νe, ν̄e, νx, ν̄x
    table(name) = read_table(IOBuffer(read_tar_member(path, name)); skip = 2)
    L, E, R = table("luminosity.dat"), table("menergy.dat"), table("rmsenergy.dat")
    cols = Dict{String, Vector{Float64}}("TIME" => L[:, 1])
    for (k, f) in enumerate(("NU_E", "NU_E_BAR", "NU_X", "NU_X_BAR"))
        cols["L_$f"] = L[:, 1+k]
        cols["E_$f"] = E[:, 1+k]
        cols["ALPHA_$f"] = alpha_from_rms(E[:, 1+k], R[:, 1+k])
    end
    return Fischer_2020(PinchedModel(cols), OrderedDict{String, Any}("EOS" => "HS(DD2)", "Progenitor mass" => 18.0),
                        "Fischer_2020.tar.gz")
end

# ----------------------------------------------------------------------------

"""
    Analytic3Species(filename)

Pinched model from a table with the columns `TIME`, `L_NU_E`, `L_NU_E_BAR`,
`L_NU_X`, `E_…` and `ALPHA_…` (luminosities in erg/s, mean energies in MeV),
such as the files written by `doc/scripts/Analytic.py` of SNEWPY.
"""
struct Analytic3Species <: SupernovaModel
    base::PinchedModel
    filename::String
end
Analytic3Species(filename::AbstractString) = Analytic3Species(PinchedModel(read_named_table(filename)), filename)
spectra(m::Analytic3Species, t::AbstractVector, E::AbstractVector) = spectra(m.base, t, E)
times(m::Analytic3Species) = times(m.base)
