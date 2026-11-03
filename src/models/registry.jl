# Models of the SNEWS model repositories: their physics parameters, the valid
# parameter combinations and the download of their data files.

"""
    Parameter(values; label = "", description = "", unit = "")

A physics parameter of a model: its allowed `values`, the `label` used in the
metadata, a `description` and the `unit` of numeric values.
"""
struct Parameter
    values::Vector{Any}
    label::String
    description::String
    unit::String
end
Parameter(values; label::AbstractString = "", description::AbstractString = "", unit::AbstractString = "") =
    Parameter(collect(Any, values), label, description, unit)

const _LABELS = Dict(:eos => "EOS", :magnetic_field_exponent => "B_0 Exponent")
const _DESCRIPTIONS = Dict(
    :eos => "Equation of state",
    :progenitor_mass => "Mass of the progenitor star",
    :revival_time => "Time of shock revival",
    :metallicity => "Progenitor metallicity",
    :turbmixing_param => "Turbulent mixing parameter α_Λ",
    :rotational_velocity => "Rotational velocity of the progenitor",
    :magnetic_field_exponent => "Exponent of the magnetic field",
)

_label(name::Symbol, p::Parameter) = !isempty(p.label) ? p.label : get(_LABELS, name, uppercasefirst(replace(String(name), '_' => ' ')))
_description(name::Symbol, p::Parameter) = !isempty(p.description) ? p.description : get(_DESCRIPTIONS, name, _label(name, p))

"""
    RegistryModel

A model of the SNEWS model repositories, created from its physics parameters
(e.g. `Bollig_2016(progenitor_mass = 27)`), which downloads its data files
the first time. See [`parameters`](@ref) and [`param_combinations`](@ref).
"""
abstract type RegistryModel <: SupernovaModel end

spectra(m::RegistryModel, t::AbstractVector, E::AbstractVector) = spectra(m.base, t, E)
times(m::RegistryModel) = times(m.base)
metadata(m::RegistryModel) = m.metadata

# Define the type of a registry model: its data model, metadata and data file
macro registry_model(name)
    return esc(quote
        Base.@__doc__ struct $name <: RegistryModel
            base::SupernovaModel
            metadata::OrderedDict{String, Any}
            filename::String
        end
    end)
end

"""
    parameters(M) -> NamedTuple

The physics parameters of the model type `M` and their allowed values.

```julia
parameters(Nakazato_2013).eos.values   # ["LS220", "shen", "togashi"]
```
"""
function parameters end

# Whether a combination of allowed parameter values has a simulation
valid_combination(::Type{<:RegistryModel}, p::NamedTuple) = true

"""
    param_combinations(M) -> Vector{NamedTuple}

Every valid combination of the parameters of the model type `M`, each of which
can create the model: `M(; p...)`.
"""
function param_combinations(M::Type{<:RegistryModel})
    ps = parameters(M)
    names = keys(ps)
    combos = [NamedTuple{names}(c) for c in Iterators.product((p.values for p in values(ps))...)]
    return filter(p -> valid_combination(M, p), vec(combos))
end

_same(a::Real, b::Real) = isapprox(a, b; rtol = 1e-6, atol = 1e-9)
_same(a, b) = isequal(a, b)

# Check the user parameters and return them in canonical form, with the
# single-valued parameters filled in
function validate(M::Type{<:RegistryModel}; kwargs...)
    ps = parameters(M)
    for k in keys(kwargs)
        k in keys(ps) || throw(SupernovaNeutrinosError(
            "$(nameof(M)) has no parameter $k; its parameters are $(join(keys(ps), ", "))"))
    end
    vals = map(collect(pairs(ps))) do (name, p)
        if !haskey(kwargs, name)
            length(p.values) == 1 && return only(p.values)
            throw(SupernovaNeutrinosError("$(nameof(M)) needs the parameter $name; allowed values: $(_allowed(p))"))
        end
        x = kwargs[name]
        i = findfirst(v -> _same(v, x), p.values)
        isnothing(i) && throw(SupernovaNeutrinosError(
            "$(nameof(M)): $name = $(repr(x)) is not allowed; allowed values: $(_allowed(p))"))
        p.values[i]
    end
    params = NamedTuple{keys(ps)}(Tuple(vals))
    valid_combination(M, params) || throw(SupernovaNeutrinosError(
        "$(nameof(M)) has no simulation for the parameters $params; see param_combinations($(nameof(M)))"))
    return params
end

function _allowed(p::Parameter)
    v = length(p.values) > 12 ? [repr.(p.values[1:6]); "…"; repr.(p.values[end-2:end])] : repr.(p.values)
    return join(v, ", ") * (isempty(p.unit) ? "" : " ($(p.unit))")
end

_metadata(M::Type{<:RegistryModel}, p::NamedTuple) =
    OrderedDict{String, Any}(_label(k, parameters(M)[k]) => v for (k, v) in pairs(p))

# ----------------------------------------------------------------------------
# Model data files

const CCSN_REPOSITORY = "https://github.com/SNEWS2/snewpy-models-ccsn/raw/v0.4/models"
const PRESN_REPOSITORY = "https://github.com/SNEWS2/snewpy-models-presn/raw/v0.2/models"
const PRESN_REPOSITORY_MAIN = "https://github.com/SNEWS2/snewpy-models-presn/raw/main/models"

"URL of the directory with the data files of the model type `M`."
repository(::Type{<:RegistryModel}) = CCSN_REPOSITORY

"""
    model_cache() -> String

Directory where the model data files are downloaded: the environment variable
`SUPERNOVANEUTRINOS_MODELS` if set, otherwise a scratch directory in the
first Julia depot.
"""
model_cache() = get(ENV, "SUPERNOVANEUTRINOS_MODELS",
    joinpath(first(DEPOT_PATH), "scratchspaces", "aa937d89-9fa3-4b83-8f41-d835ba782bde", "models"))

# Local path of a data file of model M, downloaded if missing
function model_file(M::Type{<:RegistryModel}, filename::AbstractString)
    path = joinpath(model_cache(), String(nameof(M)), filename)
    isfile(path) && return path
    mkpath(dirname(path))
    url = repository(M) * "/" * String(nameof(M)) * "/" * join(HTTP_escape.(split(filename, '/')), "/")
    @info "Downloading $url"
    tmp = path * ".part"
    try
        Downloads.download(url, tmp)
    catch e
        rm(tmp; force = true)
        throw(SupernovaNeutrinosError("cannot download $url: $(sprint(showerror, e))"))
    end
    mv(tmp, path; force = true)
    return path
end

HTTP_escape(s::AbstractString) = replace(s, " " => "%20", "#" => "%23")

"""
    data_files(M, p::NamedTuple) -> Vector{String}

Names of the data files of the model type `M` with parameters `p`, relative
to the model directory of its repository.
"""
function data_files end

"""
    download_models(models...)

Download the data files of the model types `models` (all registry models if
none are given) for all their parameter combinations.

```julia
download_models(Bollig_2016, Nakazato_2013)
```
"""
function download_models(models::Type{<:RegistryModel}...)
    isempty(models) && (models = all_models())
    for M in models, p in param_combinations(M), file in data_files(M, p)
        model_file(M, file)
    end
end
