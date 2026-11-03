# Flavor schemes: ordered sets of neutrino states such as
# (νe, ν̄e, νμ, ν̄μ, ντ, ν̄τ) or the mass states (ν1, ν2, ν3, ν̄1, ν̄2, ν̄3).

"""
    FlavorScheme(name, states)
    FlavorScheme(name, leptons; interleaved = true)

An ordered set of neutrino states, each named `:NU_<lepton>` or
`:NU_<lepton>_BAR`. The predefined schemes are [`TwoFlavor`](@ref)
(e, x), [`ThreeFlavor`](@ref) (e, μ, τ), [`FourFlavor`](@ref) (e, μ, τ, s)
and the mass and flavor bases used by the mixing matrices.

Look up a state by name, in any case and with or without the `NU_` prefix:

```julia
ThreeFlavor.NU_E            # the νe of ThreeFlavor
ThreeFlavor[:mu_bar]        # ν̄μ
ThreeFlavor["NU"]           # all the neutrinos (not antineutrinos)
ThreeFlavor[[:e, :e_bar]]   # a vector of flavors
```
"""
struct FlavorScheme
    name::Symbol
    states::Vector{Symbol}
end

function FlavorScheme(name::Symbol, leptons::AbstractVector{<:AbstractString}; interleaved::Bool = true)
    states = if interleaved
        [Symbol("NU_", L, bar) for L in leptons for bar in ("", "_BAR")]
    else
        [Symbol("NU_", L, bar) for bar in ("", "_BAR") for L in leptons]
    end
    return FlavorScheme(name, states)
end

Base.:(==)(a::FlavorScheme, b::FlavorScheme) = getfield(a, :name) == getfield(b, :name) && getfield(a, :states) == getfield(b, :states)
Base.hash(s::FlavorScheme, h::UInt) = hash(getfield(s, :states), hash(getfield(s, :name), h))

"""
    Flavor

One state of a [`FlavorScheme`](@ref), e.g. `ThreeFlavor.NU_E`.
"""
struct Flavor
    scheme::FlavorScheme
    index::Int
end

Base.:(==)(a::Flavor, b::Flavor) = a.index == b.index && a.scheme == b.scheme
Base.hash(f::Flavor, h::UInt) = hash(f.index, hash(f.scheme, h))
Base.isless(a::Flavor, b::Flavor) = isless(a.index, b.index)
Base.Symbol(f::Flavor) = getfield(f.scheme, :states)[f.index]
Base.string(f::Flavor) = String(Symbol(f))
Base.show(io::IO, f::Flavor) = print(io, getfield(f.scheme, :name), ".", Symbol(f))

Base.length(s::FlavorScheme) = length(getfield(s, :states))
Base.firstindex(s::FlavorScheme) = 1
Base.lastindex(s::FlavorScheme) = length(s)
Base.eachindex(s::FlavorScheme) = Base.OneTo(length(s))
Base.iterate(s::FlavorScheme, i = 1) = i > length(s) ? nothing : (Flavor(s, i), i + 1)
Base.eltype(::Type{FlavorScheme}) = Flavor
Base.collect(s::FlavorScheme) = [Flavor(s, i) for i in eachindex(s)]
Base.show(io::IO, s::FlavorScheme) = print(io, getfield(s, :name))
function Base.show(io::IO, ::MIME"text/plain", s::FlavorScheme)
    print(io, getfield(s, :name), ": ", join(getfield(s, :states), ", "))
end

# ThreeFlavor.NU_E: the fields are reached with getfield
function Base.getproperty(s::FlavorScheme, name::Symbol)
    name === :name && return getfield(s, :name)
    name === :states && return getfield(s, :states)
    return s[name]
end
Base.propertynames(s::FlavorScheme) = (:name, :states, getfield(s, :states)...)

Base.getindex(s::FlavorScheme, i::Integer) =
    1 <= i <= length(s) ? Flavor(s, i) : throw(BoundsError(s, i))
Base.getindex(s::FlavorScheme, ::Colon) = collect(s)
Base.getindex(s::FlavorScheme, keys::AbstractVector) = [s[k] for k in keys]
Base.getindex(s::FlavorScheme, keys::Tuple) = [s[k] for k in keys]
function Base.getindex(s::FlavorScheme, f::Flavor)
    f.scheme == s || throw(SupernovaNeutrinosError("$f is not a state of the $(s.name) scheme"))
    return f
end
Base.getindex(s::FlavorScheme, key::Symbol) = s[String(key)]

# Look a state up by name; "NU" and "NU_BAR" give all the neutrinos or antineutrinos
function Base.getindex(s::FlavorScheme, key::AbstractString)
    key = uppercase(key)
    key == "NU" && return [f for f in s if is_neutrino(f)]
    key == "NU_BAR" && return [f for f in s if is_antineutrino(f)]
    startswith(key, "NU_") || (key = "NU_" * key)
    i = findfirst(==(Symbol(key)), getfield(s, :states))
    isnothing(i) && throw(SupernovaNeutrinosError(
        "no state $key in the $(s.name) scheme; the states are $(join(s.states, ", "))"))
    return Flavor(s, i)
end

"Index of each state in the scheme, for a key accepted by `getindex(::FlavorScheme, key)`."
flavor_index(s::FlavorScheme, key) = _index(s[key])
flavor_index(s::FlavorScheme, ::Colon) = Colon()
_index(f::Flavor) = f.index
_index(fs::AbstractVector{Flavor}) = [f.index for f in fs]

"`true` for an antineutrino state (name ending with `_BAR`)."
is_antineutrino(f::Flavor) = endswith(string(f), "_BAR")
"`true` for a neutrino state (not an antineutrino)."
is_neutrino(f::Flavor) = !is_antineutrino(f)
"The lepton of a state: `\"E\"`, `\"MU\"`, `\"TAU\"`, `\"S\"`, `\"X\"`, or the mass index `\"1\"`, `\"2\"`…"
lepton(f::Flavor) = split(string(f), '_')[2]
"`true` for νe and ν̄e."
is_electron(f::Flavor) = lepton(f) == "E"
"`true` for νμ and ν̄μ."
is_muon(f::Flavor) = lepton(f) == "MU"
"`true` for ντ and ν̄τ."
is_tauon(f::Flavor) = lepton(f) == "TAU"
"`true` for sterile states."
is_sterile(f::Flavor) = lepton(f) == "S"

"""
    to_tex(f::Flavor) -> String

LaTeX symbol of the state, e.g. `\$\\overline{\\nu}_{\\mu}\$`.
"""
function to_tex(f::Flavor)
    base = is_antineutrino(f) ? "\\overline{\\nu}" : "\\nu"
    l = lowercase(lepton(f))
    (is_muon(f) || is_tauon(f)) && (l = "\\" * l)
    return "\$$(base)_{$l}\$"
end

"Two flavors: νe, ν̄e, νx, ν̄x (x stands for μ and τ)."
const TwoFlavor = FlavorScheme(:TwoFlavor, ["E", "X"])
"Three flavors: νe, ν̄e, νμ, ν̄μ, ντ, ν̄τ. The scheme of every supernova model."
const ThreeFlavor = FlavorScheme(:ThreeFlavor, ["E", "MU", "TAU"])
"Four flavors: νe, ν̄e, νμ, ν̄μ, ντ, ν̄τ, νs, ν̄s."
const FourFlavor = FlavorScheme(:FourFlavor, ["E", "MU", "TAU", "S"])

"Mass basis of three-flavor mixing: ν1, ν2, ν3, ν̄1, ν̄2, ν̄3."
const ThreeFlavorMassBasis = FlavorScheme(:ThreeFlavorMassBasis, ["1", "2", "3"]; interleaved = false)
"Flavor basis of three-flavor mixing: νe, νμ, ντ, ν̄e, ν̄μ, ν̄τ."
const ThreeFlavorFlavorBasis = FlavorScheme(:ThreeFlavorFlavorBasis, ["E", "MU", "TAU"]; interleaved = false)
"Mass basis of four-flavor mixing: ν1…ν4, ν̄1…ν̄4."
const FourFlavorMassBasis = FlavorScheme(:FourFlavorMassBasis, ["1", "2", "3", "4"]; interleaved = false)
"Flavor basis of four-flavor mixing: νe, νμ, ντ, νs, ν̄e, ν̄μ, ν̄τ, ν̄s."
const FourFlavorFlavorBasis = FlavorScheme(:FourFlavorFlavorBasis, ["E", "MU", "TAU", "S"]; interleaved = false)
