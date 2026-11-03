# Flavor transformations between the production of the neutrinos in the
# supernova and their detection on Earth.
#
# A TransformationChain multiplies three steps:
#   P(να → νβ) = Σij P_fm[β, j] · P_mm[j, i] · P_mf[i, α]
# in the supernova (flavor → matter state, P_mf), in vacuum (mass → mass
# state, P_mm) and in the Earth (mass state → flavor, P_fm).

"""
    FlavorTransformation

A transformation of the neutrino flavors, giving the transition probability
matrix `P_ff(x, t, E)` from the flavors at the source to those at the
detector. Apply it to a [`Container`](@ref) with [`apply`](@ref).
"""
abstract type FlavorTransformation end

"""
    SNTransformation

Transformation inside the supernova: `P_mf(x, params, t, E)` from the flavor
basis to the matter states leaving the star.
"""
abstract type SNTransformation end

"""
    VacuumTransformation

Transformation in vacuum: `P_mm(x, params, t, E)` between mass states.
"""
abstract type VacuumTransformation end

"""
    EarthTransformation

Transformation in the Earth: `P_fm(x, params, t, E)` from mass states to the
flavor basis.
"""
abstract type EarthTransformation end

"""
    P_ff(x, t, E) -> FlavorMatrix

Probability matrix `P[β, α]` of the transition να → νβ for times `t` (s) and
energies `E` (MeV).
"""
function P_ff end

"""
    P_mf(x::SNTransformation, params, t, E) -> FlavorMatrix

Probability that a flavor state at production leaves the supernova as a given
matter (mass) state.
"""
function P_mf end

"""
    P_mm(x::VacuumTransformation, params, t, E) -> FlavorMatrix

Transition probability between mass states during propagation in vacuum.
"""
function P_mm end

"""
    P_fm(x::EarthTransformation, params, t, E) -> FlavorMatrix

Probability that a mass state reaching the Earth is detected as a given flavor.
"""
function P_fm end

"""
    apply(x, c::Container) -> Container

Apply the flavor transformation `x` to the spectrum or flux `c`. The transition
probabilities are evaluated at the time and energy points of `c` (bin centers
on binned axes).
"""
function apply(x::FlavorTransformation, c::Container)
    M = P_ff(x, _points(c, :time), _points(c, :energy))
    return (c.scheme * M * c.scheme) * c
end

"""
    NoTransformation()

No flavor transformation: the flux at Earth has the flavors of the source.
"""
struct NoTransformation <: FlavorTransformation end
P_ff(::NoTransformation, t, E) = FlavorMatrix(I, ThreeFlavor)
apply(::NoTransformation, c::Container) = c

"""
    CompleteExchange()

Each electron flavor goes half to μ and half to τ, and each μ or τ flavor goes
half to each of the other two flavors.
"""
struct CompleteExchange <: FlavorTransformation end
P_ff(::CompleteExchange, t, E) = FlavorMatrix((β, α) -> 0.5 * (is_neutrino(β) == is_neutrino(α)) * (β != α), ThreeFlavor)

# ----------------------------------------------------------------------------
# Transformations in the supernova

"""
    AdiabaticMSW()

Adiabatic MSW effect: each flavor leaves the supernova in the matter state it
occupies at very high density ([`pmf_high_density_limit`](@ref)).
"""
struct AdiabaticMSW <: SNTransformation end
P_mf(::AdiabaticMSW, p::AbstractMixingParameters, t, E) = pmf_high_density_limit(p)

function _swap_rows!(P::FlavorMatrix, from::Vector{Symbol}, to::Vector{Symbol})
    P.array[flavor_index(P.flavor_out, from), :, :, :] = P.array[flavor_index(P.flavor_out, to), :, :, :]
    return P
end

"""
    NonAdiabaticMSWH()

Non-adiabatic MSW H resonance, e.g. when a shock passes through the resonance
(Schirato & Fuller 2002). The resonance mixes ν2 and ν3 for the normal
ordering, ν̄1 and ν̄3 for the inverted ordering.
"""
struct NonAdiabaticMSWH <: SNTransformation end
function P_mf(::NonAdiabaticMSWH, p::ThreeFlavorMixingParameters, t, E)
    P = pmf_high_density_limit(p)
    p.mass_order == NormalOrdering ? _swap_rows!(P, [:NU_2, :NU_3], [:NU_3, :NU_2]) :
                                     _swap_rows!(P, [:NU_1_BAR, :NU_3_BAR], [:NU_3_BAR, :NU_1_BAR])
end

"""
    TwoFlavorDecoherence()

Equal mixing of the two matter states of the H resonance, for density
fluctuations ≲ 10% around the resonance density (Kneller 2010; Kneller &
Mauney 2013): ν2 and ν3 (normal ordering) or ν̄1 and ν̄3 (inverted ordering).
"""
struct TwoFlavorDecoherence <: SNTransformation end
function P_mf(::TwoFlavorDecoherence, p::ThreeFlavorMixingParameters, t, E)
    P = pmf_high_density_limit(p)
    a, b = p.mass_order == NormalOrdering ? (:NU_2, :NU_3) : (:NU_1_BAR, :NU_3_BAR)
    i, j = flavor_index(P.flavor_out, a), flavor_index(P.flavor_out, b)
    avg = (P.array[i, :, :, :] .+ P.array[j, :, :, :]) ./ 2
    P.array[i, :, :, :] = avg
    P.array[j, :, :, :] = avg
    return P
end

"""
    ThreeFlavorDecoherence()

Equal mixing of the three neutrino states and of the three antineutrino
states, for density fluctuations ≳ 10% around the H resonance density.
"""
struct ThreeFlavorDecoherence <: SNTransformation end
P_mf(::ThreeFlavorDecoherence, ::ThreeFlavorMixingParameters, t, E) =
    FlavorMatrix((i, α) -> (is_neutrino(i) == is_neutrino(α)) / 3, ThreeFlavorMassBasis, ThreeFlavorFlavorBasis)

"""
    AdiabaticMSWes()

Four-neutrino mixing with ν4 the heaviest state (but with the electron–sterile
resonances outside the neutrinosphere), an adiabatic outer electron–sterile
resonance and a non-adiabatic inner one (where the electron fraction is 1/3).
See Esmaili, Peres and Serpico, PRD 90, 033013 (2014).
"""
struct AdiabaticMSWes <: SNTransformation end
P_mf(::AdiabaticMSWes, p::FourFlavorMixingParameters, t, E) = pmf_high_density_limit(p)

"""
    NonAdiabaticMSWes()

As [`AdiabaticMSWes`](@ref), with both electron–sterile resonances
non-adiabatic.
"""
struct NonAdiabaticMSWes <: SNTransformation end
function P_mf(::NonAdiabaticMSWes, p::FourFlavorMixingParameters, t, E)
    P = pmf_high_density_limit(p)
    if p.mass_order == NormalOrdering
        _swap_rows!(P, [:NU_3, :NU_4], [:NU_4, :NU_3])
        _swap_rows!(P, [:NU_2_BAR, :NU_3_BAR, :NU_4_BAR], [:NU_3_BAR, :NU_4_BAR, :NU_2_BAR])
    else
        _swap_rows!(P, [:NU_2, :NU_4], [:NU_4, :NU_2])
        _swap_rows!(P, [:NU_1_BAR, :NU_2_BAR, :NU_4_BAR], [:NU_2_BAR, :NU_4_BAR, :NU_1_BAR])
    end
    return P
end

"""
    FourFlavorDecoherence()

Equal mixing of the four neutrino states and of the four antineutrino states.
"""
struct FourFlavorDecoherence <: SNTransformation end
P_mf(::FourFlavorDecoherence, ::FourFlavorMixingParameters, t, E) =
    FlavorMatrix((i, α) -> (is_neutrino(i) == is_neutrino(α)) / 4, FourFlavorMassBasis, FourFlavorFlavorBasis)

"Default mixing parameters for a transformation in the supernova."
default_mixing(::SNTransformation) = MixingParameters()
default_mixing(::Union{AdiabaticMSWes, NonAdiabaticMSWes, FourFlavorDecoherence}) = FourFlavorMixingParameters()

# ----------------------------------------------------------------------------
# Transformations in vacuum

"""
    NoVacuumTransformation()

Mass states propagate unchanged in vacuum.
"""
struct NoVacuumTransformation <: VacuumTransformation end
P_mm(::NoVacuumTransformation, p::AbstractMixingParameters, t, E) = FlavorMatrix(I, mass_basis(p))

"""
    NeutrinoDecay(; mass = 1.0, lifetime = 86400.0, distance = 10.0)

Decay of the heaviest neutrino to the lightest one: `mass` of the heaviest
state in eV, its `lifetime` in s and the `distance` travelled in kpc. See
A. de Gouvêa et al., PRD 101, 043013 (2020) [arXiv:1910.01127].
"""
Base.@kwdef struct NeutrinoDecay <: VacuumTransformation
    mass::Float64 = 1.0
    lifetime::Float64 = DAY
    distance::Float64 = 10.0
end

"Survival probability of the decaying state at energies `E` (MeV)."
function decay_factor(x::NeutrinoDecay, E)
    # Γ·d = (m c² / E) · d / (c τ)
    return @. exp(-(x.mass / (E * 1e6)) * x.distance * KPC_IN_CM / (SPEED_OF_LIGHT * x.lifetime))
end

function P_mm(x::NeutrinoDecay, p::ThreeFlavorMixingParameters, t, E)
    d = reshape(decay_factor(x, E), 1, :)
    P = FlavorMatrix(I, ThreeFlavorMassBasis; ne = length(E))
    for bar in ("", "_BAR")
        st(i) = Symbol("NU_", i, bar)
        if p.mass_order == NormalOrdering
            # ν3 → ν1
            P[st(1), st(3)] = 1 .- d
            P[st(3), st(3)] = d
        else
            # ν2 → ν3
            P[st(2), st(2)] = d
            P[st(3), st(2)] = 1 .- d
        end
    end
    return P
end

"""
    QuantumDecoherence(; Gamma3 = 1e-27, Gamma8 = 1e-27, distance = 10.0, n = 0, E0 = 10.0)

Quantum decoherence in vacuum, which drives the mass states to equipartition:
decoherence parameters `Gamma3`, `Gamma8` in eV, `distance` in kpc, and the
energy dependence Γ(E) = Γ·(E/E0)ⁿ with `E0` in MeV. See M. V. dos Santos et
al. (2023) [arXiv:2306.17591].
"""
Base.@kwdef struct QuantumDecoherence <: VacuumTransformation
    Gamma3::Float64 = 1e-27
    Gamma8::Float64 = 1e-27
    distance::Float64 = 10.0
    n::Float64 = 0.0
    E0::Float64 = 10.0
end

function P_mm(x::QuantumDecoherence, p::ThreeFlavorMixingParameters, t, E)
    L = x.distance * KPC_IN_CM / HBARC_EV_CM            # distance in 1/eV
    xL = @. (E / x.E0)^x.n * L
    e38 = @. exp(-(x.Gamma3 + x.Gamma8 / 3) * xL)
    e8 = @. exp(-x.Gamma8 * xL)
    row(v) = reshape(v, 1, :)
    P11 = row(1 / 3 .+ e38 ./ 2 .+ e8 ./ 6)
    P12 = row(1 / 3 .- e38 ./ 2 .+ e8 ./ 6)
    P13 = row(1 / 3 .- e8 ./ 3)
    P33 = row(1 / 3 .+ 2 .* e8 ./ 3)
    P = zero_matrix(ThreeFlavorMassBasis; ne = length(E))
    for bar in ("", "_BAR")
        st(i) = Symbol("NU_", i, bar)
        P[st(1), st(1)] = P11; P[st(2), st(2)] = P11
        P[st(1), st(2)] = P12; P[st(2), st(1)] = P12
        for i in (1, 2)
            P[st(i), st(3)] = P13; P[st(3), st(i)] = P13
        end
        P[st(3), st(3)] = P33
    end
    return P
end

# ----------------------------------------------------------------------------
# Transformations in the Earth

"""
    NoEarthMatter()

No Earth matter effect: a mass state is detected as flavor α with
probability |Uαi|².
"""
struct NoEarthMatter <: EarthTransformation end
P_fm(::NoEarthMatter, p::AbstractMixingParameters, t, E) = abs2(vacuum_mixing_matrix(p))

# ----------------------------------------------------------------------------

"""
    TransformationChain(steps...; mixing_params)

Flavor transformation made of one step in the supernova
([`SNTransformation`](@ref)) and optionally one in vacuum
([`VacuumTransformation`](@ref), default [`NoVacuumTransformation`](@ref))
and one in the Earth ([`EarthTransformation`](@ref), default
[`NoEarthMatter`](@ref)), in any order. All steps use `mixing_params`, by
default [`MixingParameters`](@ref)`()` (NuFIT 6.0, normal ordering), or its
four-flavor extension for the sterile-neutrino steps.

```julia
TransformationChain(AdiabaticMSW(); mixing_params = MixingParameters(:inverted))
TransformationChain(AdiabaticMSW(), NeutrinoDecay(lifetime = 1e5))
```

A step in the supernova used on its own is the chain with the default steps
and parameters: `apply(AdiabaticMSW(), c)`.
"""
struct TransformationChain{S <: SNTransformation, V <: VacuumTransformation, E <: EarthTransformation, M <: AbstractMixingParameters} <: FlavorTransformation
    in_sn::S
    in_vacuum::V
    in_earth::E
    mixing_params::M
end

function TransformationChain(steps...; mixing_params = nothing)
    pick(T, default) = (found = filter(s -> s isa T, steps); length(found) > 1 ?
        throw(SupernovaNeutrinosError("a chain has one step of type $T, not $(length(found))")) :
        isempty(found) ? default : only(found))
    sn = pick(SNTransformation, nothing)
    isnothing(sn) && throw(SupernovaNeutrinosError("a transformation chain needs a step in the supernova (SNTransformation)"))
    vacuum = pick(VacuumTransformation, NoVacuumTransformation())
    earth = pick(EarthTransformation, NoEarthMatter())
    all(s -> s isa Union{SNTransformation, VacuumTransformation, EarthTransformation}, steps) ||
        throw(SupernovaNeutrinosError("the steps of a chain are SN, vacuum and Earth transformations"))
    params = isnothing(mixing_params) ? default_mixing(sn) : mixing_params
    return TransformationChain(sn, vacuum, earth, params)
end

P_ff(x::TransformationChain, t, E) =
    P_fm(x.in_earth, x.mixing_params, t, E) * P_mm(x.in_vacuum, x.mixing_params, t, E) * P_mf(x.in_sn, x.mixing_params, t, E)

P_ff(x::SNTransformation, t, E) = P_ff(TransformationChain(x), t, E)
apply(x::SNTransformation, c::Container) = apply(TransformationChain(x), c)

function Base.show(io::IO, x::TransformationChain)
    print(io, nameof(typeof(x.in_sn)), "+", nameof(typeof(x.in_vacuum)), "+", nameof(typeof(x.in_earth)), "_",
          x.mixing_params.mass_order)
end
