# Neutrino mixing parameters, the vacuum mixing matrix and the flavor
# content of the matter states at very high density.

"""
    MassOrdering

`NormalOrdering` (m1 < m2 < m3) or `InvertedOrdering` (m3 < m1 < m2).
"""
@enum MassOrdering NormalOrdering = 1 InvertedOrdering = 2

function MassOrdering(s::Union{Symbol, AbstractString})
    s = lowercase(String(s))
    s in ("normal", "nmo", "no", "nh", "normalordering") && return NormalOrdering
    s in ("inverted", "imo", "io", "ih", "invertedordering") && return InvertedOrdering
    throw(SupernovaNeutrinosError("unknown mass ordering $s; use :normal or :inverted"))
end
MassOrdering(o::MassOrdering) = o

"Mass ordering implied by the signs of the mass-squared differences (eV²)."
function MassOrdering(dm21_2::Real, dm31_2::Real, dm32_2::Real)
    dm21_2 > 0 || throw(SupernovaNeutrinosError("dm21_2 = $dm21_2 must be positive"))
    dm31_2 * dm32_2 >= 0 || throw(SupernovaNeutrinosError("dm31_2 = $dm31_2 and dm32_2 = $dm32_2 must have the same sign"))
    return dm32_2 >= 0 ? NormalOrdering : InvertedOrdering
end

abstract type AbstractMixingParameters end

"""
    ThreeFlavorMixingParameters(; theta12, theta13, theta23, deltaCP, dm21_2, dm31_2 = nothing, dm32_2 = nothing, mass_order = nothing)

Three-flavor mixing: angles and CP phase in degrees, mass-squared differences
in eV². Give `dm31_2` or `dm32_2` (or both, consistent); the other is derived,
and so is the mass ordering from their sign.

Copy with changes with `ThreeFlavorMixingParameters(p; theta23 = 45)`.
"""
struct ThreeFlavorMixingParameters <: AbstractMixingParameters
    theta12::Float64
    theta13::Float64
    theta23::Float64
    deltaCP::Float64
    dm21_2::Float64
    dm31_2::Float64
    dm32_2::Float64
    mass_order::MassOrdering
end

function ThreeFlavorMixingParameters(;
        theta12::Real, theta13::Real, theta23::Real, deltaCP::Real, dm21_2::Real,
        dm31_2::Union{Real, Nothing} = nothing, dm32_2::Union{Real, Nothing} = nothing,
        mass_order = nothing,
    )
    dm31_2, dm32_2 = _complete_dm2(dm21_2, dm31_2, dm32_2)
    order = isnothing(mass_order) ? MassOrdering(dm21_2, dm31_2, dm32_2) : MassOrdering(mass_order)
    if order == NormalOrdering
        dm31_2 > 0 && dm32_2 > 0 || throw(SupernovaNeutrinosError("dm31_2 and dm32_2 must be positive for the normal ordering"))
    else
        dm31_2 < 0 && dm32_2 < 0 || throw(SupernovaNeutrinosError("dm31_2 and dm32_2 must be negative for the inverted ordering"))
    end
    return ThreeFlavorMixingParameters(theta12, theta13, theta23, deltaCP, dm21_2, dm31_2, dm32_2, order)
end

function _complete_dm2(dm21_2, dm31_2, dm32_2)
    isnothing(dm31_2) && isnothing(dm32_2) && throw(SupernovaNeutrinosError("give dm31_2 or dm32_2"))
    isnothing(dm31_2) && (dm31_2 = dm32_2 + dm21_2)
    isnothing(dm32_2) && (dm32_2 = dm31_2 - dm21_2)
    isapprox(dm32_2 + dm21_2, dm31_2; rtol = 1e-8, atol = 1e-15) ||
        throw(SupernovaNeutrinosError("dm32_2 + dm21_2 = $(dm32_2 + dm21_2) differs from dm31_2 = $dm31_2"))
    return float(dm31_2), float(dm32_2)
end

function ThreeFlavorMixingParameters(p::ThreeFlavorMixingParameters; kwargs...)
    fields = (; theta12 = p.theta12, theta13 = p.theta13, theta23 = p.theta23, deltaCP = p.deltaCP, dm21_2 = p.dm21_2)
    # dm32_2 (and the mass ordering) follow from dm31_2 unless one of them is given
    dm = haskey(kwargs, :dm31_2) || haskey(kwargs, :dm32_2) ? (;) : (; dm31_2 = p.dm31_2)
    return ThreeFlavorMixingParameters(; fields..., dm..., kwargs...)
end

"""
    FourFlavorMixingParameters(p::ThreeFlavorMixingParameters; theta14 = 0, theta24 = 0, theta34 = 0, delta12 = 0, delta24 = 0, dm41_2 = 0)

Three-flavor mixing extended by a sterile neutrino: angles and phases in
degrees, `dm41_2` in eV². `dm42_2` and `dm43_2` follow from the others.
"""
struct FourFlavorMixingParameters <: AbstractMixingParameters
    three::ThreeFlavorMixingParameters
    theta14::Float64
    theta24::Float64
    theta34::Float64
    delta12::Float64
    delta24::Float64
    dm41_2::Float64
    dm42_2::Float64
    dm43_2::Float64
end

function FourFlavorMixingParameters(p::ThreeFlavorMixingParameters = MixingParameters();
        theta14::Real = 0, theta24::Real = 0, theta34::Real = 0, delta12::Real = 0, delta24::Real = 0, dm41_2::Real = 0)
    return FourFlavorMixingParameters(p, theta14, theta24, theta34, delta12, delta24,
                                      dm41_2, dm41_2 - p.dm21_2, dm41_2 - p.dm31_2)
end

# a FourFlavorMixingParameters also answers p.theta12, p.mass_order, …
function Base.getproperty(p::FourFlavorMixingParameters, name::Symbol)
    name in fieldnames(FourFlavorMixingParameters) && return getfield(p, name)
    return getproperty(getfield(p, :three), name)
end
Base.propertynames(::FourFlavorMixingParameters) =
    (fieldnames(ThreeFlavorMixingParameters)..., filter(!=(:three), fieldnames(FourFlavorMixingParameters))...)

"Angles of the mixing matrix in degrees: (θ12, θ13, θ23[, θ14, θ24, θ34])."
mixing_angles(p::ThreeFlavorMixingParameters) = (p.theta12, p.theta13, p.theta23)
mixing_angles(p::FourFlavorMixingParameters) = (mixing_angles(p.three)..., p.theta14, p.theta24, p.theta34)

"Mass-squared differences in eV²: (Δm²21, Δm²31, Δm²32[, Δm²41, Δm²42, Δm²43])."
mass_squared_differences(p::ThreeFlavorMixingParameters) = (p.dm21_2, p.dm31_2, p.dm32_2)
mass_squared_differences(p::FourFlavorMixingParameters) = (mass_squared_differences(p.three)..., p.dm41_2, p.dm42_2, p.dm43_2)

"The mass basis of the mixing parameters (`ThreeFlavorMassBasis` or `FourFlavorMassBasis`)."
mass_basis(::ThreeFlavorMixingParameters) = ThreeFlavorMassBasis
mass_basis(::FourFlavorMixingParameters) = FourFlavorMassBasis
"The flavor basis of the mixing parameters (`ThreeFlavorFlavorBasis` or `FourFlavorFlavorBasis`)."
flavor_basis(::ThreeFlavorMixingParameters) = ThreeFlavorFlavorBasis
flavor_basis(::FourFlavorMixingParameters) = FourFlavorFlavorBasis

# Rotation by θ in the (i, j) plane with phase δ (degrees, PDG sign), on the
# n neutrinos and, complex conjugated, on the n antineutrinos
function _rotation(n::Int, i::Int, j::Int, θ::Real, δ::Real)
    V = Matrix{ComplexF64}(I, 2n, 2n)
    V[i, i] = V[j, j] = cosd(θ)
    V[i, j] = sind(θ) * cis(-deg2rad(δ))
    V[j, i] = -conj(V[i, j])
    V[n+1:end, n+1:end] .= conj.(V[1:n, 1:n])
    return V
end

"""
    vacuum_mixing_matrix(p) -> FlavorMatrix

The PMNS matrix from the mass basis to the flavor basis, for neutrinos and
(complex conjugated) antineutrinos: a 6×6 matrix, or 8×8 with a sterile state.
"""
function vacuum_mixing_matrix(p::ThreeFlavorMixingParameters)
    U = _rotation(3, 2, 3, p.theta23, 0) * _rotation(3, 1, 3, p.theta13, p.deltaCP) * _rotation(3, 1, 2, p.theta12, 0)
    return FlavorMatrix(U, ThreeFlavorFlavorBasis, ThreeFlavorMassBasis)
end

function vacuum_mixing_matrix(p::FourFlavorMixingParameters)
    U = _rotation(4, 3, 4, p.theta34, 0) * _rotation(4, 2, 4, p.theta24, p.delta24) * _rotation(4, 1, 4, p.theta14, 0) *
        _rotation(4, 2, 3, p.theta23, 0) * _rotation(4, 1, 3, p.theta13, p.deltaCP) * _rotation(4, 1, 2, p.theta12, p.delta12)
    return FlavorMatrix(U, FourFlavorFlavorBasis, FourFlavorMassBasis)
end

# Vacuum Hamiltonian (up to a constant and a factor 1/2E) in the flavor basis
function _vacuum_hamiltonian(p::AbstractMixingParameters)
    m2 = [0.0; collect(mass_squared_differences(p))[[1, 2]]]
    p isa FourFlavorMixingParameters && push!(m2, p.dm41_2)
    M2 = FlavorMatrix(Diagonal(ComplexF64[m2; m2]) |> Matrix, mass_basis(p))
    U = vacuum_mixing_matrix(p)
    return U * M2 * U'
end

# Eigenvalues (low, high) of the μ–τ block of H and the flavor content of the
# matter states they belong to: |⟨νμ|ν_low⟩|², |⟨ντ|ν_low⟩|², and the same for ν_high
function _mutau_states(H::FlavorMatrix, bar::Bool)
    μ, τ = bar ? (:mu_bar, :tau_bar) : (:mu, :tau)
    a, c, b = real(H[μ, μ]), real(H[τ, τ]), H[μ, τ]
    T = (a + c) / 2
    # eigenvalues of [[a, b], [b*, c]]: T ± √(|b|² + (a - c)²/4)
    D = sqrt(abs2(b) + (a - c)^2 / 4)
    low, high = T - D, T + D
    content(k) = ((c - k) / (high + low - 2k), (a - k) / (high + low - 2k))
    return (; low = content(low), high = content(high))
end

"""
    pmf_high_density_limit(p) -> FlavorMatrix

Probability that a state of the flavor basis is a given matter state in the
infinite density limit of the supernova: `P[i, α]`, from the flavor basis to
the mass basis. The electron (anti)neutrino is the heaviest (lightest) matter
state of its ordering; the other two share the μ–τ states.
"""
function pmf_high_density_limit(p::ThreeFlavorMixingParameters)
    H = _vacuum_hamiltonian(p)
    ν, ν̄ = _mutau_states(H, false), _mutau_states(H, true)
    P = zero_matrix(ThreeFlavorMassBasis, ThreeFlavorFlavorBasis)
    set!(state, (pμ, pτ), bar = false) = (P[state, bar ? :mu_bar : :mu] = pμ; P[state, bar ? :tau_bar : :tau] = pτ)
    if p.mass_order == NormalOrdering
        # νe is matter state 3, ν̄e is matter state 1̄
        set!(:NU_1, ν.low); set!(:NU_2, ν.high); P[:NU_3, :e] = 1
        P[:NU_1_BAR, :e_bar] = 1; set!(:NU_2_BAR, ν̄.low, true); set!(:NU_3_BAR, ν̄.high, true)
    else
        # νe is matter state 2, ν̄e is matter state 3̄
        set!(:NU_3, ν.low); set!(:NU_1, ν.high); P[:NU_2, :e] = 1
        set!(:NU_1_BAR, ν̄.low, true); set!(:NU_2_BAR, ν̄.high, true); P[:NU_3_BAR, :e_bar] = 1
    end
    return P
end

function pmf_high_density_limit(p::FourFlavorMixingParameters)
    # assumes ν4 is the heaviest state and the electron fraction stays above 1/3
    H = _vacuum_hamiltonian(p)
    ν, ν̄ = _mutau_states(H, false), _mutau_states(H, true)
    P = zero_matrix(FourFlavorMassBasis, FourFlavorFlavorBasis)
    set!(state, (pμ, pτ), bar = false) = (P[state, bar ? :mu_bar : :mu] = pμ; P[state, bar ? :tau_bar : :tau] = pτ)
    if p.mass_order == NormalOrdering
        # νe is 4, νs is 3; ν̄e is 1̄, ν̄s is 2̄
        set!(:NU_1, ν.low); set!(:NU_2, ν.high); P[:NU_3, :s] = 1; P[:NU_4, :e] = 1
        P[:NU_1_BAR, :e_bar] = 1; P[:NU_2_BAR, :s_bar] = 1; set!(:NU_3_BAR, ν̄.low, true); set!(:NU_4_BAR, ν̄.high, true)
    else
        # νe is 4, νs is 2; ν̄e is 3̄, ν̄s is 1̄
        set!(:NU_3, ν.low); set!(:NU_1, ν.high); P[:NU_2, :s] = 1; P[:NU_4, :e] = 1
        P[:NU_1_BAR, :s_bar] = 1; set!(:NU_2_BAR, ν̄.low, true); P[:NU_3_BAR, :e_bar] = 1; set!(:NU_4_BAR, ν̄.high, true)
    end
    return P
end

# Mixing parameter presets: (θ12, θ13, θ23, δCP) in degrees, Δm²21 and Δm²31 or Δm²32 in eV²
const MIXING_PRESETS = Dict(
    # NuFIT 5.0, JHEP 09 (2020) 178 [arXiv:2007.14792]
    "NuFIT5.0" => Dict(
        NormalOrdering => (theta12 = 33.44, theta13 = 8.57, theta23 = 49.2, deltaCP = 197, dm21_2 = 7.42e-5, dm31_2 = 2.517e-3),
        InvertedOrdering => (theta12 = 33.45, theta13 = 8.60, theta23 = 49.3, deltaCP = 282, dm21_2 = 7.42e-5, dm32_2 = -2.498e-3),
    ),
    # NuFIT 5.2 (2022), www.nu-fit.org
    "NuFIT5.2" => Dict(
        NormalOrdering => (theta12 = 33.41, theta13 = 8.58, theta23 = 42.2, deltaCP = 232, dm21_2 = 7.41e-5, dm31_2 = 2.507e-3),
        InvertedOrdering => (theta12 = 33.41, theta13 = 8.57, theta23 = 49.0, deltaCP = 276, dm21_2 = 7.41e-5, dm32_2 = -2.486e-3),
    ),
    # NuFIT 6.0, JHEP 12 (2024) 216 [arXiv:2410.05380]
    "NuFIT6.0" => Dict(
        NormalOrdering => (theta12 = 33.68, theta13 = 8.56, theta23 = 43.3, deltaCP = 212, dm21_2 = 7.49e-5, dm31_2 = 2.513e-3),
        InvertedOrdering => (theta12 = 33.68, theta13 = 8.59, theta23 = 47.9, deltaCP = 274, dm21_2 = 7.49e-5, dm32_2 = -2.484e-3),
    ),
    # R.L. Workman et al. (Particle Data Group), PTEP 2022, 083C01
    "PDG2022" => Dict(
        NormalOrdering => (theta12 = 33.65, theta13 = 8.53, theta23 = 47.64, deltaCP = 245, dm21_2 = 7.53e-5, dm32_2 = 2.453e-3),
        InvertedOrdering => (theta12 = 33.65, theta13 = 8.53, theta23 = 47.24, deltaCP = 245, dm21_2 = 7.53e-5, dm32_2 = -2.536e-3),
    ),
    # S. Navas et al. (Particle Data Group), PRD 110, 030001 (2024)
    "PDG2024" => Dict(
        NormalOrdering => (theta12 = 33.65, theta13 = 8.51, theta23 = 48.33, deltaCP = 214, dm21_2 = 7.53e-5, dm32_2 = 2.455e-3),
        InvertedOrdering => (theta12 = 33.65, theta13 = 8.51, theta23 = 48.04, deltaCP = 214, dm21_2 = 7.53e-5, dm32_2 = -2.529e-3),
    ),
)

"""
    MixingParameters(mass_order = :normal; version = "NuFIT6.0") -> ThreeFlavorMixingParameters

Best-fit three-flavor mixing parameters of a global fit: `version` is one of
`"NuFIT5.0"`, `"NuFIT5.2"`, `"NuFIT6.0"`, `"PDG2022"`, `"PDG2024"`;
`mass_order` is `:normal`/`:inverted` or a [`MassOrdering`](@ref).
"""
function MixingParameters(mass_order = NormalOrdering; version::AbstractString = "NuFIT6.0")
    haskey(MIXING_PRESETS, version) || throw(SupernovaNeutrinosError(
        "unknown mixing parameter version $version; choose from $(join(sort(collect(keys(MIXING_PRESETS))), ", "))"))
    order = MassOrdering(mass_order)
    return ThreeFlavorMixingParameters(; MIXING_PRESETS[version][order]..., mass_order = order)
end

function Base.show(io::IO, ::MIME"text/plain", p::ThreeFlavorMixingParameters)
    println(io, "ThreeFlavorMixingParameters (", p.mass_order, ")")
    @printf(io, "  θ12 = %g°, θ13 = %g°, θ23 = %g°, δCP = %g°\n", p.theta12, p.theta13, p.theta23, p.deltaCP)
    @printf(io, "  Δm²21 = %g eV², Δm²31 = %g eV², Δm²32 = %g eV²", p.dm21_2, p.dm31_2, p.dm32_2)
end

function Base.show(io::IO, ::MIME"text/plain", p::FourFlavorMixingParameters)
    show(io, MIME"text/plain"(), p.three)
    @printf(io, "\n  θ14 = %g°, θ24 = %g°, θ34 = %g°, δ12 = %g°, δ24 = %g°\n", p.theta14, p.theta24, p.theta34, p.delta12, p.delta24)
    @printf(io, "  Δm²41 = %g eV², Δm²42 = %g eV², Δm²43 = %g eV²", p.dm41_2, p.dm42_2, p.dm43_2)
end
