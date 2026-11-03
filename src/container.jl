# A container for a quantity that depends on flavor, time and energy: the
# spectrum emitted by the supernova, the flux at Earth, the fluence, event
# rates…

"""
    ContainerUnit(MeV, s, cm)

Unit of a [`Container`](@ref) as powers of MeV, s and cm: the flux is
`ContainerUnit(-1, -1, -2)`, i.e. 1/(MeV s cm²).
"""
struct ContainerUnit
    MeV::Int
    s::Int
    cm::Int
end

Base.:*(a::ContainerUnit, b::ContainerUnit) = ContainerUnit(a.MeV + b.MeV, a.s + b.s, a.cm + b.cm)

const _SUPERSCRIPTS = Dict(zip("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
_power(name, n) = n == 0 ? "" : n == 1 ? name : name * map(c -> _SUPERSCRIPTS[c], string(n))

function Base.show(io::IO, u::ContainerUnit)
    s = join(filter(!isempty, [_power("MeV", u.MeV), _power("s", u.s), _power("cm", u.cm)]), " ")
    print(io, isempty(s) ? "1" : s)
end

const SPECTRUM_UNIT = ContainerUnit(-1, -1, 0)
const FLUX_UNIT = ContainerUnit(-1, -1, -2)
const _KINDS = Dict(
    ContainerUnit(-1, -1, 0) => "spectrum or differential event rate",
    ContainerUnit(-1, -1, -2) => "flux",
    ContainerUnit(-1, 0, -2) => "fluence",
    ContainerUnit(0, -1, -2) => "integrated flux",
    ContainerUnit(0, 0, -2) => "integrated fluence",
    ContainerUnit(-1, 0, 0) => "event spectrum",
    ContainerUnit(0, -1, 0) => "event rate",
    ContainerUnit(0, 0, 0) => "events",
)

const _AXES = (flavor = 1, time = 2, energy = 3)
_axis(axis::Symbol) = haskey(_AXES, axis) ? _AXES[axis] :
    throw(SupernovaNeutrinosError("unknown axis $axis; the axes are :flavor, :time and :energy"))
_axis(axis::AbstractString) = _axis(Symbol(axis))

"""
    Container(array, flavor, time, energy; unit, scheme, integrable)

A quantity that depends on flavor, time (s) and energy (MeV): `array` has
dimensions `(flavor, time, energy)`.

Along time and energy, the axis values are either sampling points
(`length(time) == size(array, 2)`) or bin edges
(`length(time) == size(array, 2) + 1`). Sampled axes can be integrated with
[`integrate`](@ref), binned axes summed with `sum`.

`unit` is a [`ContainerUnit`](@ref); `scheme` the [`FlavorScheme`](@ref) of
the flavors (by default the scheme of `flavor`).

Select flavors with `c[ThreeFlavor.NU_E]`, `c[:e_bar]` or `c[:NU]`; convert to
another scheme with `TwoFlavor * c` or [`convert_flavor`](@ref).
"""
struct Container
    array::Array{Float64, 3}
    flavor::Vector{Flavor}
    scheme::FlavorScheme
    time::Vector{Float64}
    energy::Vector{Float64}
    integrable::NTuple{2, Bool}
    unit::ContainerUnit
    function Container(array::AbstractArray{<:Real, 3}, flavor::AbstractVector{Flavor}, time::AbstractVector{<:Real},
                       energy::AbstractVector{<:Real}; unit::ContainerUnit,
                       scheme::FlavorScheme = _scheme(flavor), integrable = nothing)
        nf, nt, ne = size(array)
        nf in (length(flavor), length(flavor) - 1) || throw(SupernovaNeutrinosError(
            "array with $nf flavors does not match the $(length(flavor)) flavors"))
        nt in (length(time), length(time) - 1) || throw(SupernovaNeutrinosError(
            "array with $nt times does not match $(length(time)) time values (sampling points or bin edges)"))
        ne in (length(energy), length(energy) - 1) || throw(SupernovaNeutrinosError(
            "array with $ne energies does not match $(length(energy)) energy values (sampling points or bin edges)"))
        integ = isnothing(integrable) ? (nt == length(time), ne == length(energy)) : Tuple(integrable)
        return new(Array{Float64}(array), collect(flavor), scheme, collect(float.(time)), collect(float.(energy)), integ, unit)
    end
end

function _scheme(flavor::AbstractVector{Flavor})
    schemes = unique(f.scheme for f in flavor)
    length(schemes) == 1 || throw(SupernovaNeutrinosError("flavors must come from one scheme, not $(join(schemes, ", "))"))
    return only(schemes)
end

"""
    Container(dict, time, energy; unit, scheme = ThreeFlavor)

Container from a dictionary `flavor => (time, energy) array` with one entry per
state of `scheme`.
"""
function Container(dict::AbstractDict{Flavor}, time::AbstractVector, energy::AbstractVector;
                   unit::ContainerUnit, scheme::FlavorScheme = ThreeFlavor, integrable = nothing)
    A = zeros(length(scheme), size(first(values(dict)))...)
    for f in scheme
        A[f.index, :, :] = dict[f]
    end
    return Container(A, collect(scheme), time, energy; unit, scheme, integrable)
end

_with(c::Container; array = c.array, flavor = c.flavor, time = c.time, energy = c.energy,
      integrable = c.integrable, unit = c.unit, scheme = c.scheme) =
    Container(array, flavor, time, energy; unit, scheme, integrable)

Base.size(c::Container) = size(c.array)
Base.size(c::Container, d) = size(c.array, d)

"`true` if the axis (`:time` or `:energy`) holds sampling points and can be integrated."
can_integrate(c::Container, axis) = (a = _axis(axis); a != 1 && c.integrable[a-1])
"`true` if the axis holds bins (or flavors) and can be summed."
can_sum(c::Container, axis) = !can_integrate(c, axis)

"Centers of the bins with edges `x`."
bin_centers(x::AbstractVector) = (x[1:end-1] .+ x[2:end]) ./ 2

# Values of the axis at which the quantity is given: the samples, or the bin centers
_points(c::Container, axis::Symbol) = (x = getfield(c, axis); can_integrate(c, axis) ? x : bin_centers(x))

Base.getindex(c::Container, key) = _select_flavors(c, flavor_index(c.scheme, key))
_select_flavors(c::Container, i::Integer) = _select_flavors(c, [i])
function _select_flavors(c::Container, idx::AbstractVector{<:Integer})
    _is_full_flavor(c) || throw(SupernovaNeutrinosError("flavors can be selected only in a container with every flavor of $(c.scheme.name)"))
    return _with(c; array = c.array[idx, :, :], flavor = c.flavor[idx])
end

_is_full_flavor(c::Container) = c.flavor == collect(c.scheme) && size(c, 1) == length(c.flavor)

"""
    sum(c::Container, axis) -> Container

Sum over `:flavor` or over a binned `:time` or `:energy` axis. The summed axis
keeps one entry and only its first and last values.
"""
function Base.sum(c::Container, axis::Union{Symbol, AbstractString})
    a = _axis(axis)
    can_sum(c, Symbol(axis)) || throw(SupernovaNeutrinosError("cannot sum over the sampled $axis axis; integrate it"))
    array = sum(c.array; dims = a)
    if a == 1
        return _with(c; array, flavor = c.flavor[[1, end]])
    end
    x = getfield(c, Symbol(axis))[[1, end]]
    integrable = a == 2 ? (false, c.integrable[2]) : (c.integrable[1], false)
    return a == 2 ? _with(c; array, time = x, integrable) : _with(c; array, energy = x, integrable)
end

"""
    integrate(c::Container, axis; limits = nothing) -> Container

Integrate over the sampled `:time` or `:energy` axis, the quantity being
linear between the samples and zero outside them. `limits` are bin edges: the
result holds the integral in each bin. Without `limits`, integrate over the
whole axis.

```julia
fluence = integrate(F, :time)                                  # one time bin
binned = integrate(F, :time; limits = range(0, 1, length = 21)) # 20 bins of 50 ms
```
"""
function integrate(c::Container, axis::Union{Symbol, AbstractString}; limits = nothing)
    axis = Symbol(axis)
    a = _axis(axis)
    can_integrate(c, axis) || throw(SupernovaNeutrinosError("cannot integrate over the $axis axis; it is binned or the flavor axis"))
    x = getfield(c, axis)
    length(x) == 1 && return c
    limits = isnothing(limits) ? [first(x), last(x)] : collect(float.(limits))
    issorted(limits) || throw(SupernovaNeutrinosError("integration limits must be sorted"))
    issorted(x) || throw(SupernovaNeutrinosError("the $axis axis must be sorted to integrate"))
    A = c.array
    shape = collect(size(A))
    shape[a] = length(limits) - 1
    result = zeros(shape...)
    for I in CartesianIndices(size(A)[[i for i in 1:3 if i != a]])
        y = a == 2 ? view(A, I[1], :, I[2]) : view(A, I[1], I[2], :)
        F = [_integral_up_to(x, y, L) for L in limits]
        r = a == 2 ? view(result, I[1], :, I[2]) : view(result, I[1], I[2], :)
        r .= diff(F)
    end
    unit = c.unit * (a == 2 ? ContainerUnit(0, 1, 0) : ContainerUnit(1, 0, 0))
    integrable = a == 2 ? (false, c.integrable[2]) : (c.integrable[1], false)
    return a == 2 ? _with(c; array = result, time = limits, integrable, unit) :
                    _with(c; array = result, energy = limits, integrable, unit)
end

# ∫ from x[1] to L of the piecewise-linear function through (x, y), zero outside x
function _integral_up_to(x::AbstractVector, y::AbstractVector, L::Real)
    L <= x[1] && return 0.0
    total = 0.0
    for k in 1:length(x)-1
        x0, x1 = x[k], x[k+1]
        if L >= x1
            total += (x1 - x0) * (y[k] + y[k+1]) / 2
        else
            yL = y[k] + (y[k+1] - y[k]) * (L - x0) / (x1 - x0)
            return total + (L - x0) * (y[k] + yL) / 2
        end
    end
    return total
end

"Integrate over `axis` if it is sampled, sum over it if it is binned."
integrate_or_sum(c::Container, axis) = can_integrate(c, axis) ? integrate(c, axis) : sum(c, axis)

Base.:*(c::Container, x::Real) = _with(c; array = c.array .* x)
Base.:*(x::Real, c::Container) = c * x

# multiply by a function of energy given at the energy points of c
function _scale_energy(c::Container, v::AbstractVector, unit::ContainerUnit)
    length(v) == size(c, 3) || throw(DimensionMismatch("$(length(v)) values for $(size(c, 3)) energies"))
    return _with(c; array = c.array .* reshape(v, 1, 1, :), unit = c.unit * unit)
end

function Base.:*(M::FlavorMatrix, c::Container)
    _is_full_flavor(c) || throw(SupernovaNeutrinosError("a flavor matrix applies to a container with every flavor of $(c.scheme.name)"))
    M.flavor_in == c.scheme || throw(SupernovaNeutrinosError(
        "cannot apply $(M.flavor_out.name) ← $(M.flavor_in.name) to a container of $(c.scheme.name)"))
    nf, nt, ne = size(c)
    mt, me = size(M, 3), size(M, 4)
    (mt in (1, nt) && me in (1, ne)) || throw(DimensionMismatch(
        "flavor matrix with (time, energy) size ($mt, $me) does not fit a container of size ($nt, $ne)"))
    out = zeros(length(M.flavor_out), nt, ne)
    for e in 1:ne, t in 1:nt
        m = view(M.array, :, :, _at(mt, t), _at(me, e))
        out[:, t, e] = real(m * view(c.array, :, t, e))
    end
    return _with(c; array = out, flavor = collect(M.flavor_out), scheme = M.flavor_out)
end

"""
    convert_flavor(c::Container, scheme::FlavorScheme) -> Container

The container in another flavor scheme ([`conversion_matrix`](@ref)); also
written `scheme * c`.
"""
convert_flavor(c::Container, scheme::FlavorScheme) = c.scheme == scheme ? c : conversion_matrix(scheme, c.scheme) * c
Base.:*(s::FlavorScheme, c::Container) = convert_flavor(c, s)

Base.:(==)(a::Container, b::Container) =
    a.array == b.array && a.flavor == b.flavor && a.scheme == b.scheme && a.time == b.time &&
    a.energy == b.energy && a.integrable == b.integrable && a.unit == b.unit
Base.isapprox(a::Container, b::Container; kwargs...) =
    size(a) == size(b) && a.flavor == b.flavor && a.scheme == b.scheme && a.unit == b.unit &&
    isapprox(a.time, b.time) && isapprox(a.energy, b.energy) && isapprox(a.array, b.array; kwargs...)

_range(x) = isempty(x) ? "" : @sprintf("%g … %g", first(x), last(x))

function Base.show(io::IO, c::Container)
    print(io, "Container(", get(_KINDS, c.unit, "quantity"), ", ", join(size(c), "×"), ")")
end

function Base.show(io::IO, ::MIME"text/plain", c::Container)
    kind = get(_KINDS, c.unit, "quantity")
    println(io, "Container: ", kind, " [", c.unit, "], size ", join(size(c), "×"))
    println(io, "  flavor: ", join(string.(c.flavor), ", "), " (", c.scheme.name, ")")
    println(io, "  time:   ", length(c.time), c.integrable[1] ? " points " : " bin edges ", _range(c.time), " s")
    print(io, "  energy: ", length(c.energy), c.integrable[2] ? " points " : " bin edges ", _range(c.energy), " MeV")
end
