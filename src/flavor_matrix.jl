# Matrices between flavor schemes, such as transition probabilities
# P(να → νβ) or the mixing matrix from the mass to the flavor basis.

"""
    FlavorMatrix(array, flavor_out, flavor_in = flavor_out)
    FlavorMatrix(f, flavor_out, flavor_in = flavor_out)

A matrix from the states of `flavor_in` to the states of `flavor_out`:
`M[β, α]` is the element from state `α` to state `β`. The array has four
dimensions `(flavor_out, flavor_in, time, energy)`; the time and energy
dimensions have length 1 when the matrix does not depend on them. A
two-dimensional `array` gets length-1 time and energy dimensions.

The second form fills the matrix with `f(β, α)` for each pair of states.

Index with any key of a [`FlavorScheme`](@ref):

```julia
P = P_ff(AdiabaticMSW(), [0.0], [10.0])
P[:e, :e]            # survival probability of νe
P[:NU, :NU]          # the neutrino block
```

Matrices multiply with `*` (over the flavor dimensions, elementwise over time
and energy). A flavor scheme on either side converts the matrix to that
scheme ([`conversion_matrix`](@ref)): `TwoFlavor * P * TwoFlavor`.
"""
struct FlavorMatrix{T <: Number}
    array::Array{T, 4}
    flavor_out::FlavorScheme
    flavor_in::FlavorScheme
    function FlavorMatrix(array::AbstractArray{T}, flavor_out::FlavorScheme, flavor_in::FlavorScheme = flavor_out) where {T <: Number}
        A = if ndims(array) == 2
            reshape(Array(array), size(array)..., 1, 1)
        elseif ndims(array) == 4
            Array(array)
        else
            throw(SupernovaNeutrinosError(
                "a FlavorMatrix array has 2 or 4 dimensions (flavor_out, flavor_in[, time, energy]), not $(ndims(array))"))
        end
        expected = (length(flavor_out), length(flavor_in))
        size(A)[1:2] == expected || throw(SupernovaNeutrinosError(
            "FlavorMatrix array of size $(size(A)) does not match the $(expected) states of $(flavor_out.name) ← $(flavor_in.name)"))
        return new{T}(A, flavor_out, flavor_in)
    end
end

function FlavorMatrix(f::Function, flavor_out::FlavorScheme, flavor_in::FlavorScheme = flavor_out)
    return FlavorMatrix(Float64[f(β, α) for β in flavor_out, α in flavor_in], flavor_out, flavor_in)
end

"Identity matrix on a scheme, with optional time and energy dimensions."
FlavorMatrix(::UniformScaling, s::FlavorScheme; nt::Int = 1, ne::Int = 1) =
    FlavorMatrix(repeat(Matrix{Float64}(I, length(s), length(s)), 1, 1, nt, ne), s, s)

"Zero matrix from `flavor_in` to `flavor_out`, with optional time and energy dimensions."
zero_matrix(flavor_out::FlavorScheme, flavor_in::FlavorScheme = flavor_out; nt::Int = 1, ne::Int = 1) =
    FlavorMatrix(zeros(length(flavor_out), length(flavor_in), nt, ne), flavor_out, flavor_in)

Base.size(M::FlavorMatrix) = size(M.array)
Base.size(M::FlavorMatrix, d) = size(M.array, d)
Base.eltype(::FlavorMatrix{T}) where {T} = T
Base.copy(M::FlavorMatrix) = FlavorMatrix(copy(M.array), M.flavor_out, M.flavor_in)

_indices(M::FlavorMatrix, out, in) = (flavor_index(M.flavor_out, out), flavor_index(M.flavor_in, in))

# M[β, α]: a number if the matrix does not depend on time and energy, otherwise
# the (time, energy) matrix; vector or colon keys keep the flavor dimensions
function Base.getindex(M::FlavorMatrix, out, in)
    i, j = _indices(M, out, in)
    # without time and energy dimensions, a number or a plain matrix
    size(M, 3) == 1 && size(M, 4) == 1 && return M.array[i, j, 1, 1]
    return M.array[i, j, :, :]
end

function Base.setindex!(M::FlavorMatrix, value, out, in)
    i, j = _indices(M, out, in)
    M.array[i, j, :, :] .= value
    return M
end

Base.:(==)(A::FlavorMatrix, B::FlavorMatrix) =
    A.flavor_out == B.flavor_out && A.flavor_in == B.flavor_in && A.array == B.array
Base.isapprox(A::FlavorMatrix, B::FlavorMatrix; kwargs...) =
    A.flavor_out == B.flavor_out && A.flavor_in == B.flavor_in && size(A) == size(B) &&
    isapprox(A.array, B.array; kwargs...)

for f in (:real, :imag, :abs, :abs2, :conj)
    @eval Base.$f(M::FlavorMatrix) = FlavorMatrix($f.(M.array), M.flavor_out, M.flavor_in)
end
Base.transpose(M::FlavorMatrix) = FlavorMatrix(permutedims(M.array, (2, 1, 3, 4)), M.flavor_in, M.flavor_out)
Base.adjoint(M::FlavorMatrix) = conj(transpose(M))
Base.:*(x::Number, M::FlavorMatrix) = FlavorMatrix(x .* M.array, M.flavor_out, M.flavor_in)
Base.:*(M::FlavorMatrix, x::Number) = x * M
Base.:+(A::FlavorMatrix, B::FlavorMatrix) = FlavorMatrix(_combine(+, A, B), A.flavor_out, A.flavor_in)
Base.:-(A::FlavorMatrix, B::FlavorMatrix) = FlavorMatrix(_combine(-, A, B), A.flavor_out, A.flavor_in)

function _combine(op, A::FlavorMatrix, B::FlavorMatrix)
    (A.flavor_out == B.flavor_out && A.flavor_in == B.flavor_in) || throw(SupernovaNeutrinosError(
        "cannot combine matrices $(A.flavor_out.name) ← $(A.flavor_in.name) and $(B.flavor_out.name) ← $(B.flavor_in.name)"))
    return op.(A.array, B.array)
end

# Length of a time or energy dimension after broadcasting lengths a and b
function _broadcast_length(a::Int, b::Int)
    a == b && return a
    a == 1 && return b
    b == 1 && return a
    throw(DimensionMismatch("time or energy dimensions of length $a and $b do not broadcast"))
end
_at(n::Int, i::Int) = n == 1 ? 1 : i

function Base.:*(A::FlavorMatrix, B::FlavorMatrix)
    A.flavor_in == B.flavor_out || throw(SupernovaNeutrinosError(
        "cannot multiply $(A.flavor_out.name) ← $(A.flavor_in.name) by $(B.flavor_out.name) ← $(B.flavor_in.name)"))
    nt = _broadcast_length(size(A, 3), size(B, 3))
    ne = _broadcast_length(size(A, 4), size(B, 4))
    C = zeros(promote_type(eltype(A), eltype(B)), size(A, 1), size(B, 2), nt, ne)
    for e in 1:ne, t in 1:nt
        a = view(A.array, :, :, _at(size(A, 3), t), _at(size(A, 4), e))
        b = view(B.array, :, :, _at(size(B, 3), t), _at(size(B, 4), e))
        mul!(view(C, :, :, t, e), a, b)
    end
    return FlavorMatrix(C, A.flavor_out, B.flavor_in)
end

"""
    conversion_matrix(to::FlavorScheme, from::FlavorScheme) -> FlavorMatrix

Matrix converting quantities per state of `from` to quantities per state of
`to`. States with the same name map to each other; `X` of [`TwoFlavor`](@ref)
becomes μ and τ each (factor 1), and μ and τ become their average `X`
(factor 1/2).
"""
function conversion_matrix(to::FlavorScheme, from::FlavorScheme)
    M = FlavorMatrix(to, from) do β, α
        Symbol(β) == Symbol(α) && return 1.0
        is_neutrino(β) == is_neutrino(α) || return 0.0
        lepton(α) == "X" && lepton(β) in ("MU", "TAU") && return 1.0
        lepton(β) == "X" && lepton(α) in ("MU", "TAU") && return 0.5
        return 0.0
    end
    all(iszero, M.array) && throw(SupernovaNeutrinosError("no conversion from $(from.name) to $(to.name)"))
    return M
end

"""
    convert_flavor(M::FlavorMatrix; out = M.flavor_out, in = M.flavor_in)

`M` with its output states converted to `out` and its input states to `in`.
`TwoFlavor * M * TwoFlavor` is a shorthand.
"""
function convert_flavor(M::FlavorMatrix; out::FlavorScheme = M.flavor_out, in::FlavorScheme = M.flavor_in)
    in == M.flavor_in || (M = M * conversion_matrix(M.flavor_in, in))
    out == M.flavor_out || (M = conversion_matrix(out, M.flavor_out) * M)
    return M
end

Base.:*(s::FlavorScheme, M::FlavorMatrix) = convert_flavor(M; out = s)
Base.:*(M::FlavorMatrix, s::FlavorScheme) = convert_flavor(M; in = s)

function Base.show(io::IO, M::FlavorMatrix)
    print(io, "FlavorMatrix(", M.flavor_out.name, " ← ", M.flavor_in.name, ", size ", size(M), ")")
end

function Base.show(io::IO, ::MIME"text/plain", M::FlavorMatrix)
    show(io, M)
    size(M, 3) == size(M, 4) == 1 || return
    println(io)
    width = maximum(length ∘ string, M.flavor_out.states) + 1
    print(io, rpad("", width))
    for α in M.flavor_in
        print(io, lpad(string(α), 13))
    end
    for (i, β) in enumerate(M.flavor_out)
        print(io, "\n", rpad(string(β), width))
        for j in eachindex(M.flavor_in)
            x = M.array[i, j, 1, 1]
            print(io, lpad(x isa Real ? @sprintf("%.4g", x) : @sprintf("%.3g%+.3gi", real(x), imag(x)), 13))
        end
    end
end
