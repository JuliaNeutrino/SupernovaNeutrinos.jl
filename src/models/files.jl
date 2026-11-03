# Readers for the text, tar and zip files of the model repositories.

_parse_floats(line::AbstractString) = [parse(Float64, replace(w, r"[dD]" => "e")) for w in split(line)]
_is_numeric(line::AbstractString) = (w = split(line); !isempty(w) && all(x -> !isnothing(tryparse(Float64, x)), w))

"""
    read_table(source; skip = 0, comment = "#") -> Matrix{Float64}

Whitespace-separated numeric table of a file path or `IO`: one row per line,
skipping the first `skip` lines, blank lines and lines starting with `comment`.
Rows shorter than the first row are dropped (truncated last lines).
"""
read_table(path::AbstractString; kwargs...) = open(io -> read_table(io; kwargs...), path)

function read_table(io::IO; skip::Int = 0, comment::AbstractString = "#")
    rows = Vector{Float64}[]
    for (n, line) in enumerate(eachline(io))
        n <= skip && continue
        s = strip(line)
        (isempty(s) || startswith(s, comment)) && continue
        push!(rows, _parse_floats(s))
    end
    isempty(rows) && throw(SupernovaNeutrinosError("no data in the table"))
    ncol = length(first(rows))
    rows = filter(r -> length(r) == ncol, rows)
    return permutedims(reduce(hcat, rows))
end

"""
    read_named_table(source; comment = "#") -> OrderedDict{String, Vector{Float64}}

Numeric table whose column names are in the last line before the data: a
comment line (`# TIME L_NU_E …`) or a plain line of names.
"""
read_named_table(path::AbstractString; kwargs...) = open(io -> read_named_table(io; kwargs...), path)

function read_named_table(io::IO; comment::AbstractString = "#")
    header = String[]
    rows = Vector{Float64}[]
    for line in eachline(io)
        s = strip(line)
        isempty(s) && continue
        if !isempty(rows) || _is_numeric(s)
            startswith(s, comment) || push!(rows, _parse_floats(s))
        else
            header = split(lstrip(s, collect(comment)...))
        end
    end
    isempty(rows) && throw(SupernovaNeutrinosError("no data in the table"))
    length(header) == length(first(rows)) || throw(SupernovaNeutrinosError(
        "$(length(header)) column names for $(length(first(rows))) columns"))
    A = permutedims(reduce(hcat, filter(r -> length(r) == length(header), rows)))
    return OrderedDict(String(h) => A[:, i] for (i, h) in enumerate(header))
end

"Contents of the member `name` of a `.tar.gz` archive."
function read_tar_member(path::AbstractString, name::AbstractString)
    return mktempdir() do dir
        open(path) do io
            Tar.extract(h -> h.path == name || h.path == "./" * name, CodecZlib.GzipDecompressorStream(io), dir)
        end
        file = isfile(joinpath(dir, name)) ? joinpath(dir, name) : joinpath(dir, ".", name)
        isfile(file) || throw(SupernovaNeutrinosError("no file $name in $path"))
        read(file, String)
    end
end

"Contents of every member of a zip archive whose name satisfies `keep`, as `name => text`."
function read_zip_members(keep, path::AbstractString)
    r = ZipFile.Reader(path)
    try
        return [f.name => read(f, String) for f in r.files if keep(f.name)]
    finally
        close(r)
    end
end

# Times and columns of a pinched model in the order TIME, L, E, RMS of the
# flavors NU_E, NU_E_BAR, NU_X: the format of O'Connor, Bugli…
function _pinched_from_rms(A::AbstractMatrix; L_scale = 1.0, rms_is_e2 = false)
    cols = Dict{String, Vector{Float64}}("TIME" => A[:, 1])
    for (k, f) in enumerate(("NU_E", "NU_E_BAR", "NU_X"))
        L, E, R = A[:, 1+k] .* L_scale, A[:, 4+k], A[:, 7+k]
        cols["L_$f"] = L
        cols["E_$f"] = E
        cols["ALPHA_$f"] = rms_is_e2 ? alpha_from_e2(E, R) : alpha_from_rms(E, R)
    end
    return cols
end
