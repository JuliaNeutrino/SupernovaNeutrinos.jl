module SupernovaNeutrinosHDF5Ext

using HDF5
import SupernovaNeutrinos

# Every dataset of an HDF5 file in nested dictionaries, the attributes of
# each group or dataset under "@attrs" (arrays transposed by HDF5.jl:
# column-major order of a row-major file)
function SupernovaNeutrinos.read_hdf5(path::AbstractString)
    return h5open(_read, path, "r")
end

_attrs(obj) = Dict{String, Any}(k => read(attributes(obj)[k]) for k in keys(attributes(obj)))

function _read(g::Union{HDF5.File, HDF5.Group})
    d = Dict{String, Any}("@attrs" => _attrs(g))
    for k in keys(g)
        obj = g[k]
        if obj isa HDF5.Dataset
            d[k] = read(obj)
            a = _attrs(obj)
            isempty(a) || (d[k*"@attrs"] = a)
        else
            d[k] = _read(obj)
        end
    end
    return d
end

end
