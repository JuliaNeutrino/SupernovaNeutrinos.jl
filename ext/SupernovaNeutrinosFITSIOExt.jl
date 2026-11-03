module SupernovaNeutrinosFITSIOExt

using FITSIO
import SupernovaNeutrinos

# The columns of the first binary table of a FITS file, by name
function SupernovaNeutrinos.read_fits_table(path::AbstractString)
    return FITS(path, "r") do f
        hdu = first(h for h in f if h isa TableHDU)
        Dict{String, Vector{Float64}}(name => Float64.(read(hdu, name)) for name in FITSIO.colnames(hdu))
    end
end

end
