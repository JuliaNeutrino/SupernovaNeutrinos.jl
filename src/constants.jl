# Physical constants and unit conversions (CODATA 2018 and IAU 2015 values).
# Every quantity in the package is a plain number in the units listed in the
# module docstring; these convert from other units to those.

"One MeV in erg."
const MEV_IN_ERG = 1.602176634e-6

"One kiloparsec in cm."
const KPC_IN_CM = 3.0856775814913673e21

"Speed of light in cm/s."
const SPEED_OF_LIGHT = 2.99792458e10

"ħc in eV·cm."
const HBARC_EV_CM = 1.973269804e-5

"One hour in seconds."
const HOUR = 3600.0

"One day in seconds."
const DAY = 86400.0

"Error thrown for invalid input to this package."
struct SupernovaNeutrinosError <: Exception
    msg::String
end
Base.showerror(io::IO, e::SupernovaNeutrinosError) = print(io, "SupernovaNeutrinosError: ", e.msg)
