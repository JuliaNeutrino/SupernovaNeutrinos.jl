"""
    SupernovaNeutrinos

Neutrino fluxes from supernova simulations and the flavor transformations they
undergo between the star and the detector. A Julia port of SNEWPY
(https://github.com/SNEWS2/snewpy).

Quantities are plain numbers in fixed units: energies in MeV, times in
seconds, distances in kpc, areas in cm², luminosities in erg/s, masses in M☉,
neutrino mixing angles in degrees and mass-squared differences in eV².
"""
module SupernovaNeutrinos

include("constants.jl")
include("flavors.jl")

# flavors
export FlavorScheme, Flavor, TwoFlavor, ThreeFlavor, FourFlavor
export ThreeFlavorMassBasis, ThreeFlavorFlavorBasis, FourFlavorMassBasis, FourFlavorFlavorBasis
export is_neutrino, is_antineutrino, is_electron, is_muon, is_tauon, is_sterile, lepton, to_tex

end
