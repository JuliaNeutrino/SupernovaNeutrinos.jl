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

using LinearAlgebra
using Printf
using OrderedCollections: OrderedDict
using SpecialFunctions: loggamma

include("constants.jl")
include("flavors.jl")
include("flavor_matrix.jl")
include("mixing.jl")
include("container.jl")
include("transformations.jl")
include("models/base.jl")
include("models/extended.jl")

# flavors
export FlavorScheme, Flavor, TwoFlavor, ThreeFlavor, FourFlavor
export ThreeFlavorMassBasis, ThreeFlavorFlavorBasis, FourFlavorMassBasis, FourFlavorFlavorBasis
export is_neutrino, is_antineutrino, is_electron, is_muon, is_tauon, is_sterile, lepton, to_tex
export FlavorMatrix, conversion_matrix, convert_flavor

# mixing
export MassOrdering, NormalOrdering, InvertedOrdering
export ThreeFlavorMixingParameters, FourFlavorMixingParameters, MixingParameters
export mixing_angles, mass_squared_differences, vacuum_mixing_matrix, pmf_high_density_limit
export mass_basis, flavor_basis

# flux containers
export Container, integrate, integrate_or_sum, can_integrate, can_sum, bin_centers

# flavor transformations
export FlavorTransformation, SNTransformation, VacuumTransformation, EarthTransformation
export P_ff, P_mf, P_mm, P_fm, apply
export NoTransformation, CompleteExchange, TransformationChain
export AdiabaticMSW, NonAdiabaticMSWH, TwoFlavorDecoherence, ThreeFlavorDecoherence
export AdiabaticMSWes, NonAdiabaticMSWes, FourFlavorDecoherence
export NoVacuumTransformation, NeutrinoDecay, QuantumDecoherence
export NoEarthMatter

# supernova models
export SupernovaModel, PinchedModel
export initial_spectra, transformed_spectra, flux, times, metadata
export ExtendedModel, extended_time_dependence, is_extended_tail

end
