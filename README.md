# SupernovaNeutrinos.jl

Neutrino fluxes from core-collapse and pre-supernova simulations and the
flavor transformations they undergo on the way to Earth. A Julia port of
[SNEWPY](https://github.com/SNEWS2/snewpy). Event rates in neutrino detectors
come from [SNOwGLoBES.jl](https://github.com/JuliaNeutrino/SNOwGLoBES.jl),
which takes the fluences of this package directly.

```julia
using SupernovaNeutrinos

model = Bollig_2016(progenitor_mass = 27)
E = range(0, 100, length = 201)                     # MeV
F = flux(model, times(model), E; distance = 10,      # kpc
         transformation = AdiabaticMSW())
fluence = integrate(F, :time)                       # 1/(MeV cm²)

using SNOwGLoBES
rates = event_rates(SupernovaDetector(:icecube, :water), fluence)
total_events(rates; interaction = :ibd)
```

Julia 1.11 or later. `SNOwGLoBES.Fluence(fluence)` converts a fluence; a
fluence with several time bins (`integrate(F, :time; limits = ...)`) gives one
set of rates per bin.

The notebook [`notebooks/01_supernova_neutrinos.ipynb`](notebooks/01_supernova_neutrinos.ipynb) goes through the package: the models, luminosity curves, flavor transformations, event rates in DUNE, Hyper-K and IceCube, the neutronization burst and pre-supernova neutrinos ([open in Colab](https://colab.research.google.com/github/JuliaNeutrino/SupernovaNeutrinos.jl/blob/main/notebooks/01_supernova_neutrinos.ipynb); choose **Runtime → Change runtime type → Julia** first).

## Units

Quantities are plain `Float64` in fixed units: energies in MeV, times in s,
distances in kpc, areas in cm², luminosities in erg/s, masses in M☉, mixing
angles in degrees and mass-squared differences in eV². A `Container` (the
`(flavor, time, energy)` array of a spectrum, flux, fluence or rate) records
its powers of MeV, s and cm.

## Models

`all_models()` lists the 24 models of the SNEWS model repositories
(`Bollig_2016`, `Nakazato_2013`, `Fornax_2022`, `Warren_2020`, …, and the
pre-supernova models `Patton_2017`, `Kato_2017`, `Odrzywolek_2010`,
`Yoshida_2016`, `Myers_2026`). A model is selected by its parameters:

```julia
parameters(Nakazato_2013)            # allowed values of each parameter
param_combinations(Nakazato_2013)    # the simulations that exist
m = Nakazato_2013(progenitor_mass = 13, revival_time = 100, metallicity = 0.02, eos = "shen")
```

Model files are downloaded on first use to a scratch space, or to the
directory in `SUPERNOVANEUTRINOS_MODELS`. Models stored as FITS or HDF5 files
need `using FITSIO` or `using HDF5`. `PinchedModel`, `Analytic3Species`,
`ExtendedModel` (a cooling tail after the end of a simulation) and the
`SupernovaModel` interface let you add your own.

## Flavor transformations

In the supernova: `AdiabaticMSW`, `NonAdiabaticMSWH`, `TwoFlavorDecoherence`,
`ThreeFlavorDecoherence`, `AdiabaticMSWes`, `NonAdiabaticMSWes`,
`FourFlavorDecoherence`, `CompleteExchange`, `NoTransformation`. In vacuum:
`NeutrinoDecay`, `QuantumDecoherence`. They combine with
`TransformationChain(sn_step, vacuum_step; mixing_params = MixingParameters(:inverted))`;
the mixing presets are NuFIT 5.0, 5.2, 6.0 (default) and PDG 2022, 2024.

## Differences from SNEWPY

These are bugs in SNEWPY that are fixed here:

1. The μ–τ content of the matter eigenstates at high density used
   `D = √(|b|² − (a−c)²/4)` instead of `√(|b|² + (a−c)²/4)`; it now agrees
   with the eigenvectors of the Hamiltonian.
2. `Bugli_2021` squared the ⟨E²⟩ column of the files a second time, which
   gave pinch parameters α ≈ −1 instead of α ≈ 3–5.
3. Spectra are zero outside the simulated times. SNEWPY held the first and
   last values for all earlier and later times.
4. Pre-supernova tables are zero outside their time range (SNEWPY's
   log-interpolation gave about 1 per MeV·s there).
5. `ExtendedModel` starts its tail from the spectrum at the end of the
   simulation, not at the last requested time.
6. `Fornax_2019` evaluates its spherical-harmonic expansion exactly in the
   requested direction rather than at the center of a HEALPix pixel. The
   direction is given when the model is built.
7. The Fornax models interpolate linearly in time between snapshots.
8. Transformations of binned axes are evaluated at the bin centers.
9. Energy-dependent vacuum transformations (`NeutrinoDecay`,
   `QuantumDecoherence`) work in chains with energy-dependent supernova
   steps, where SNEWPY's broadcasting failed.
10. Simulations whose time column steps back on a restart (`Warren_2020`)
    keep the rows of the restart.
11. Integration over time or energy is the exact integral of the linear
    interpolation, including partial bins, and nothing is modified in place.

The event totals for the `Bollig_2016` 11.2 M☉ model agree with SNEWPY's
cross-check values to 1 % (see `test/network.jl`).

Not ported: the SNOSHEWS and BEMEWS (Earth matter) transformations, plotting,
and saving fluxes to files. SNEWPY's rate calculator is replaced by
SNOwGLoBES.jl.

## Tests

```julia
using Pkg; Pkg.test("SupernovaNeutrinos")
```

The tests download one model and, with SNOwGLoBES.jl, the SNOwGLoBES data; set
`SUPERNOVANEUTRINOS_TEST_NETWORK=false` to skip them.

## License

BSD 3-Clause, like SNEWPY. SNOwGLoBES.jl (GPL-2.0) is an optional dependency,
loaded only when you load it. The model files keep their own licenses.
