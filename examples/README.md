# PROPOSAL.jl Examples

Julia equivalents of the [PROPOSAL Python examples](https://github.com/tudo-astroparticlephysics/PROPOSAL/tree/master/examples).

## Examples (matching Python notebooks)

| # | Julia File | Python Equivalent | Concepts |
|---|-----------|-------------------|----------|
| 1 | `01_energy_cuts.jl` | `EnergyCut.ipynb` | EnergyCutSettings, ecut/vcut separation |
| 2 | `02_crosssections.jl` | `CrossSection.ipynb` | Standard cross sections, dEdx, dNdx, stochastic loss sampling |
| 3 | `03_parametrizations.jl` | `Parametrization.ipynb` | Differential cross sections, kinematic limits, comparing parametrizations |
| 4 | `04_displacement.jl` | `Displacement.ipynb` | Track integral, energy loss vs distance |
| 5 | `05_interaction.jl` | `Interaction.ipynb` | Stochastic interaction sampling, mean free path, rates |
| 6 | `06_decay.jl` | `Decay.ipynb` | Decay energy sampling |
| 7 | `07_time.jl` | `Time.ipynb` | Exact vs approximate time calculation |
| 8 | `08_continuous_randomization.jl` | `ContRand.ipynb` | Energy fluctuations in continuous losses |
| 9 | `09_secondaries.jl` | `Secondaries.ipynb` | SecondariesCalculator, secondary particle production |
| 10 | `10_advanced_propagator.jl` | `AdvancedPropagator.ipynb` | Building propagator from components, track inspection |

## Additional examples (Julia-only features)

| # | Julia File | Concepts |
|---|-----------|----------|
| 11 | `11_density_distributions.jl` | Homogeneous, exponential, polynomial, spline densities |
| 12 | `12_math_types.jl` | Polynom, LinearSpline, CubicSpline |
| 13 | `13_particle_lookup.jl` | ParticleDef, Builder, type lookup, constants |
| 14 | `14_geometry.jl` | Sphere, Cylinder, Box queries, hit detection, entry/exit points |
| 15 | `15_scattering.jl` | Highland, Moliere, stochastic deflection, ScatteringMultiplier |
| 16 | `16_multi_sector_propagator.jl` | Multi-sector propagation with different media/densities |
| 17 | `17_decay_channels.jl` | LeptonicDecayChannel, TwoBodyPhaseSpace, ManyBodyPhaseSpace, DecayTable |
| 18 | `18_shadow_effects.jl` | ShadowDuttaRenoSarcevicSeckel, ShadowButkevichMikheyev, PhotoQ2 with shadow |
| 19 | `19_lpm_effects.jl` | BremsLPM, EpairLPM, PhotoPairLPM suppression factors (detailed) |
| 20 | `20_medium_components.jl` | Medium properties, Component details, hash lookups, Sternheimer parameters |
| 21 | `21_config_propagator.jl` | JSON configuration-based propagator, per-particle factories |
| 22 | `22_settings.jl` | InterpolationSettings, PropagationSettings, interpolated vs exact comparison |
