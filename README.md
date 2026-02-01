# PROPOSAL.jl

Julia bindings for [PROPOSAL](https://github.com/tudo-astroparticlephysics/PROPOSAL) (Propagator with Optimal Precision and Optimized Speed for All Leptons).

PROPOSAL is a library for propagating leptons and gamma rays through media, featuring up-to-date cross sections for ionization, bremsstrahlung, photonuclear interactions, electron pair production, Landau-Pomeranchuk-Migdal and Ter-Mikaelian effects, muon and tau decay, as well as Moliere scattering.

## Installation

Once registered in the Julia General registry:

```julia
using Pkg
Pkg.add("PROPOSAL")
```

### Development Install

To use PROPOSAL.jl before it is registered, clone the repository and install in development mode:

```julia
using Pkg
Pkg.develop(path="/path/to/PROPOSAL.jl")
```

Note that the native library must be available for the bindings to work. You can either:
1. Wait for `PROPOSAL_cxxwrap_jll` to be registered (provides the library automatically)
2. Build the wrapper library locally (see below)
3. Set the `PROPOSAL_JL_LIB_PATH` environment variable to a directory containing the built library

### Building the Wrapper Library

The wrapper library bridges PROPOSAL's C++ API to Julia via CxxWrap. Building it requires:

- C++17 compatible compiler
- CMake 3.15+
- [PROPOSAL](https://github.com/tudo-astroparticlephysics/PROPOSAL) (installed with CMake config files)
- [CxxWrap.jl](https://github.com/JuliaInterop/CxxWrap.jl) (`Pkg.add("CxxWrap")`)

#### 1. Install PROPOSAL C++ library

```bash
git clone https://github.com/tudo-astroparticlephysics/PROPOSAL.git
cd PROPOSAL
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$HOME/.local -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=Release
cmake --build . --parallel
cmake --install .
```

**Important:** Use `-DCMAKE_BUILD_TYPE=Release` to disable C++ assertions. Debug/default builds include assertions that can cause `SIGABRT` crashes during propagation (e.g., `energy_initial >= energy_final` in `PropagationUtilityInterpolant.cxx`).

**Known issue:** The installed `PROPOSALConfig.cmake` (at `$HOME/.local/lib/cmake/PROPOSAL/PROPOSALConfig.cmake`) may have broken paths. If the wrapper build fails to find PROPOSAL, edit the file:

1. Change the prefix calculation from `../../` to `../../../`:
   ```cmake
   get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)
   ```
2. Fix the legacy variable paths:
   ```cmake
   set_and_check(PROPOSAL_INCLUDE_DIR  "${PACKAGE_PREFIX_DIR}/include")
   set_and_check(PROPOSAL_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include")
   set_and_check(PROPOSAL_LIBRARIES    "${PACKAGE_PREFIX_DIR}/lib/libPROPOSAL.dylib")  # .so on Linux
   ```

#### 2. Get the JlCxx CMake path

```bash
JLCXX_DIR=$(julia -e 'using CxxWrap; print(CxxWrap.prefix_path())')/lib/cmake/JlCxx
```

#### 3. Build the wrapper

```bash
cd PROPOSAL.jl/deps/binarybuilder/wrapper_src
mkdir build && cd build
cmake .. \
    -DJlCxx_DIR=$JLCXX_DIR \
    -DPROPOSAL_DIR=$HOME/.local/lib/cmake/PROPOSAL \
    -DCMAKE_BUILD_TYPE=Release
cmake --build . --parallel
```

#### 4. Install the library

The Julia code expects the library to be named `libPROPOSAL_cxxwrap`, so copy and rename it:

```bash
mkdir -p PROPOSAL.jl/build/lib
cp build/libPROPOSAL_jl.dylib PROPOSAL.jl/build/lib/libPROPOSAL_cxxwrap.dylib  # macOS
# cp build/libPROPOSAL_jl.so PROPOSAL.jl/build/lib/libPROPOSAL_cxxwrap.so      # Linux
```

Or set the environment variable instead (the library must still be named `libPROPOSAL_cxxwrap`):

```bash
export PROPOSAL_JL_LIB_PATH=/path/to/directory/containing/libPROPOSAL_cxxwrap
```

After copying the library, you must clear the Julia precompile cache so the package picks up the new library:

```bash
rm -rf ~/.julia/compiled/v1.*/PROPOSAL
```

Then start Julia and verify:

```julia
using PROPOSAL
is_library_available()  # should return true
```

#### Apple Silicon Note

On Apple Silicon Macs, ensure all components use the same architecture. If you encounter `symbol(s) not found for architecture x86_64` errors, verify that Julia, CxxWrap, PROPOSAL, and all dependencies are built for arm64:

```bash
file /path/to/library.dylib  # should show "arm64"
```

See `deps/binarybuilder/README.md` for BinaryBuilder-based build instructions and Yggdrasil submission details.

## Requirements

- Julia 1.6 or higher
- CxxWrap.jl 0.14 or 0.15

## Quick Start

```julia
using PROPOSAL

set_random_seed(42)

# --- Option 1: JSON config-based propagator ---
propagator = create_propagator_muminus("path/to/config.json")

# --- Option 2: Build from components ---
pdef = MuMinusDef()
medium = create_ice()
handle = create_crosssections(pdef, medium, 500.0, 0.05, false, true)

coll = PropagationUtilityCollection()
set_displacement!(coll, make_displacement(handle, true))
set_interaction!(coll, make_interaction(handle, true))
set_time!(coll, make_time_calc(handle, pdef, true))

geo = Sphere(Cartesian3D(0, 0, 0), 1e20, 0.0)
dens = DensityHomogeneous(get_mass_density(medium))
propagator = make_propagator_single_sector(pdef, geo, dens, coll)

# --- Propagate ---
state = ParticleState(
    PARTICLE_TYPE_MUMINUS,
    0.0, 0.0, 0.0,      # position (x, y, z) in cm
    0.0, 0.0, 1.0,      # direction (dx, dy, dz)
    1e6;                 # energy in MeV (1 TeV)
)

secondaries = propagate(propagator, state; max_distance=1e7, min_energy=1e3)

final_state = get_final_state(secondaries)
println("Final energy: ", get_energy(final_state), " MeV")
println("Distance: ", get_propagated_distance(final_state), " cm")

# Inspect stochastic losses
for loss in get_stochastic_losses(secondaries)
    println("  type=", get_type(loss), " E=", get_energy(loss), " MeV")
end

free_crosssections(handle)
```

## Examples

The [`examples/`](examples/) directory contains 22 runnable Julia scripts covering all features:

| # | File | Concepts |
|---|------|----------|
| 1-10 | `01_energy_cuts.jl` .. `10_advanced_propagator.jl` | Julia equivalents of the [PROPOSAL Python examples](https://github.com/tudo-astroparticlephysics/PROPOSAL/tree/master/examples) |
| 11 | `11_density_distributions.jl` | Homogeneous, exponential, polynomial, spline densities |
| 12 | `12_math_types.jl` | Polynom, LinearSpline, CubicSpline |
| 13 | `13_particle_lookup.jl` | ParticleDef, Builder, type lookup, constants |
| 14 | `14_geometry.jl` | Sphere, Cylinder, Box queries, hit detection, entry/exit points |
| 15 | `15_scattering.jl` | Highland, Moliere, stochastic deflection, ScatteringMultiplier |
| 16 | `16_multi_sector_propagator.jl` | Multi-sector propagation with different media/densities |
| 17 | `17_decay_channels.jl` | LeptonicDecayChannel, TwoBodyPhaseSpace, ManyBodyPhaseSpace, DecayTable |
| 18 | `18_shadow_effects.jl` | Nuclear shadow effects on photonuclear cross sections |
| 19 | `19_lpm_effects.jl` | LPM suppression factors for bremsstrahlung, pair production, photo-pair |
| 20 | `20_medium_components.jl` | Medium properties, Component details, hash lookups |
| 21 | `21_config_propagator.jl` | JSON configuration-based propagator, per-particle factories |
| 22 | `22_settings.jl` | InterpolationSettings, PropagationSettings |

See [`examples/README.md`](examples/README.md) for the full table.

## Available Types and Functions

### Particle Definitions
- `MuMinusDef()`, `MuPlusDef()` — Muons
- `EMinusDef()`, `EPlusDef()` — Electrons/Positrons
- `TauMinusDef()`, `TauPlusDef()` — Taus
- `GammaDef()` — Photons
- `Pi0Def()`, `PiMinusDef()`, `PiPlusDef()` — Pions
- `K0Def()`, `KMinusDef()`, `KPlusDef()` — Kaons
- `NuEDef()`, `NuMuDef()`, `NuTauDef()` (and bar variants) — Neutrinos
- `ParticleDefBuilder` — Custom particles with `set_mass`, `set_lifetime`, `set_decay_table`, `build`

### Media
- `create_water()`, `create_ice()`, `create_air()`, `create_standard_rock()`, `create_iron_medium()`, `create_lead_medium()`, `create_liquid_argon()`, and more
- `get_mass_density`, `get_radiation_length`, `get_components`, Sternheimer parameters
- Named components: `ComponentHydrogen()`, `ComponentCarbon()`, `ComponentOxygen()`, `ComponentIron()`, etc.
- `get_component_for_hash(hash)` — lookup by hash

### Geometry
- `Sphere(center, outer_radius, inner_radius)`, `Cylinder(center, half_height, radius, inner_radius)`, `Box(center, half_x, half_y, half_z)`
- `is_inside`, `is_infront`, `is_behind`, `distance_to_border`
- `set_hierarchy` / `get_hierarchy` for multi-sector geometry

### Density Distributions
- `DensityHomogeneous(mass_density)`
- `DensityExponential(axis, sigma, d0, mass_density)`
- `DensityPolynomial(axis, polynom, mass_density)`
- `DensitySplines(axis, spline, mass_density)`
- Axes: `CartesianAxis(origin, direction)`, `RadialAxis(origin, direction)`

### Cross Sections
- `create_crosssections(pdef, medium, ecut, vcut, cont_rand, interpolate)` — standard cross sections
- `calculate_dEdx`, `calculate_dE2dx`, `calculate_dNdx`, `calculate_stochastic_loss`
- Individual parametrizations for bremsstrahlung, pair production, photonuclear, ionization, etc.

### Scattering
- `create_scattering_ms_only(name, pdef, medium)` — Highland, Moliere, MoliereInterpol
- `create_scattering_with_sd(ms_name, pdef, medium, sd_name)` — with stochastic deflection
- `create_scattering_multiplier(...)` — scale deflection angles

### Propagation
- `make_propagator_single_sector(pdef, geo, dens, coll)` — from components
- `clear_sectors` / `add_sector_*` / `make_propagator_from_sectors(pdef)` — multi-sector
- `create_propagator_muminus(config)`, `create_propagator(pdef, config)` — from JSON config
- `propagate(propagator, state; max_distance, min_energy)` — propagation
- `get_stochastic_losses(sec)`, `get_continuous_losses(sec)` — loss access
- `hit_geometry`, `get_entry_point`, `get_exit_point` — geometry intersection

### Decay
- `LeptonicDecayChannel`, `LeptonicDecayChannelApprox`, `TwoBodyPhaseSpace`, `ManyBodyPhaseSpace`
- `DecayTable` with `add_channel`, `select_channel`

### Shadow & LPM Effects
- `ShadowDuttaRenoSarcevicSeckel()`, `ShadowButkevichMikheyev()` — nuclear shadow
- `create_brems_lpm`, `create_epair_lpm`, `create_photopair_lpm` — LPM suppression
- `suppression_factor(lpm, ...)` — evaluate suppression

### Math Types
- `Cartesian3D(x, y, z)`, `Spherical3D(r, azimuth, zenith)`
- `Polynom(coefficients)`, `LinearSpline(x, y)`, `CubicSpline(x, y)`
- `evaluate`, `derive`, `antiderivative`

### Settings
- `set_tables_path(path)`, `set_upper_energy_lim(E)` — interpolation table config
- `set_nodes_dedx(n)`, `set_nodes_dndx_e(n)`, etc. — interpolation accuracy
- `set_max_steps(n)` — propagation step limit

### Constants
- `SPEED_OF_LIGHT`, `ELECTRON_MASS`, `MUON_MASS`, `TAU_MASS`, `PROTON_MASS`
- `PARTICLE_TYPE_MUMINUS`, `PARTICLE_TYPE_MUPLUS`, `PARTICLE_TYPE_EMINUS`, etc.
- `INTERACTION_TYPE_BREMS`, `INTERACTION_TYPE_EPAIR`, `INTERACTION_TYPE_PHOTONUCLEAR`, etc.

### Utility Functions
- `set_random_seed(seed)` — Set RNG seed for reproducibility
- `get_proposal_version()` — Get PROPOSAL version string
- `is_library_available()` — Check if the native library is loaded

## Configuration

PROPOSAL can also use JSON configuration files. Example minimal configuration:

```json
{
    "global": {
        "seed": 1234,
        "continous_loss_output": true
    },
    "sectors": [
        {
            "medium": "ice",
            "geometries": [
                {
                    "shape": "sphere",
                    "origin": [0, 0, 0],
                    "outer_radius": 1e20,
                    "inner_radius": 0
                }
            ],
            "cuts": {
                "e_cut": 500,
                "v_cut": 0.05,
                "cont_rand": false
            },
            "do_exact_time": true,
            "scattering": {
                "model": "highland"
            }
        }
    ]
}
```

For detailed configuration options, see the [PROPOSAL documentation](https://github.com/tudo-astroparticlephysics/PROPOSAL/blob/master/docs/config_docu.md).

## Troubleshooting

### `PROPOSALConfig.cmake` broken install paths

After installing PROPOSAL with `cmake --install .`, the generated `PROPOSALConfig.cmake` (located at `$CMAKE_INSTALL_PREFIX/lib/cmake/PROPOSAL/PROPOSALConfig.cmake`) may contain incorrect paths for `PROPOSAL_INCLUDE_DIR`, `PROPOSAL_INCLUDE_DIRS`, and `PROPOSAL_LIBRARIES`. For example:

```cmake
set_and_check(PROPOSAL_INCLUDE_DIR  "/include")  # missing prefix
```

Fix by editing the file to use the `PACKAGE_PREFIX_DIR` variable, and adjust the prefix calculation to point to the install root (three levels up from the cmake directory, not two):

```cmake
get_filename_component(PACKAGE_PREFIX_DIR "${CMAKE_CURRENT_LIST_DIR}/../../../" ABSOLUTE)
# ...
set_and_check(PROPOSAL_INCLUDE_DIR  "${PACKAGE_PREFIX_DIR}/include")
set_and_check(PROPOSAL_INCLUDE_DIRS "${PACKAGE_PREFIX_DIR}/include")
set_and_check(PROPOSAL_LIBRARIES    "${PACKAGE_PREFIX_DIR}/lib/libPROPOSAL.dylib")
```

### `StdFill` / STL wrapping errors with CxxWrap.jl 0.15

If you see an error like:

```
UndefVarError: `StdFill` not defined in `CxxWrap.StdLib`
```

This is caused by a mismatch between the CxxWrap C++ headers (which define `StdFill`) and the bundled `libcxxwrap_julia_jll` binary (v0.12.3, which does not export the required `StlWrappers::instance()` symbol). The wrapper's `CMakeLists.txt` must **not** link against `JlCxx::cxxwrap_julia_stl`, and the C++ wrapper source must avoid `#include "jlcxx/stl.hpp"` and not return `std::vector` types directly from wrapped functions. Instead, use indexed access patterns (size + element-at-index) to expose collection data.

### Precompile cache not updated after building the library

If `is_library_available()` returns `false` even though the library file exists at the expected path, the Julia precompile cache is stale. The library path is resolved at precompile time as a `const`. Clear the cache:

```bash
rm -rf ~/.julia/compiled/v1.*/PROPOSAL
```

Then restart Julia and `using PROPOSAL` again.

## License

This software is distributed under the GNU Lesser General Public License v3.0 (LGPL-3.0), the same license as the upstream PROPOSAL library.

### Citation Requirements

If you use PROPOSAL.jl, please cite the PROPOSAL paper:

> J.H. Koehne et al., "PROPOSAL: A tool for propagation of charged leptons"
> Comput. Phys. Commun. 184 (2013) 2070-2090
> DOI: [10.1016/j.cpc.2013.04.001](https://doi.org/10.1016/j.cpc.2013.04.001)

BibTeX:
```bibtex
@article{koehne2013proposal,
  title   = {PROPOSAL: A tool for propagation of charged leptons},
  author  = {Koehne, Jan-Hendrik and Frantzen, Katharina and Schmitz, Martin
             and Fuchs, Tomasz and Rhode, Wolfgang and Chirkin, Dmitry
             and Tjus, J Becker},
  journal = {Computer Physics Communications},
  volume  = {184},
  number  = {9},
  pages   = {2070--2090},
  year    = {2013},
  doi     = {10.1016/j.cpc.2013.04.001}
}
```

## Related Projects

- [PROPOSAL](https://github.com/tudo-astroparticlephysics/PROPOSAL) - The upstream C++ library
- [proposal (PyPI)](https://pypi.org/project/proposal/) - Python bindings for PROPOSAL

## Issues

For issues specific to the Julia bindings, please open an issue on this repository.
For issues with the core PROPOSAL library, please report them to the [upstream repository](https://github.com/tudo-astroparticlephysics/PROPOSAL/issues).
