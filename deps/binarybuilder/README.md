# BinaryBuilder Recipes for PROPOSAL.jl

This directory contains BinaryBuilder recipes for building PROPOSAL and its Julia bindings.

## Overview

There are two JLL packages needed:

1. **PROPOSAL_jll** - The core PROPOSAL C++ library
2. **PROPOSAL_julia_jll** - The CxxWrap Julia bindings

## Building Locally

### Prerequisites

```bash
# Install BinaryBuilder
julia -e 'using Pkg; Pkg.add("BinaryBuilder")'
```

### Build PROPOSAL_jll

```bash
cd deps/binarybuilder
julia build_PROPOSAL_jll.jl --deploy=local
```

### Build PROPOSAL_julia_jll

After PROPOSAL_jll is built:

```bash
julia build_PROPOSAL_julia_jll.jl --deploy=local
```

## Submitting to Yggdrasil

To make these packages available in the Julia General registry:

1. Fork [Yggdrasil](https://github.com/JuliaPackaging/Yggdrasil)

2. Create the PROPOSAL recipe:
   ```
   Yggdrasil/P/PROPOSAL/build_tarballs.jl
   ```

3. Create the PROPOSAL_julia recipe:
   ```
   Yggdrasil/P/PROPOSAL_julia/build_tarballs.jl
   ```

4. Submit a PR to Yggdrasil

## Dependencies

### PROPOSAL_jll Dependencies

These JLL packages are required (most already exist in Yggdrasil):

- `boost_jll` ✓ (exists)
- `spdlog_jll` ✓ (exists)
- `nlohmann_json_jll` ✓ (exists)
- `Eigen_jll` ✓ (exists)
- `CubicInterpolation_jll` ⚠️ (may need to be created)

### PROPOSAL_julia_jll Dependencies

- `PROPOSAL_jll` (created above)
- `libcxxwrap_julia_jll` ✓ (exists)
- `libjulia_jll` ✓ (exists, build dependency)

## Creating CubicInterpolation_jll

If CubicInterpolation_jll doesn't exist, you'll need to create it first:

```julia
# build_CubicInterpolation_jll.jl
using BinaryBuilder

name = "CubicInterpolation"
version = v"0.1.5"

sources = [
    GitSource("https://github.com/MaxSac/cubic_interpolation.git", "...")
]

# ... (similar pattern to PROPOSAL_jll)
```

## Testing

After building locally with `--deploy=local`, test in Julia:

```julia
using PROPOSAL_jll
using PROPOSAL_julia_jll
using PROPOSAL

# Should now work!
is_library_available()  # true
```

## Troubleshooting

### Architecture Issues

If you see architecture mismatch errors, ensure all dependencies are built for the same architecture. BinaryBuilder handles this automatically for cross-compilation.

### Missing Symbols

If you get undefined symbol errors, check that:
1. PROPOSAL_jll was built with `BUILD_SHARED_LIBS=ON`
2. All dependencies are properly linked
3. The wrapper is linking against the correct PROPOSAL version

## File Structure

```
deps/binarybuilder/
├── README.md                      # This file
├── build_PROPOSAL_jll.jl          # Recipe for C++ library
├── build_PROPOSAL_julia_jll.jl    # Recipe for Julia wrapper
└── wrapper_src/
    ├── CMakeLists.txt             # CMake for wrapper
    └── PROPOSAL_jl.cxx            # CxxWrap bindings
```
