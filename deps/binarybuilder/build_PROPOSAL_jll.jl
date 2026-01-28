# Build recipe for PROPOSAL C++ library
# This creates PROPOSAL_jll which provides the core PROPOSAL library
#
# To test locally:
#   julia build_PROPOSAL_jll.jl --deploy=local
#
# To submit to Yggdrasil, copy this to:
#   Yggdrasil/P/PROPOSAL/build_tarballs.jl

using BinaryBuilder, Pkg

name = "PROPOSAL"
version = v"7.6.2"

# Collection of sources
sources = [
    GitSource("https://github.com/tudo-astroparticlephysics/PROPOSAL.git",
              "8fa1a04e87b54f64ef60b77efde2d84f81d7ace4"),  # v7.6.2 tag
]

# Bash recipe for building
script = raw"""
cd $WORKSPACE/srcdir/PROPOSAL

# Create build directory
mkdir -p build && cd build

# Configure with CMake
cmake .. \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Release \
    -DBUILD_PYTHON=OFF \
    -DBUILD_TESTING=OFF \
    -DBUILD_DOCUMENTATION=OFF \
    -DBUILD_SHARED_LIBS=ON

# Build and install
make -j${nproc}
make install
"""

# Platforms to build for
platforms = supported_platforms()

# Filter out platforms that don't work well with C++14
platforms = expand_cxxstring_abis(platforms)

# Products that this package provides
products = [
    LibraryProduct("libPROPOSAL", :libPROPOSAL),
]

# Dependencies
dependencies = [
    # CubicInterpolation - may need to create this JLL first if it doesn't exist
    # For now, we'll build it as part of PROPOSAL or use header-only mode
    Dependency("boost_jll"; compat="1.79.0"),
    Dependency("spdlog_jll"),
    Dependency("nlohmann_json_jll"; compat="3.11"),
    Dependency("Eigen_jll"),
    BuildDependency("CubicInterpolation_jll"),  # Header-only at link time
]

# Build the tarballs
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
               julia_compat="1.6",
               preferred_gcc_version=v"9")
