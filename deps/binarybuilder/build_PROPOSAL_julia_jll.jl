# Build recipe for PROPOSAL Julia wrapper (CxxWrap bindings)
# This creates PROPOSAL_julia_jll which provides the Julia bindings
#
# To test locally:
#   julia build_PROPOSAL_julia_jll.jl --deploy=local
#
# To submit to Yggdrasil, copy this to:
#   Yggdrasil/P/PROPOSAL_julia/build_tarballs.jl

using BinaryBuilder, Pkg

name = "PROPOSAL_julia"
version = v"0.1.0"

# Collection of sources
sources = [
    # The wrapper source code
    DirectorySource("./wrapper_src"),
]

# Bash recipe for building
script = raw"""
cd $WORKSPACE/srcdir

# Create build directory
mkdir -p build && cd build

# Configure with CMake
cmake ../wrapper_src \
    -DCMAKE_INSTALL_PREFIX=${prefix} \
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TARGET_TOOLCHAIN} \
    -DCMAKE_BUILD_TYPE=Release \
    -DJlCxx_DIR=${prefix}/lib/cmake/JlCxx \
    -DPROPOSAL_DIR=${prefix}/lib/cmake/PROPOSAL

# Build and install
make -j${nproc}
make install
"""

# Platforms - must match libcxxwrap_julia platforms
platforms = supported_platforms()
platforms = expand_cxxstring_abis(platforms)

# Filter to platforms supported by libcxxwrap_julia
filter!(p -> libc(p) != "musl", platforms)  # CxxWrap doesn't support musl well

# Products
products = [
    LibraryProduct("libPROPOSAL_jl", :libPROPOSAL_jl),
]

# Dependencies
dependencies = [
    Dependency("PROPOSAL_jll"; compat="7.6"),
    Dependency("libcxxwrap_julia_jll"; compat="0.12"),
    BuildDependency("libjulia_jll"),
]

# Build the tarballs
build_tarballs(ARGS, name, version, sources, script, platforms, products, dependencies;
               julia_compat="1.6",
               preferred_gcc_version=v"9")
