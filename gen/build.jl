# Build script for generating PROPOSAL Julia wrappers
# This script:
# 1. Runs WrapIt to generate CxxWrap wrapper code
# 2. Compiles the wrapper library using CMake

using Pkg

# Ensure CxxWrap is available
if !haskey(Pkg.project().dependencies, "CxxWrap")
    Pkg.add("CxxWrap")
end

using CxxWrap

const SCRIPT_DIR = @__DIR__
const PROJECT_DIR = dirname(SCRIPT_DIR)
const BUILD_DIR = joinpath(PROJECT_DIR, "build")
const PROPOSAL_DIR = get(ENV, "PROPOSAL_DIR", joinpath(dirname(PROJECT_DIR), "PROPOSAL_upstream"))

# Get CxxWrap paths
const JLCXX_DIR = joinpath(dirname(dirname(CxxWrap.jlcxx_path)), "lib", "cmake", "JlCxx")

println("Configuration:")
println("  Script dir: $SCRIPT_DIR")
println("  Project dir: $PROJECT_DIR")
println("  Build dir: $BUILD_DIR")
println("  PROPOSAL dir: $PROPOSAL_DIR")
println("  JlCxx dir: $JLCXX_DIR")

# Step 1: Generate wrapper code using WrapIt (if available)
function run_wrapit()
    wit_config = joinpath(SCRIPT_DIR, "PROPOSAL.wit")

    # Check if wrapit is available
    wrapit_path = Sys.which("wrapit")
    if wrapit_path === nothing
        @warn "WrapIt not found in PATH. Please install WrapIt or generate wrappers manually."
        @info "Install WrapIt from: https://github.com/grasph/wrapit"
        return false
    end

    proposal_include = joinpath(PROPOSAL_DIR, "src", "PROPOSAL")

    println("\nRunning WrapIt...")
    cmd = `$wrapit_path $wit_config -I$proposal_include`
    try
        run(cmd)
        return true
    catch e
        @error "WrapIt failed" exception=e
        return false
    end
end

# Step 2: Build the wrapper library using CMake
function build_wrapper()
    mkpath(BUILD_DIR)

    cmake_args = [
        "-DCMAKE_BUILD_TYPE=Release",
        "-DJlCxx_DIR=$JLCXX_DIR",
        "-DPROPOSAL_SOURCE_DIR=$PROPOSAL_DIR",
    ]

    println("\nConfiguring CMake...")
    cd(BUILD_DIR) do
        run(`cmake $(cmake_args) $SCRIPT_DIR`)

        println("\nBuilding...")
        run(`cmake --build . --parallel`)
    end

    println("\nBuild complete!")
    println("Library location: $(joinpath(BUILD_DIR, "lib"))")
end

# Main execution
if abspath(PROGRAM_FILE) == @__FILE__
    if "--generate" in ARGS || "-g" in ARGS
        run_wrapit()
    elseif "--build" in ARGS || "-b" in ARGS
        build_wrapper()
    elseif "--all" in ARGS || "-a" in ARGS
        run_wrapit() && build_wrapper()
    else
        println("""
        PROPOSAL.jl Wrapper Build Script

        Usage: julia build.jl [options]

        Options:
          -g, --generate    Generate wrapper code using WrapIt
          -b, --build       Build the wrapper library using CMake
          -a, --all         Generate and build

        Environment variables:
          PROPOSAL_DIR      Path to PROPOSAL source directory
                           (default: ../PROPOSAL_upstream)
        """)
    end
end
