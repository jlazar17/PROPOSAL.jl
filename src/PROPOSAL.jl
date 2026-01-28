module PROPOSAL

using CxxWrap

# Try to load the JLL package if available, otherwise use local development path
const _using_jll = try
    @eval using PROPOSAL_cxxwrap_jll
    true
catch
    false
end

function get_lib_path()
    if _using_jll
        return PROPOSAL_cxxwrap_jll.libPROPOSAL_cxxwrap_path
    else
        # Development fallback - check local build
        if Sys.isapple()
            ext = "dylib"
        elseif Sys.iswindows()
            ext = "dll"
        else
            ext = "so"
        end

        dev_libpath = joinpath(@__DIR__, "..", "build", "lib")
        dev_lib = joinpath(dev_libpath, "libPROPOSAL_cxxwrap.$ext")

        if isfile(dev_lib)
            return dev_lib
        end

        # Check environment variable
        env_path = get(ENV, "PROPOSAL_JL_LIB_PATH", "")
        if !isempty(env_path)
            env_lib = joinpath(env_path, "libPROPOSAL_cxxwrap.$ext")
            if isfile(env_lib)
                return env_lib
            end
        end

        return nothing
    end
end

const _lib_path = get_lib_path()
const _library_available = _lib_path !== nothing

if _library_available
    @wrapmodule(() -> _lib_path)

    function __init__()
        @initcxx
    end
else
    function __init__()
        @warn """
        PROPOSAL wrapper library not found.

        To use PROPOSAL.jl, either:
        1. Install PROPOSAL_cxxwrap_jll (when available in registry)
        2. Build locally: see deps/binarybuilder/README.md
        3. Set PROPOSAL_JL_LIB_PATH environment variable

        The package structure is ready but bindings are not loaded.
        """
    end
end

# Check if library is available
is_library_available() = _library_available
export is_library_available

# Exports - these will be available when the library is loaded
if _library_available
    # Math types
    export Cartesian3D, Spherical3D
    export get_x, get_y, get_z, set_x, set_y, set_z, magnitude
    export get_radius, get_azimuth, get_zenith, set_radius, set_azimuth, set_zenith

    # ParticleState accessors
    export get_type, get_position, get_direction, get_energy, get_time, get_propagated_distance

    # Particle definitions
    export ParticleDef, ParticleState
    export MuMinusDef, MuPlusDef, EMinusDef, EPlusDef
    export TauMinusDef, TauPlusDef, GammaDef

    # Propagation
    export EnergyCutSettings, Medium, Secondaries, Propagator
    export track_size, get_track_state, get_track_length
    export get_initial_state, get_final_state
    export get_track_energies_array, get_track_times_array, get_track_propagated_distances_array
    export has_decay, get_decay_products_to_array, get_total_continuous_energy_loss

    # Functions
    export make_particle_state, create_medium
    export set_random_seed, random_double, get_proposal_version
    export create_propagator_muminus, create_propagator_muplus
    export create_propagator_eminus, create_propagator_eplus
    export create_propagator_tauminus, create_propagator_tauplus
    export create_propagator_gamma

    # Constants
    export SPEED_OF_LIGHT, ELECTRON_MASS, MUON_MASS, TAU_MASS, PROTON_MASS
    export PARTICLE_TYPE_NONE, PARTICLE_TYPE_EMINUS, PARTICLE_TYPE_EPLUS
    export PARTICLE_TYPE_MUMINUS, PARTICLE_TYPE_MUPLUS
    export PARTICLE_TYPE_TAUMINUS, PARTICLE_TYPE_TAUPLUS, PARTICLE_TYPE_GAMMA
end

# High-level Julia API (available when library is loaded)
if _library_available
    """
        propagate(propagator, state; max_distance=1e20, min_energy=0.0)

    Propagate a particle through the configured medium.

    # Arguments
    - `propagator`: A `Propagator` instance
    - `state`: Initial `ParticleState`
    - `max_distance`: Maximum propagation distance in cm (default: 1e20)
    - `min_energy`: Minimum energy in MeV to stop propagation (default: 0.0)

    # Returns
    - `Secondaries` object containing the track and secondary particles
    """
    function propagate(prop::Propagator, state::ParticleState;
                       max_distance::Real=1e20, min_energy::Real=0.0)
        return PROPOSAL.propagate(prop, state, Float64(max_distance), Float64(min_energy))
    end
    export propagate

    """
        ParticleState(particle_type, x, y, z, dx, dy, dz, energy; time=0.0, propagated_distance=0.0)

    Create a new particle state.

    # Arguments
    - `particle_type`: Integer particle type (use PARTICLE_TYPE_* constants)
    - `x, y, z`: Initial position in cm
    - `dx, dy, dz`: Initial direction (will be used as-is, normalize if needed)
    - `energy`: Initial energy in MeV
    - `time`: Initial time in seconds (default: 0.0)
    - `propagated_distance`: Initial propagated distance in cm (default: 0.0)

    # Example
    ```julia
    state = ParticleState(PARTICLE_TYPE_MUMINUS, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1e8)
    ```
    """
    function ParticleState(particle_type::Integer,
                          x::Real, y::Real, z::Real,
                          dx::Real, dy::Real, dz::Real,
                          energy::Real;
                          time::Real=0.0,
                          propagated_distance::Real=0.0)
        return make_particle_state(Int(particle_type),
                                   Float64(x), Float64(y), Float64(z),
                                   Float64(dx), Float64(dy), Float64(dz),
                                   Float64(energy), Float64(time),
                                   Float64(propagated_distance))
    end
end

end # module
