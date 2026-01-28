# Basic PROPOSAL.jl Example
# This example demonstrates how to propagate a muon through rock

using PROPOSAL

# Set random seed for reproducibility
set_random_seed(42)

# Create a propagator for muon minus using a configuration file
# Note: You need a valid PROPOSAL configuration JSON file
config_path = joinpath(@__DIR__, "config_minimal.json")

if !isfile(config_path)
    @warn "Configuration file not found. Creating a minimal example config."

    # Create a minimal configuration
    config = """
    {
        "global": {
            "medium": "StandardRock",
            "cuts": {
                "e_cut": 500,
                "v_cut": 0.05,
                "cont_rand": false
            },
            "do_interpolation": true,
            "do_exact_time": false
        },
        "sectors": [
            {
                "geometry": {
                    "shape": "sphere",
                    "origin": [0, 0, 0],
                    "outer_radius": 1e20,
                    "inner_radius": 0
                }
            }
        ]
    }
    """
    write(config_path, config)
end

# Create particle definition for muon minus
particle_def = MuMinusDef()

# Create propagator
prop = Propagator(particle_def, config_path)

# Create initial particle state
# Muon at origin, pointing in +z direction, with 100 TeV energy
initial_position = Cartesian3D(0.0, 0.0, 0.0)
initial_direction = Cartesian3D(0.0, 0.0, 1.0)
initial_energy = 1e8  # 100 TeV in MeV

initial_state = ParticleState(
    ParticleType_MuMinus,
    initial_position,
    initial_direction,
    initial_energy,
    0.0,  # time
    0.0   # propagated distance
)

# Propagate the particle
println("Propagating muon with E = $(initial_energy/1e6) TeV...")
secondaries = propagate(prop, initial_state, max_distance=1e7)  # 100 km

# Get the track (list of particle states along the trajectory)
track = get_track(secondaries)

println("\nTrack has $(length(track)) points")
println("Entry energy: $(get_entry_energy(secondaries)/1e6) TeV")
println("Exit energy: $(get_exit_energy(secondaries)/1e6) TeV")

# Print some track information
if length(track) > 0
    entry = get_entry_point(secondaries)
    exit_pt = get_exit_point(secondaries)

    println("\nEntry point: ($(get_x(entry)), $(get_y(entry)), $(get_z(entry))) cm")
    println("Exit point: ($(get_x(exit_pt)), $(get_y(exit_pt)), $(get_z(exit_pt))) cm")
end

# Get secondary particles produced
secs = get_secondaries(secondaries)
println("\nNumber of secondaries: $(length(secs))")

println("\nPROPOSAL version: $(get_version())")
