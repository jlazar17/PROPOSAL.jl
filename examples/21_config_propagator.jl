# Configuration-Based Propagator Example
# Demonstrates creating propagators from JSON configuration files
# using factory functions for each particle type.

using PROPOSAL

println("=== Configuration-Based Propagator ===\n")

set_tables_path("/tmp")

# --- Write a minimal JSON config ---
config_path = joinpath(tempdir(), "proposal_config.json")

# This is a minimal configuration for a muon propagator.
# PROPOSAL reads JSON configs that specify sectors with geometry,
# medium, cuts, and cross section settings.
config_json = """
{
    "global": {
        "seed": 1234,
        "continous_loss_output": true,
        "only_loss_inside_detector": false
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
"""

open(config_path, "w") do f
    write(f, config_json)
end
println("  Wrote config to: $config_path")

# --- Create propagators from config ---
println("\n--- Factory functions for each particle type ---")

# Each particle type has a dedicated factory
println("  Available factories:")
println("    create_propagator_muminus(config_path)")
println("    create_propagator_muplus(config_path)")
println("    create_propagator_eminus(config_path)")
println("    create_propagator_eplus(config_path)")
println("    create_propagator_tauminus(config_path)")
println("    create_propagator_tauplus(config_path)")
println("    create_propagator_gamma(config_path)")
println("    create_propagator(particle_def, config_path)  # generic")

# Create a mu-minus propagator
prop = create_propagator_muminus(config_path)
println("\n  Created mu-minus propagator from config")

# --- Propagate ---
state = ParticleState(PARTICLE_TYPE_MUMINUS, 0, 0, 0, 0, 0, 1, 1e7)
set_random_seed(42)
sec = propagate(prop, state)

println("\n--- Propagation result ---")
println("  Initial energy: $(get_energy(get_initial_state(sec))) MeV")
println("  Final energy: $(round(get_energy(get_final_state(sec)), sigdigits=5)) MeV")
println("  Track length: $(round(get_track_length(sec), sigdigits=5)) cm")
println("  Stochastic losses: $(get_stochastic_losses_count(sec))")
println("  Continuous losses: $(get_continuous_losses_count(sec))")

# --- Generic factory ---
println("\n--- Generic factory with ParticleDef ---")
tau = TauMinusDef()
prop_tau = create_propagator(tau, config_path)
state_tau = ParticleState(PARTICLE_TYPE_TAUMINUS, 0, 0, 0, 0, 0, 1, 1e9)
set_random_seed(42)
sec_tau = propagate(prop_tau, state_tau)
println("  Tau propagation:")
println("  Final energy: $(round(get_energy(get_final_state(sec_tau)), sigdigits=5)) MeV")
println("  Track length: $(round(get_track_length(sec_tau), sigdigits=5)) cm")
println("  Has decay? $(has_decay(sec_tau))")

# Cleanup
rm(config_path)
println("\nDone.")
