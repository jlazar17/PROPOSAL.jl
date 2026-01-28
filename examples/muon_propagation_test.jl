#!/usr/bin/env julia
"""
PROPOSAL.jl Muon Propagation Test

Propagate 1000 muons at logarithmically spaced energies through ice
and report mean and median propagation distances.
"""

using PROPOSAL
using Statistics
using Printf

# Check if library is available
if !is_library_available()
    error("PROPOSAL library not available. Please build or install PROPOSAL_cxxwrap_jll.")
end

# Configuration for ice medium
const ICE_CONFIG = """
{
    "global": {
        "cuts": {
            "e_cut": 500,
            "v_cut": 0.05,
            "cont_rand": true
        }
    },
    "sectors": [
        {
            "medium": "ice",
            "geometries": [
                {
                    "hierarchy": 0,
                    "shape": "sphere",
                    "origin": [0, 0, 0],
                    "outer_radius": 1e20
                }
            ]
        }
    ]
}
"""

function setup_propagator()
    # Create temporary config file
    config_path = joinpath(tempdir(), "proposal_ice_config.json")
    write(config_path, ICE_CONFIG)

    # Create propagator for muon minus
    prop = create_propagator_muminus(config_path)
    return prop
end

function propagate_muon(prop, energy_mev)
    """
    Propagate a single muon with given energy.
    Returns the total propagation distance in cm.
    """
    # Create initial state: muon at origin, pointing in +z direction
    initial_state = ParticleState(
        PARTICLE_TYPE_MUMINUS,
        0.0, 0.0, 0.0,      # position (cm)
        0.0, 0.0, 1.0,      # direction
        energy_mev;         # energy (MeV)
        time=0.0,
        propagated_distance=0.0
    )

    # Propagate until muon stops or reaches max distance (1000 km)
    max_distance = 1e8  # 1000 km in cm
    min_energy = MUON_MASS  # Stop at rest mass

    secondaries = propagate(prop, initial_state, max_distance=max_distance, min_energy=min_energy)

    # Get the final propagation distance
    final_state = get_final_state(secondaries)
    return get_propagated_distance(final_state)
end

function run_propagation_study(; n_muons_total=1000, e_min=1e3, e_max=1e9, n_energies=20, seed=42)
    """
    Run muon propagation study.

    Parameters:
    - n_muons_total: total number of muons to propagate (distributed across energies)
    - e_min: minimum energy in MeV (default: 1 GeV)
    - e_max: maximum energy in MeV (default: 1 PeV)
    - n_energies: number of energy points
    - seed: random seed for reproducibility
    """

    n_muons_per_energy = max(1, n_muons_total ÷ n_energies)

    println("=" ^ 70)
    println("PROPOSAL.jl Muon Propagation Test")
    println("=" ^ 70)
    println()
    println("PROPOSAL version: $(get_proposal_version())")
    println("Medium: Ice")
    println("Particle: Muon (mu-)")
    println()
    println("Parameters:")
    println("  Total muons: $(n_muons_per_energy * n_energies)")
    println("  Muons per energy: $n_muons_per_energy")
    println("  Energy range: $(e_min/1e3) GeV - $(e_max/1e6) TeV")
    println("  Number of energy points: $n_energies")
    println("  Random seed: $seed")
    println()

    # Set random seed
    set_random_seed(seed)

    # Setup propagator
    println("Setting up propagator...")
    prop = setup_propagator()
    println("Propagator ready.")
    println()

    # Generate logarithmically spaced energies
    energies = 10 .^ range(log10(e_min), log10(e_max), length=n_energies)

    # Store results
    results = Dict{Float64, Vector{Float64}}()

    println("-" ^ 70)
    @printf("%-15s %15s %15s %15s\n", "Energy", "Mean Distance", "Median Distance", "Std Dev")
    @printf("%-15s %15s %15s %15s\n", "(GeV)", "(km)", "(km)", "(km)")
    println("-" ^ 70)

    for energy in energies
        distances = Float64[]

        for i in 1:n_muons_per_energy
            dist = propagate_muon(prop, energy)
            push!(distances, dist)
        end

        results[energy] = distances

        # Convert cm to km for display
        distances_km = distances ./ 1e5

        mean_dist = mean(distances_km)
        median_dist = median(distances_km)
        std_dist = std(distances_km)

        @printf("%-15.2e %15.3f %15.3f %15.3f\n",
                energy/1e3,  # Convert MeV to GeV
                mean_dist, median_dist, std_dist)
    end

    println("-" ^ 70)
    println()

    # Summary statistics across all energies
    println("=" ^ 70)
    println("Summary")
    println("=" ^ 70)

    all_distances = vcat(values(results)...)
    all_distances_km = all_distances ./ 1e5

    println()
    @printf("Total muons propagated: %d\n", length(all_distances))
    @printf("Overall mean distance:  %.3f km\n", mean(all_distances_km))
    @printf("Overall median distance: %.3f km\n", median(all_distances_km))
    @printf("Min distance: %.3f km\n", minimum(all_distances_km))
    @printf("Max distance: %.3f km\n", maximum(all_distances_km))
    println()

    return results
end

# Run the study if this file is executed directly
if abspath(PROGRAM_FILE) == @__FILE__
    # Parse command line arguments
    n_muons_total = 1000
    if length(ARGS) >= 1
        n_muons_total = parse(Int, ARGS[1])
    end

    results = run_propagation_study(n_muons_total=n_muons_total)
end
