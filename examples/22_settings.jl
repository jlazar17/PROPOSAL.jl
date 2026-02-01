# InterpolationSettings and PropagationSettings Example
# Demonstrates tuning interpolation nodes, table paths, energy limits,
# and max propagation steps.

using PROPOSAL

println("=== Interpolation & Propagation Settings ===\n")

# =========================================================================
# InterpolationSettings
# =========================================================================
# PROPOSAL uses interpolation tables to speed up cross section evaluation.
# These settings control the accuracy and storage of those tables.

println("--- InterpolationSettings (defaults) ---")
println("  Tables path: \"$(get_tables_path())\"")
println("  Upper energy limit: $(get_upper_energy_lim()) MeV")
println("  Nodes dEdx:     $(get_nodes_dedx())")
println("  Nodes dE2dx:    $(get_nodes_de2dx())")
println("  Nodes dNdx (E): $(get_nodes_dndx_e())")
println("  Nodes dNdx (v): $(get_nodes_dndx_v())")
println("  Nodes utility:  $(get_nodes_utility())")
println("  Nodes rate:     $(get_nodes_rate_interpolant())")

# --- Modify settings ---
println("\n--- Modifying settings ---")

# Tables path: where interpolation tables are cached
original_path = get_tables_path()
set_tables_path("/tmp/proposal_tables")
println("  Set tables path to: $(get_tables_path())")

# Upper energy limit for interpolation tables
original_upper = get_upper_energy_lim()
set_upper_energy_lim(1e14)
println("  Set upper energy limit to: $(get_upper_energy_lim()) MeV")

# Number of interpolation nodes (higher = more accurate but slower to build)
original_dedx = get_nodes_dedx()
set_nodes_dedx(200)
println("  Set dEdx nodes to: $(get_nodes_dedx())")

original_dndx_e = get_nodes_dndx_e()
set_nodes_dndx_e(150)
println("  Set dNdx energy nodes to: $(get_nodes_dndx_e())")

original_dndx_v = get_nodes_dndx_v()
set_nodes_dndx_v(80)
println("  Set dNdx v nodes to: $(get_nodes_dndx_v())")

# Restore defaults
set_tables_path(original_path)
set_upper_energy_lim(original_upper)
set_nodes_dedx(original_dedx)
set_nodes_dndx_e(original_dndx_e)
set_nodes_dndx_v(original_dndx_v)
println("\n  Restored all interpolation settings to defaults.")

# =========================================================================
# PropagationSettings
# =========================================================================
println("\n--- PropagationSettings ---")
println("  Max steps (default): $(get_max_steps())")

# The max_steps setting limits the number of propagation steps
# to prevent infinite loops. Increase for very long tracks.
original_steps = get_max_steps()
set_max_steps(100000)
println("  Set max steps to: $(get_max_steps())")

# Restore
set_max_steps(original_steps)
println("  Restored to: $(get_max_steps())")

# =========================================================================
# Practical example: effect of interpolation on performance
# =========================================================================
println("\n--- Interpolated vs non-interpolated cross sections ---")
set_tables_path("/tmp")

pdef = MuMinusDef()
medium = create_ice()

# Non-interpolated (exact, slower)
handle_exact = create_crosssections(pdef, medium, 500.0, 0.05, false, false)
# Interpolated (fast, approximate)
handle_interp = create_crosssections(pdef, medium, 500.0, 0.05, false, true)

E = 1e6
println("  dEdx at E=1e6 MeV:")
for i in 0:(crosssection_count(handle_exact) - 1)
    xs_e = crosssection_at(handle_exact, i)
    xs_i = crosssection_at(handle_interp, i)
    dedx_e = calculate_dEdx(xs_e, E)
    dedx_i = calculate_dEdx(xs_i, E)
    itype = get_interaction_type(xs_e)
    if dedx_e != 0.0
        rel_diff = abs(dedx_i - dedx_e) / abs(dedx_e) * 100
        println("  type=$itype: exact=$(round(dedx_e, sigdigits=6)), interp=$(round(dedx_i, sigdigits=6)), diff=$(round(rel_diff, sigdigits=3))%")
    end
end

free_crosssections(handle_exact)
free_crosssections(handle_interp)
println("\nDone.")
