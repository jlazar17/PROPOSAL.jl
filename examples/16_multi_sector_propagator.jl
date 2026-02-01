# Multi-Sector Propagation Example
# Demonstrates building a propagator with multiple geometry regions,
# each with its own medium/density.

using PROPOSAL

println("=== Multi-Sector Propagation ===\n")

set_tables_path("/tmp")

pdef = MuMinusDef()

# --- Setup two media ---
ice = create_ice()
rock = create_standard_rock()

# Cross sections for each medium
handle_ice = create_crosssections(pdef, ice, 500.0, 0.05, false, true)
handle_rock = create_crosssections(pdef, rock, 500.0, 0.05, false, true)

# Build utility collections for each sector
coll_ice = PropagationUtilityCollection()
set_displacement!(coll_ice, make_displacement(handle_ice, true))
set_interaction!(coll_ice, make_interaction(handle_ice, true))
set_time!(coll_ice, make_time_calc(handle_ice, pdef, true))

coll_rock = PropagationUtilityCollection()
set_displacement!(coll_rock, make_displacement(handle_rock, true))
set_interaction!(coll_rock, make_interaction(handle_rock, true))
set_time!(coll_rock, make_time_calc(handle_rock, pdef, true))

# --- Define geometry: inner sphere of ice, outer shell of rock ---
println("--- Two-sector setup: ice sphere (r<5000) + rock shell (5000<r<1e7) ---")

inner_geo = Sphere(Cartesian3D(0, 0, 0), 5000.0, 0.0)
set_hierarchy(inner_geo, 1)

outer_geo = Sphere(Cartesian3D(0, 0, 0), 1e7, 0.0)
set_hierarchy(outer_geo, 0)

dens_ice = DensityHomogeneous(get_mass_density(ice))
dens_rock = DensityHomogeneous(get_mass_density(rock))

# Build multi-sector propagator
clear_sectors()
add_sector_sphere_homogeneous(inner_geo, dens_ice, coll_ice)
add_sector_sphere_homogeneous(outer_geo, dens_rock, coll_rock)
prop = make_propagator_from_sectors(pdef)

# --- Propagate ---
println("\n--- Propagating 1e7 MeV muon from origin along +z ---")
state = ParticleState(PARTICLE_TYPE_MUMINUS, 0, 0, 0, 0, 0, 1, 1e7)

set_random_seed(42)
n_samples = 20
final_energies = Float64[]
final_distances = Float64[]
for i in 1:n_samples
    sec = propagate(prop, state, max_distance=1e8, min_energy=get_mass(pdef))
    fs = get_final_state(sec)
    push!(final_energies, get_energy(fs))
    push!(final_distances, get_propagated_distance(fs))
end

mean_d = sum(final_distances) / length(final_distances)
println("  Mean propagated distance: $(round(mean_d / 1e5, sigdigits=5)) km")
println("  Min distance: $(round(minimum(final_distances) / 1e5, sigdigits=5)) km")
println("  Max distance: $(round(maximum(final_distances) / 1e5, sigdigits=5)) km")

# --- Compare with single ice sector ---
println("\n--- Comparison: pure ice propagation ---")
big_sphere = Sphere(Cartesian3D(0, 0, 0), 1e20, 0.0)
prop_ice = make_propagator_single_sector(pdef, big_sphere, dens_ice, coll_ice)

set_random_seed(42)
distances_ice = Float64[]
for i in 1:n_samples
    sec = propagate(prop_ice, state, max_distance=1e8, min_energy=get_mass(pdef))
    push!(distances_ice, get_propagated_distance(get_final_state(sec)))
end
mean_d_ice = sum(distances_ice) / length(distances_ice)
println("  Pure ice mean distance: $(round(mean_d_ice / 1e5, sigdigits=5)) km")
println("  Two-sector mean distance: $(round(mean_d / 1e5, sigdigits=5)) km")
println("  (Rock is denser, so muons stop sooner in the two-sector case)")

free_crosssections(handle_ice)
free_crosssections(handle_rock)
println("\nDone.")
