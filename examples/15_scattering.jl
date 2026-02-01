# Scattering Example
# Demonstrates multiple scattering parametrizations, stochastic deflection,
# and ScatteringMultiplier for scaling deflection angles.

using PROPOSAL

println("=== Scattering ===\n")

set_tables_path("/tmp")

pdef = MuMinusDef()
medium = create_ice()

# --- Create scattering objects ---
println("--- Scattering types ---")

# MS-only scattering (Highland)
scat_highland = create_scattering_ms_only("highland", pdef, medium)
println("  Created Highland MS-only scattering")
println("  MS random numbers needed: $(scattering_ms_random_numbers(scat_highland))")

# MS-only Moliere
scat_moliere = create_scattering_ms_only("moliere", pdef, medium)
println("  Created Moliere MS-only scattering")

# MS + stochastic deflection for bremsstrahlung
scat_with_sd = create_scattering_with_sd("highland", pdef, medium, "bremsginneken")
println("  Created Highland + BremsGinneken stochastic deflection")

# Scattering with multiple stochastic deflection types
sd_types = Int32[Int32(INTERACTION_TYPE_BREMS), Int32(INTERACTION_TYPE_EPAIR)]
scat_multi = create_scattering_by_types("highland", pdef, medium, sd_types)
println("  Created Highland + Brems+Epair stochastic deflection")

# --- ScatteringMultiplier ---
println("\n--- ScatteringMultiplier ---")
sd_types = Int32[Int32(INTERACTION_TYPE_BREMS), Int32(INTERACTION_TYPE_EPAIR)]
sd_mults = Float64[2.0, 1.5]
scat_scaled = create_scattering_multiplier("highland", pdef, medium, 0.5, sd_types, sd_mults)
println("  Created ScatteringMultiplier: MS×0.5, Brems SD×2.0, Epair SD×1.5")

# --- Use scattering in a propagator ---
println("\n--- Propagation with scattering ---")
handle = create_crosssections(pdef, medium, 500.0, 0.05, false, true)

coll = PropagationUtilityCollection()
set_displacement!(coll, make_displacement(handle, true))
set_interaction!(coll, make_interaction(handle, true))
set_time!(coll, make_time_calc(handle, pdef, true))
set_scattering!(coll, scat_highland)

geo = Sphere(Cartesian3D(0, 0, 0), 1e20, 0.0)
dens = DensityHomogeneous(get_mass_density(medium))
prop = make_propagator_single_sector(pdef, geo, dens, coll)

# Propagate and check angular deflection
state = ParticleState(PARTICLE_TYPE_MUMINUS, 0, 0, 0, 0, 0, 1, 1e6)
set_random_seed(42)

n_samples = 50
final_dirs = []
for i in 1:n_samples
    sec = propagate(prop, state, max_distance=1e5, min_energy=get_mass(pdef))
    fs = get_final_state(sec)
    push!(final_dirs, get_direction(fs))
end

# Compute angular spread
angles = [acos(clamp(get_z(d), -1.0, 1.0)) for d in final_dirs]
mean_angle = sum(angles) / length(angles)
println("  Mean deflection angle after 1e5 cm: $(round(rad2deg(mean_angle), sigdigits=4))°")
println("  Max deflection: $(round(rad2deg(maximum(angles)), sigdigits=4))°")
println("  Min deflection: $(round(rad2deg(minimum(angles)), sigdigits=4))°")

# --- Compare Highland vs Moliere ---
println("\n--- Highland vs Moliere comparison ---")

coll2 = PropagationUtilityCollection()
set_displacement!(coll2, make_displacement(handle, true))
set_interaction!(coll2, make_interaction(handle, true))
set_time!(coll2, make_time_calc(handle, pdef, true))
set_scattering!(coll2, scat_moliere)

prop2 = make_propagator_single_sector(pdef, geo, dens, coll2)

set_random_seed(42)
angles_moliere = Float64[]
for i in 1:n_samples
    sec = propagate(prop2, state, max_distance=1e5, min_energy=get_mass(pdef))
    fs = get_final_state(sec)
    d = get_direction(fs)
    push!(angles_moliere, acos(clamp(get_z(d), -1.0, 1.0)))
end
mean_moliere = sum(angles_moliere) / length(angles_moliere)
println("  Highland mean angle: $(round(rad2deg(mean_angle), sigdigits=4))°")
println("  Moliere mean angle: $(round(rad2deg(mean_moliere), sigdigits=4))°")

free_crosssections(handle)
println("\nDone.")
