# Advanced Propagator Example
# Julia equivalent of Python AdvancedPropagator.ipynb
#
# Build a propagator from individual components: cross sections,
# displacement, interaction, time, and geometry/density.

using PROPOSAL

println("=== Advanced Propagator ===\n")

set_tables_path("/tmp")

# --- Setup ---
pdef = MuMinusDef()
medium = create_ice()
handle = create_crosssections(pdef, medium, Inf, 0.05, false, true)

# Build the propagation utility collection
coll = PropagationUtilityCollection()
set_displacement!(coll, make_displacement(handle, true))
set_interaction!(coll, make_interaction(handle, true))
set_time!(coll, make_time_calc(handle, pdef, true))

# Define geometry and density
geo = Sphere(Cartesian3D(0, 0, 0), 1e20, 0.0)
dens = DensityHomogeneous(get_mass_density(medium))

# Create the propagator
prop = make_propagator_single_sector(pdef, geo, dens, coll)

# --- Propagate with max distance ---
println("--- Propagate 1e9 MeV muon for 1e4 cm ---")
init_state = ParticleState(PARTICLE_TYPE_MUMINUS, 0, 0, 0, 0, 0, 1, 1e9)

set_random_seed(42)
n_samples = 20
final_energies = Float64[]
for i in 1:n_samples
    track = propagate(prop, init_state, max_distance=1e4)
    push!(final_energies, get_energy(get_final_state(track)))
end

mean_E = sum(final_energies) / length(final_energies)
println("Final energies after 1e4 cm ($(n_samples) samples):")
println("  Mean: $(round(mean_E / 1e6, sigdigits=5)) TeV")
println("  Min:  $(round(minimum(final_energies) / 1e6, sigdigits=5)) TeV")
println("  Max:  $(round(maximum(final_energies) / 1e6, sigdigits=5)) TeV")

# --- Propagate until min energy ---
println("\n--- Propagate 1e6 MeV muon until E < 1e4 MeV ---")
init_state2 = ParticleState(PARTICLE_TYPE_MUMINUS, 0, 0, 0, 0, 0, 1, 1e6)

set_random_seed(42)
distances = Float64[]
for i in 1:n_samples
    track = propagate(prop, init_state2, max_distance=1e20, min_energy=1e4)
    push!(distances, get_propagated_distance(get_final_state(track)))
end

mean_d = sum(distances) / length(distances)
println("Propagated distances ($(n_samples) samples):")
println("  Mean:   $(round(mean_d / 1e5, sigdigits=5)) km")
println("  Min:    $(round(minimum(distances) / 1e5, sigdigits=5)) km")
println("  Max:    $(round(maximum(distances) / 1e5, sigdigits=5)) km")

# --- Detailed track inspection ---
println("\n--- Detailed track for a single 1e7 MeV muon ---")
init_state3 = ParticleState(PARTICLE_TYPE_MUMINUS, 0, 0, 0, 0, 0, 1, 1e7)
set_random_seed(1234)
track = propagate(prop, init_state3)

println("Track size: $(track_size(track))")
println("Initial energy: $(get_energy(get_initial_state(track))) MeV")
println("Final energy: $(get_energy(get_final_state(track))) MeV")
println("Track length: $(get_track_length(track)) cm")

# Stochastic losses
losses = get_stochastic_losses(track)
println("\nStochastic losses ($(length(losses))):")
for loss in losses
    E_loss = get_energy(loss)
    E_parent = get_parent_particle_energy(loss)
    ltype = get_type(loss)
    target = get_target_hash(loss)
    println("  type=$ltype, E_loss=$(round(E_loss, sigdigits=5)) MeV, E_parent=$(round(E_parent, sigdigits=5)) MeV, target_hash=$target")
end

# Continuous losses
cont_losses = get_continuous_losses(track)
println("\nContinuous losses ($(length(cont_losses))):")
println("  Total continuous energy loss: $(round(sum(get_energy(cl) for cl in cont_losses), sigdigits=6)) MeV")

free_crosssections(handle)
println("\nDone.")
