# Secondaries Calculator Example
# Julia equivalent of Python Secondaries.ipynb
#
# The SecondariesCalculator computes the secondary particles produced
# in a stochastic interaction (e.g., the electron-positron pair from
# pair production, the photon from bremsstrahlung).

using PROPOSAL

println("=== Secondaries Calculator ===\n")

pdef = MuMinusDef()
medium = create_air()

# Create cross sections to initialize the secondaries calculator
handle = create_crosssections(pdef, medium, 500.0, 0.05, false, false)
sec_handle = create_secondaries_calculator(handle, pdef, medium)

# --- Check required random numbers for each interaction type ---
println("--- Required random numbers per interaction type ---")
for (name, itype) in [
    ("Brems", INTERACTION_TYPE_BREMS),
    ("Epair", INTERACTION_TYPE_EPAIR),
    ("Ioniz", INTERACTION_TYPE_IONIZ),
    ("Photonuclear", INTERACTION_TYPE_PHOTONUCLEAR),
]
    n_rnd = sec_calc_required_random_numbers(sec_handle, itype)
    println("  $name: $n_rnd random numbers")
end

println("\nNote: The SecondariesCalculator is used internally by the Propagator")
println("to compute daughter particles. StochasticLoss objects come from propagation.")

# --- Demonstrate with a full propagation ---
println("\n--- Secondaries from a full propagation ---")
medium_ice = create_ice()
handle2 = create_crosssections(pdef, medium_ice, 500.0, 0.05, false, true)

# Build propagator
coll = PropagationUtilityCollection()
set_displacement!(coll, make_displacement(handle2, true))
set_interaction!(coll, make_interaction(handle2, true))
set_time!(coll, make_time_calc(handle2, pdef, true))

geo = Sphere(Cartesian3D(0, 0, 0), 1e20, 0.0)
dens = DensityHomogeneous(get_mass_density(medium_ice))
prop = make_propagator_single_sector(pdef, geo, dens, coll)

state = ParticleState(PARTICLE_TYPE_MUMINUS, 0, 0, 0, 0, 0, 1, 1e7)
set_random_seed(1234)
sec = propagate(prop, state)

# Print stochastic losses
losses = get_stochastic_losses(sec)
println("Number of stochastic losses: $(length(losses))")
for (i, loss) in enumerate(losses)
    E_loss = get_energy(loss)
    E_parent = get_parent_particle_energy(loss)
    ltype = get_type(loss)
    println("  Loss $i: type=$ltype, E_loss=$(round(E_loss, sigdigits=5)) MeV, E_parent=$(round(E_parent, sigdigits=5)) MeV")
end

free_secondaries_calculator(sec_handle)
free_crosssections(handle)
free_crosssections(handle2)
println("\nDone.")
