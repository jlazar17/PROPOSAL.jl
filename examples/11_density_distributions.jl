# Density Distribution Example
# Demonstrates homogeneous, exponential, polynomial, and spline density distributions.
# This covers features from Phases 10-11 not in the Python examples.

using PROPOSAL

println("=== Density Distributions ===\n")

# --- Homogeneous density ---
println("--- Homogeneous density ---")
dens_h = DensityHomogeneous(0.917)  # ice density
pos = Cartesian3D(0, 0, 0)
dir = Cartesian3D(0, 0, 1)
println("  Evaluate at origin: $(evaluate(dens_h, pos))")
println("  Integrate 0->100cm: $(integrate(dens_h, pos, dir, 100.0))")

# --- Exponential density ---
println("\n--- Exponential density ---")
ax = CartesianAxis(Cartesian3D(0, 0, 0), Cartesian3D(0, 0, 1))
# rho(z) = massDensity * exp(z/sigma - d0)
# Evaluate returns the relative correction factor
sigma = 100.0
d0 = 0.0
dens_exp = DensityExponential(ax, sigma, d0, 1.0)
println("  At z=0:   $(evaluate(dens_exp, Cartesian3D(0, 0, 0)))")
println("  At z=50:  $(evaluate(dens_exp, Cartesian3D(0, 0, 50)))")
println("  At z=100: $(evaluate(dens_exp, Cartesian3D(0, 0, 100)))")
println("  At z=200: $(evaluate(dens_exp, Cartesian3D(0, 0, 200)))")

# --- Polynomial density ---
println("\n--- Polynomial density ---")
# Linear density profile: rho(x) = 1.0 + 0.001 * x
poly = Polynom([1.0, 0.001])
dens_poly = DensityPolynomial(ax, poly, 1.0)
println("  At z=0:    $(evaluate(dens_poly, Cartesian3D(0, 0, 0)))")
println("  At z=500:  $(evaluate(dens_poly, Cartesian3D(0, 0, 500)))")
println("  At z=1000: $(evaluate(dens_poly, Cartesian3D(0, 0, 1000)))")

# --- Spline density ---
println("\n--- Spline density ---")
# Tabulated density profile
zs = [0.0, 1000.0, 2000.0, 3000.0, 4000.0, 5000.0]
rhos = [1.0, 0.95, 0.85, 0.7, 0.5, 0.3]
spline = LinearSpline(zs, rhos)
dens_spline = DensitySplines(ax, spline, 1.0)
println("  At z=0:    $(evaluate(dens_spline, Cartesian3D(0, 0, 0)))")
println("  At z=1500: $(evaluate(dens_spline, Cartesian3D(0, 0, 1500)))")
println("  At z=3000: $(evaluate(dens_spline, Cartesian3D(0, 0, 3000)))")
println("  At z=5000: $(evaluate(dens_spline, Cartesian3D(0, 0, 5000)))")

# --- Use homogeneous density in a propagator (for comparison) ---
println("\n--- Propagation through homogeneous vs exponential density ---")
pdef = MuMinusDef()
medium = create_ice()
handle = create_crosssections(pdef, medium, 500.0, 0.05, false, true)

coll = PropagationUtilityCollection()
set_displacement!(coll, make_displacement(handle, true))
set_interaction!(coll, make_interaction(handle, true))
set_time!(coll, make_time_calc(handle, pdef, true))

geo = Sphere(Cartesian3D(0, 0, 0), 1e20, 0.0)

# Homogeneous density propagation
prop_h = make_propagator_single_sector(pdef, geo, DensityHomogeneous(get_mass_density(medium)), coll)
state = ParticleState(PARTICLE_TYPE_MUMINUS, 0, 0, 0, 0, 0, 1, 1e6)
set_random_seed(42)
sec_h = propagate(prop_h, state, max_distance=1e7, min_energy=get_mass(pdef))
println("  Homogeneous: final energy = $(round(get_energy(get_final_state(sec_h)), sigdigits=5)) MeV, dist = $(round(get_propagated_distance(get_final_state(sec_h)), sigdigits=5)) cm")

free_crosssections(handle)
println("\nDone.")
