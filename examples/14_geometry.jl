# Geometry Operations Example
# Demonstrates Sphere, Cylinder, Box construction and spatial queries.

using PROPOSAL

println("=== Geometry Operations ===\n")

# --- Sphere ---
println("--- Sphere ---")
origin = Cartesian3D(0, 0, 0)
sphere = Sphere(origin, 1000.0, 0.0)  # outer_radius=1000, inner_radius=0
println("  Radius: $(get_radius(sphere))")
println("  Inner radius: $(get_inner_radius(sphere))")

# Spatial queries
pos_inside = Cartesian3D(0, 0, 500)
pos_outside = Cartesian3D(0, 0, 1500)
dir = Cartesian3D(0, 0, 1)

println("  Point at z=500 inside? $(is_inside(sphere, pos_inside, dir))")
println("  Point at z=1500 inside? $(is_inside(sphere, pos_outside, dir))")
println("  Point at z=500 is_infront? $(is_infront(sphere, pos_inside, dir))")
println("  Point at z=1500 is_behind? $(is_behind(sphere, pos_outside, dir))")

# Distance to border
d1, d2 = distance_to_border(sphere, pos_inside, dir)
println("  From z=500 going +z: first border=$(round(d1, sigdigits=5)), second border=$(round(d2, sigdigits=5))")

d1, d2 = distance_to_border(sphere, Cartesian3D(0, 0, -1500), dir)
println("  From z=-1500 going +z: first border=$(round(d1, sigdigits=5)), second border=$(round(d2, sigdigits=5))")

# Hollow sphere
println("\n--- Hollow sphere ---")
hollow = Sphere(origin, 1000.0, 500.0)
println("  Outer radius: $(get_radius(hollow)), Inner radius: $(get_inner_radius(hollow))")
println("  Point at z=300 inside? $(is_inside(hollow, Cartesian3D(0, 0, 300), dir))")
println("  Point at z=700 inside? $(is_inside(hollow, Cartesian3D(0, 0, 700), dir))")

# --- Cylinder ---
println("\n--- Cylinder ---")
cyl = Cylinder(origin, 500.0, 200.0, 0.0)  # z_half=500, radius=200, inner_radius=0
println("  Radius: $(get_radius(cyl))")
println("  Half-height: $(get_z(cyl))")

pos_in = Cartesian3D(100, 0, 0)
pos_out = Cartesian3D(300, 0, 0)
dir_x = Cartesian3D(1, 0, 0)
println("  Point at x=100 inside? $(is_inside(cyl, pos_in, dir_x))")
println("  Point at x=300 inside? $(is_inside(cyl, pos_out, dir_x))")
d1, d2 = distance_to_border(cyl, pos_in, dir_x)
println("  From x=100 going +x: first=$(round(d1, sigdigits=5)), second=$(round(d2, sigdigits=5))")

# --- Box ---
println("\n--- Box ---")
box = Box(origin, 200.0, 300.0, 400.0)  # half-widths in x, y, z
println("  Half-widths: x=$(get_x(box)), y=$(get_y(box)), z=$(get_z(box))")

pos_in = Cartesian3D(50, 50, 50)
pos_out = Cartesian3D(250, 0, 0)
println("  Point at (50,50,50) inside? $(is_inside(box, pos_in, dir))")
println("  Point at (250,0,0) inside? $(is_inside(box, pos_out, dir))")

# --- Hierarchy ---
println("\n--- Geometry hierarchy ---")
set_hierarchy(sphere, 1)
println("  Sphere hierarchy: $(get_hierarchy(sphere))")
println("  Geometry name: $(get_geometry_name(sphere))")

# --- Hit detection with propagation ---
println("\n--- Hit detection with propagation ---")
pdef = MuMinusDef()
medium = create_ice()
handle = create_crosssections(pdef, medium, 500.0, 0.05, false, true)

coll = PropagationUtilityCollection()
set_displacement!(coll, make_displacement(handle, true))
set_interaction!(coll, make_interaction(handle, true))
set_time!(coll, make_time_calc(handle, pdef, true))

big_sphere = Sphere(origin, 1e20, 0.0)
dens = DensityHomogeneous(get_mass_density(medium))
prop = make_propagator_single_sector(pdef, big_sphere, dens, coll)

state = ParticleState(PARTICLE_TYPE_MUMINUS, 0, 0, 0, 0, 0, 1, 1e7)
set_random_seed(42)
sec = propagate(prop, state)

# Check if track hits a detector sphere
detector = Sphere(Cartesian3D(0, 0, 5000), 1000.0, 0.0)
does_hit = hit_geometry(sec, detector)
println("  Track hits detector sphere at z=5000, r=1000? $does_hit")
if does_hit
    if has_entry_point(sec, detector)
        ep = get_entry_point(sec, detector)
        pos = get_position(ep)
        println("  Entry point: ($(round(get_x(pos), sigdigits=5)), $(round(get_y(pos), sigdigits=5)), $(round(get_z(pos), sigdigits=5)))")
        println("  Entry energy: $(round(get_energy(ep), sigdigits=5)) MeV")
    end
    if has_exit_point(sec, detector)
        xp = get_exit_point(sec, detector)
        pos = get_position(xp)
        println("  Exit point: ($(round(get_x(pos), sigdigits=5)), $(round(get_y(pos), sigdigits=5)), $(round(get_z(pos), sigdigits=5)))")
        println("  Exit energy: $(round(get_energy(xp), sigdigits=5)) MeV")
    end
end

free_crosssections(handle)
println("\nDone.")
