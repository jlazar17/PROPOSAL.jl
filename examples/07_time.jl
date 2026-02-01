# Time Calculation Example
# Julia equivalent of Python Time.ipynb
#
# Computes the time elapsed during propagation, comparing exact and
# approximate time calculations.

using PROPOSAL

println("=== Time Calculation ===\n")

set_tables_path("/tmp")

pdef = MuMinusDef()
medium = create_standard_rock()
handle = create_crosssections(pdef, medium, 500.0, 1.0, true, true)

# Exact time calculator
time_calc = make_time_calc(handle, pdef, true)

# Approximate time calculator (assumes v ≈ c)
time_approx = make_time_approximate()

# Displacement calculator (needed to get grammage)
disp = make_displacement(handle, true)

density = get_mass_density(medium)
E_f = get_mass(pdef)  # muon mass

println("--- Time to propagate from E_i to rest mass ---")
println("Medium: StandardRock, density=$(density) g/cm^3")
println("E_f = $(E_f) MeV (muon mass)\n")

println("  E_i (MeV)     | Exact time (s)    | Approx time (s)   | Ratio")
println("  " * "-"^70)

for E_i in [1e3, 3e3, 1e4, 3e4, 1e5, 3e5, 1e6]
    grammage = solve_track_integral(disp, E_i, E_f)
    t_exact = time_elapsed(time_calc, E_i, E_f, grammage, density)
    t_approx = time_elapsed(time_approx, E_i, E_f, grammage, density)
    ratio = t_exact / t_approx
    println("  $(lpad(E_i, 12))  | $(lpad(round(t_exact, sigdigits=6), 17)) | $(lpad(round(t_approx, sigdigits=6), 17)) | $(round(ratio, sigdigits=6))")
end

println("\nNote: The ratio deviates from 1 at low energies where v < c.")

free_crosssections(handle)
println("\nDone.")
