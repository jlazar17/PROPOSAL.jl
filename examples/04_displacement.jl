# Displacement Example
# Julia equivalent of Python Displacement.ipynb
#
# The displacement utility integrates the continuous energy loss to compute
# how far a particle travels to lose a given amount of energy, and vice versa.

using PROPOSAL

println("=== Displacement ===\n")

set_tables_path("/tmp")

pdef = MuMinusDef()
medium = create_ice()
handle = create_crosssections(pdef, medium, Inf, 1.0, false, true)
disp = make_displacement(handle, false)

# --- Solve track integral: energy_initial -> energy_final gives grammage ---
println("--- Grammage to lose half the energy ---")
for E_i in [1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9]
    E_f = 0.5 * E_i
    grammage = solve_track_integral(disp, E_i, E_f)  # g/cm^2
    println("  E_i=$(E_i) MeV -> grammage = $(round(grammage, sigdigits=5)) g/cm^2 = $(round(grammage/100, sigdigits=5)) m.w.e.")
end

# --- Upper limit track integral: energy + distance -> final energy ---
println("\n--- Energy after 1000 g/cm^2 ---")
distance = 1e3  # g/cm^2
for E_i in [1e4, 1e5, 1e6, 1e7, 1e8, 1e9]
    E_f = upper_limit_track_integral(disp, E_i, distance)
    loss = E_i - E_f
    println("  E_i=$(E_i) MeV -> E_f=$(round(E_f, sigdigits=5)) MeV, loss=$(round(loss, sigdigits=5)) MeV")
end

# --- Specific example from Python notebook ---
energy = 1e6
distance_half = solve_track_integral(disp, energy, 0.5 * energy)
println("\nA $(energy/1e3) GeV muon propagates $(round(distance_half/100, sigdigits=4)) m.w.e. to lose half its energy.")

energy = 1e6
distance = 1e3
E_after = upper_limit_track_integral(disp, energy, distance)
println("A $(energy/1e3) GeV muon loses $(round(energy - E_after, sigdigits=5)) MeV after $distance g/cm^2.")

# --- Lower limit ---
lower = get_lower_lim(disp)
println("\nDisplacement lower energy limit: $(lower) MeV")

free_crosssections(handle)
println("\nDone.")
