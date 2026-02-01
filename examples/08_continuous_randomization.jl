# Continuous Randomization Example
# Julia equivalent of Python ContRand.ipynb
#
# Continuous randomization adds statistical fluctuations to continuous
# energy losses, improving the physical accuracy of the simulation.

using PROPOSAL

println("=== Continuous Randomization ===\n")

set_tables_path("/tmp")

pdef = MuMinusDef()
medium = create_standard_rock()
handle = create_crosssections(pdef, medium, 500.0, 1.0, true, true)
contrand = make_contrand(handle, true)

# --- Variance of continuous energy loss ---
println("--- Variance of continuous energy loss ---")
E_i = 1e5
for E_f in [0.999e5, 0.99e5, 0.9e5, 0.5e5]
    var = variance(contrand, E_i, E_f)
    std_dev = sqrt(var)
    println("  E_i=$(E_i), E_f=$(E_f): variance=$(round(var, sigdigits=5)), std=$(round(std_dev, sigdigits=5)) MeV")
end

# --- Randomized final energies ---
println("\n--- Randomized final energies (E_i=1e8, E_f=1e6, 20 samples) ---")
E_i = 1e8
E_f = 1e6
min_energy = get_mass(pdef)

set_random_seed(123)
energies = Float64[]
for i in 1:20
    rnd = random_double()
    E_randomized = energy_randomize(contrand, E_i, E_f, rnd, min_energy)
    push!(energies, E_randomized)
end

mean_E = sum(energies) / length(energies)
std_E = sqrt(sum((energies .- mean_E).^2) / (length(energies) - 1))
println("  Mean randomized energy: $(round(mean_E, sigdigits=6)) MeV (nominal E_f=$(E_f))")
println("  Std dev: $(round(std_E, sigdigits=5)) MeV")
println("  Min: $(round(minimum(energies), sigdigits=6)) MeV")
println("  Max: $(round(maximum(energies), sigdigits=6)) MeV")

free_crosssections(handle)
println("\nDone.")
