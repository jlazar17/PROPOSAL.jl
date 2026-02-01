# Interaction Example
# Julia equivalent of Python Interaction.ipynb
#
# The interaction utility determines when and what type of stochastic
# interaction occurs, based on total cross sections.

using PROPOSAL

println("=== Interaction ===\n")

set_tables_path("/tmp")

pdef = MuMinusDef()
medium = create_ice()
handle = create_crosssections(pdef, medium, 500.0, 1.0, false, true)
inter = make_interaction(handle, false)

# --- Energy of next stochastic interaction ---
println("--- Energy of next stochastic interaction (E_i = 1e5 MeV) ---")
energy = 1e5
for rnd in [0.01, 0.1, 0.3, 0.5, 0.7, 0.9, 0.99]
    E_next = energy_interaction(inter, energy, rnd)
    println("  rnd=$(rnd): E_interaction = $(round(E_next, sigdigits=6)) MeV")
end

# --- Mean free path ---
println("\n--- Mean free path at different energies ---")
for E in [1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9]
    mfp = mean_free_path(inter, E)
    println("  E=$(E) MeV: lambda = $(round(mfp, sigdigits=5)) g/cm^2")
end

# --- Energy integral ---
println("\n--- Energy integral (interaction probability between E_i and E_f) ---")
E_i = 1e6
for E_f in [9e5, 5e5, 1e5, 1e4, 1e3]
    integral = energy_integral(inter, E_i, E_f)
    println("  E_i=$(E_i), E_f=$(E_f): integral = $(round(integral, sigdigits=5))")
end

# --- Interaction rates per component ---
println("\n--- Interaction rates at E=1e5 MeV ---")
energy = 1e5
rate_arr = zeros(Float64, 20)
hash_arr = zeros(Int64, 20)
n_rates = interaction_rates(inter[], energy, rate_arr, hash_arr)
println("Number of rate components: $n_rates")
for i in 1:n_rates
    println("  component hash=$(hash_arr[i]), rate=$(rate_arr[i])")
end

# --- Sample a loss ---
println("\n--- Sampling stochastic losses at E=1e5 MeV ---")
result = zeros(Float64, 3)
for rnd in [0.1, 0.3, 0.5, 0.7, 0.9]
    interaction_sample_loss(inter[], energy, rnd, rate_arr, hash_arr, n_rates, result)
    loss_type = Int(result[1])
    v_loss = result[3]
    println("  rnd=$(rnd): type=$(loss_type), v_loss=$(round(v_loss, sigdigits=5))")
end

free_crosssections(handle)
println("\nDone.")
