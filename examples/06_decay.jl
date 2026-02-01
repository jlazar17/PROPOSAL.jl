# Decay Example
# Julia equivalent of Python Decay.ipynb
#
# The decay utility samples the energy at which a particle decays,
# given its initial energy and the medium density.

using PROPOSAL

println("=== Decay ===\n")

set_tables_path("/tmp")

pdef = MuMinusDef()
medium = create_standard_rock()
handle = create_crosssections(pdef, medium, 500.0, 1.0, false, true)
decay = make_decay_calc(handle, pdef, true)

# --- Decay energy sampling ---
println("--- Decay energy for muon at E=1e8 MeV in standard rock ---")
density = get_mass_density(medium)
E_init = 1e8
println("Initial energy: $(E_init) MeV")
println("Medium density: $(density) g/cm^3")

set_random_seed(42)
for i in 1:10
    rnd = random_double()
    E_decay = energy_decay(decay, E_init, rnd, density)
    println("  trial $i: rnd=$(round(rnd, sigdigits=4)), E_decay=$(round(E_decay, sigdigits=6)) MeV")
end

# --- Show that higher energy muons decay at higher energies ---
println("\n--- Decay energy vs initial energy ---")
rnd = 0.5
for E_init in [1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9, 1e12, 1e14]
    E_decay = energy_decay(decay, E_init, rnd, density)
    println("  E_init=$(E_init) MeV -> E_decay=$(round(E_decay, sigdigits=6)) MeV")
end

free_crosssections(handle)
println("\nDone.")
