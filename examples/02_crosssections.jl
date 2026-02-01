# Cross Section Example
# Julia equivalent of Python CrossSection.ipynb
#
# Demonstrates how to create and use cross sections, compute dEdx, dNdx,
# and sample stochastic losses.

using PROPOSAL

println("=== Cross Sections ===\n")

# --- Standard cross sections for different particles ---
println("--- Standard cross sections per particle ---")
medium = create_ice()
ecut = 500.0; vcut = 1.0; cont_rand = false; interpolate = false

for (name, pdef_fn) in [("mu-", MuMinusDef), ("e+", EPlusDef), ("gamma", GammaDef)]
    pdef = pdef_fn()
    handle = create_crosssections(pdef, medium, ecut, vcut, cont_rand, interpolate)
    n = crosssection_count(handle)
    println("Particle $name has $n standard cross section types:")
    for i in 0:(n-1)
        cs = crosssection_at(handle, i)
        println("  - type=$(get_interaction_type(cs)), name=$(get_parametrization_name(cs))")
    end
    free_crosssections(handle)
end

# --- dEdx calculation at different energies ---
println("\n--- dEdx for mu- in ice (e+e- pair production) ---")
pdef = MuMinusDef()
handle = create_crosssections(pdef, medium, 500.0, 0.05, false, false)
n_cross = crosssection_count(handle)

energies = [1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9]
for i in 0:(n_cross-1)
    cs = crosssection_at(handle, i)
    itype = get_interaction_type(cs)
    pname = get_parametrization_name(cs)
    println("\nInteraction type $itype ($pname):")
    for E in energies
        dedx = calculate_dEdx(cs, E)
        dndx = calculate_dNdx(cs, E)
        println("  E=$(E) MeV: dEdx=$(dedx), dNdx=$(dndx)")
    end
end

# --- Stochastic loss sampling ---
println("\n--- Stochastic loss sampling ---")
cs = crosssection_at(handle, 0)
energy = 1e6
dndx = calculate_dNdx(cs, energy)
println("At E=$(energy) MeV, dNdx = $(dndx)")
for rnd in [0.1, 0.3, 0.5, 0.7, 0.9]
    # get_hash returns the medium/component hash for sampling
    h = get_crosssection_hash(cs)
    v = calculate_stochastic_loss(cs, 0, energy, rnd * dndx)  # target=0 for first component
    println("  rnd=$(rnd) -> v_loss=$(v), E_loss=$(v * energy) MeV")
end

free_crosssections(handle)
println("\nDone.")
