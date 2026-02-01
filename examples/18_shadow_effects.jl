# Shadow Effects Example
# Demonstrates nuclear shadow effects and their impact on photonuclear
# cross sections at high energies.
#
# Shadow effects account for the reduced effective nuclear cross section
# at small Bjorken-x (high energy), where partons overlap and screen
# each other. Two parametrizations are provided:
# - DuttaRenoSarcevicSeckel (DRSS)
# - ButkevichMikheyev (BM)

using PROPOSAL

println("=== Shadow Effects ===\n")

# --- Create shadow effect objects ---
shadow_drss = ShadowDuttaRenoSarcevicSeckel()
shadow_bm = ShadowButkevichMikheyev()

# --- Evaluate shadow factors for different kinematic regimes ---
# Shadow factor depends on: Component, Bjorken-x, nu (energy fraction)
# Use iron for more pronounced shadow effect (higher Z)
iron_medium = create_iron_medium()
iron_comp = get_component(iron_medium, 0)
println("  Using iron component: Z=$(get_nuc_charge(iron_comp)), A=$(get_atomic_num(iron_comp))")

# Also show oxygen for comparison
medium = create_ice()
oxygen = get_component(medium, 0)
println("  Using oxygen component: Z=$(get_nuc_charge(oxygen)), A=$(get_atomic_num(oxygen))")

println("\n--- Shadow factor vs Bjorken-x (nu=1e4) ---")
println("  x          Iron DRSS     Iron BM      O DRSS       O BM")
for log_x in range(-5, 0, length=11)
    x = 10.0^log_x
    sf_drss_fe = calculate_shadow_effect(shadow_drss, iron_comp, x, 1e4)
    sf_bm_fe = calculate_shadow_effect(shadow_bm, iron_comp, x, 1e4)
    sf_drss_o = calculate_shadow_effect(shadow_drss, oxygen, x, 1e4)
    sf_bm_o = calculate_shadow_effect(shadow_bm, oxygen, x, 1e4)
    println("  $(round(x, sigdigits=3))\t$(round(sf_drss_fe, sigdigits=5))\t$(round(sf_bm_fe, sigdigits=5))\t$(round(sf_drss_o, sigdigits=5))\t$(round(sf_bm_o, sigdigits=5))")
end

println("\n--- Shadow factor vs nu (iron, x=1e-3) ---")
println("  nu         DRSS          BM")
for nu in [1e1, 1e2, 1e3, 1e4, 1e5, 1e6]
    sf_drss = calculate_shadow_effect(shadow_drss, iron_comp, 1e-3, nu)
    sf_bm = calculate_shadow_effect(shadow_bm, iron_comp, 1e-3, nu)
    println("  $(nu)\t$(round(sf_drss, sigdigits=5))\t$(round(sf_bm, sigdigits=5))")
end

# --- Impact on photonuclear cross sections ---
println("\n--- PhotoQ2 parametrizations with different shadow effects ---")
pdef = MuMinusDef()
cuts = EnergyCutSettings(500.0, 0.05, false)

# ALLM97 with DRSS shadow
allm97_drss = create_photo_allm97_drss(shadow_drss)
# ALLM97 with BM shadow
allm97_bm = create_photo_allm97_bm(shadow_bm)

println("  Comparing ALLM97 with DRSS vs BM shadow:")
println("  E (MeV)        dσ/dv (DRSS)    dσ/dv (BM)")
for log_E in [4, 5, 6, 7, 8, 9]
    E = 10.0^log_E
    limits_drss = get_kinematic_limits(allm97_drss, pdef, oxygen, E)
    v_mid = 0.5 * (get_v_min(limits_drss) + get_v_max(limits_drss))
    if v_mid > 0
        dsigma_drss = differential_crosssection(allm97_drss, pdef, oxygen, E, v_mid)
        dsigma_bm = differential_crosssection(allm97_bm, pdef, oxygen, E, v_mid)
        println("  1e$log_E\t$(round(dsigma_drss, sigdigits=5))\t$(round(dsigma_bm, sigdigits=5))")
    end
end

# --- Available PhotoQ2 parametrizations ---
println("\n--- Available PhotoQ2 parametrizations with shadow ---")
println("  Each can be created with either DRSS or BM shadow:")
println("  - ALLM91:  create_photo_allm91_drss / create_photo_allm91_bm")
println("  - ALLM97:  create_photo_allm97_drss / create_photo_allm97_bm")
println("  - Butkevich: create_photo_butkevich_drss / create_photo_butkevich_bm")
println("  - Reno:    create_photo_reno_drss / create_photo_reno_bm")
println("  - AbtFT:   create_photo_abtft_drss / create_photo_abtft_bm")
println("  - BlockDurandHa: create_photo_blockdurandha_drss / create_photo_blockdurandha_bm")

println("\nDone.")
