# LPM Effect Example
# Demonstrates the Landau-Pomeranchuk-Migdal (LPM) effect and its
# suppression of bremsstrahlung, pair production, and photon pair
# production cross sections at ultra-high energies.
#
# == Physical Background ==
#
# The LPM effect arises because at very high energies, the longitudinal
# momentum transfer in bremsstrahlung and pair production becomes so
# small that the formation length of the process (the distance over
# which the radiation is coherently emitted) exceeds the mean free path
# for multiple scattering. When this happens, multiple scatterings
# during the formation zone disrupt the coherence, suppressing the
# cross section.
#
# The formation length scales as:
#   l_f ~ E / (m^2 * v * (1-v))   for bremsstrahlung
#   l_f ~ E / (m^2 * v)           for pair production
# where E is the particle energy and v is the fractional energy loss.
#
# The critical energy (LPM energy) above which suppression becomes
# significant is:
#   E_LPM ~ (m^2 * X_0 * alpha) / (4 * pi)
# where X_0 is the radiation length and alpha is the fine structure
# constant. For ice, E_LPM ~ 2e15 eV = 2 PeV.
#
# PROPOSAL implements three LPM correction objects:
# - BremsLPM:     suppression factor for bremsstrahlung
# - EpairLPM:     suppression factor for electron-positron pair production
# - PhotoPairLPM: suppression factor for photon-induced pair production
#
# Each provides a `suppression_factor(...)` method that returns a
# multiplicative correction to the differential cross section.
#
# == References ==
# - L.D. Landau, I.Ya. Pomeranchuk, Dokl. Akad. Nauk SSSR 92, 535 (1953)
# - A.B. Migdal, Phys. Rev. 103, 1811 (1956)
# - S. Klein, Rev. Mod. Phys. 71, 1501 (1999)
# - PROPOSAL: Koehne et al., Comp. Phys. Comm. 184, 2070 (2013)

using PROPOSAL

println("=== LPM (Landau-Pomeranchuk-Migdal) Effect ===\n")

pdef = MuMinusDef()
medium = create_ice()
iron_medium = create_iron_medium()
iron_comp = get_component(iron_medium, 0)

# =========================================================================
# 1. Bremsstrahlung LPM suppression
# =========================================================================
println("--- Bremsstrahlung LPM ---")
println("  The BremsLPM object computes the suppression factor for")
println("  bremsstrahlung as a function of (energy, v, component, density_correction).\n")

# Create a bremsstrahlung parametrization (needed as input)
brems = BremsKelnerKokoulinPetrukhin()
brems_lpm = create_brems_lpm(pdef, iron_medium, brems)
println("  Created BremsLPM from KelnerKokoulinPetrukhin parametrization")
println("  Hash: $(get_hash(brems_lpm))")

# The density_correction accounts for the actual density of the medium
# relative to the reference density. For homogeneous media at standard
# density, this is 1.0.
density_correction = 1.0

# density_correction = 1.0 for homogeneous media at standard density

# Show suppression factor vs energy at fixed v
println("\n  Suppression factor vs energy (v=0.5, iron component):")
println("  E (MeV)         suppression")
for log_E in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    E = 10.0^log_E
    sf = suppression_factor(brems_lpm, E, 0.5, iron_comp, density_correction)
    label = if log_E <= 5
        "  (no suppression)"
    elseif log_E >= 9
        "  (strong suppression)"
    else
        ""
    end
    println("  1e$(log_E)\t\t$(round(sf, sigdigits=5))$label")
end

# Show suppression factor vs v at fixed high energy
println("\n  Suppression factor vs v (E=1e10 MeV = 10 PeV, oxygen):")
println("  v              suppression")
for v in [0.01, 0.05, 0.1, 0.2, 0.3, 0.5, 0.7, 0.9, 0.95, 0.99]
    sf = suppression_factor(brems_lpm, 1e10, v, iron_comp, density_correction)
    println("  $(v)\t\t$(round(sf, sigdigits=5))")
end

println("\n  Physical interpretation:")
println("  - At low energies (E << E_LPM), suppression ≈ 1 (no effect)")
println("  - At high energies (E >> E_LPM), suppression < 1")
println("  - Suppression is strongest for soft radiation (small v)")
println("  - For hard radiation (v→1), the formation length is shorter")
println("    so the LPM effect is weaker")

# =========================================================================
# 2. Electron pair production LPM suppression
# =========================================================================
println("\n\n--- Electron Pair Production LPM ---")
println("  The EpairLPM object computes the suppression factor for e+e-")
println("  pair production. It takes more kinematic variables:\n")
println("  suppression_factor(E, v, r2, beta, xi, density_correction)")
println("    E     = parent particle energy")
println("    v     = fractional energy loss")
println("    r2    = (asymmetry parameter)^2, related to energy sharing")
println("    beta  = velocity parameter of the produced pair")
println("    xi    = screening parameter")

epair_lpm = create_epair_lpm(pdef, iron_medium)
println("\n  Created EpairLPM for muon in ice")
println("  Hash: $(get_hash(epair_lpm))")

# Use representative kinematic values
println("\n  Suppression factor vs energy (v=0.1, r2=0.25, beta=1.0, xi=1.0):")
println("  E (MeV)         suppression")
for log_E in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    E = 10.0^log_E
    sf = suppression_factor(epair_lpm, E, 0.1, 0.25, 1.0, 1.0, density_correction)
    println("  1e$(log_E)\t\t$(round(sf, sigdigits=5))")
end

# =========================================================================
# 3. Photon pair production LPM suppression
# =========================================================================
println("\n\n--- Photon Pair Production LPM ---")
println("  The PhotoPairLPM object computes the suppression factor for")
println("  gamma -> e+e- pair production (photo-pair production).\n")
println("  suppression_factor(E, x, component, density_correction)")
println("    E     = photon energy")
println("    x     = energy fraction of one member of the pair")

# Need a PhotoPairProduction parametrization
pp = PhotoPairTsai()
photopair_lpm = create_photopair_lpm(pdef, iron_medium, pp)
println("  Created PhotoPairLPM from Tsai parametrization")
println("  Hash: $(get_hash(photopair_lpm))")

println("\n  Suppression factor vs energy (x=0.5, oxygen):")
println("  E (MeV)         suppression")
for log_E in [3, 4, 5, 6, 7, 8, 9, 10, 11, 12]
    E = 10.0^log_E
    sf = suppression_factor(photopair_lpm, E, 0.5, iron_comp, density_correction)
    println("  1e$(log_E)\t\t$(round(sf, sigdigits=5))")
end

# =========================================================================
# 4. Impact on cross sections: with vs without LPM
# =========================================================================
println("\n\n--- Impact on total cross sections ---")
println("  Comparing dEdx for muon bremsstrahlung with and without LPM effect.")
println("  (Using standard cross section factories)\n")

set_tables_path("/tmp")

# Standard cross sections include LPM by default at high energies
handle_std = create_crosssections(pdef, medium, 500.0, 0.05, false, false)

println("  E (MeV)         dEdx_brems (MeV·cm²/g)")
for log_E in [3, 5, 7, 9, 11]
    E = 10.0^log_E
    # Sum dEdx over all cross section components
    n_xs = crosssection_count(handle_std)
    total_dedx = 0.0
    for i in 0:(n_xs - 1)
        xs = crosssection_at(handle_std, i)
        if get_interaction_type(xs) == INTERACTION_TYPE_BREMS
            total_dedx += calculate_dEdx(xs, E)
        end
    end
    println("  1e$(log_E)\t\t$(round(total_dedx, sigdigits=5))")
end

free_crosssections(handle_std)

# =========================================================================
# 5. Summary
# =========================================================================
println("\n\n--- Summary ---")
println("  The LPM effect suppresses electromagnetic radiation processes at")
println("  ultra-high energies (E >> E_LPM ≈ few PeV in ice):")
println("  ")
println("  Process              LPM Object    Key signature")
println("  ─────────────────────────────────────────────────────────")
println("  Bremsstrahlung       BremsLPM      Soft photon suppression")
println("  e+e- pair prod.      EpairLPM      Asymmetric pair suppression")
println("  γ→e+e- (photo-pair)  PhotoPairLPM  Pair creation suppression")
println("  ")
println("  In PROPOSAL, the suppression factors are applied multiplicatively")
println("  to the differential cross sections. The standard cross section")
println("  factories (create_crosssections) include LPM corrections")
println("  automatically when appropriate.")

println("\nDone.")
