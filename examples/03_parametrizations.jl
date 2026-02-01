# Parametrization Example
# Julia equivalent of Python Parametrization.ipynb
#
# Shows how to use individual cross-section parametrizations to compute
# differential cross sections and kinematic limits.

using PROPOSAL

println("=== Parametrizations ===\n")

# Create a bremsstrahlung parametrization
param = BremsKelnerKokoulinPetrukhin(false)  # lpm=false
pdef = MuMinusDef()
comp = ComponentHydrogen(1.0)

# --- Kinematic limits at different energies ---
println("--- Kinematic limits for Brems (KKP) on Hydrogen ---")
energies = [1e3, 1e4, 1e5, 1e6, 1e7, 1e8, 1e9]
for E in energies
    limits = get_kinematic_limits(param, pdef, comp, E)
    println("  E=$(E) MeV: v_min=$(get_v_min(limits)), v_max=$(get_v_max(limits))")
end

# --- Differential cross section ---
println("\n--- Differential cross section at E=1e6 MeV ---")
E = 1e6
limits = get_kinematic_limits(param, pdef, comp, E)
vmin = get_v_min(limits)
vmax = get_v_max(limits)
for v in range(vmin, vmax, length=10)
    dsigma = differential_crosssection(param, pdef, comp, E, v)
    println("  v=$(round(v, sigdigits=4)): dsigma/dv = $(dsigma)")
end

# --- Lower energy limit ---
elim = get_lower_energy_lim_param(param, pdef)
println("\nLower energy limit: $(elim) MeV")

# --- Compare multiple bremsstrahlung parametrizations ---
println("\n--- dEdx comparison at E=1e6 MeV (different brems parametrizations) ---")
params = [
    ("KKP", BremsKelnerKokoulinPetrukhin(false)),
    ("PS",  BremsPetrukhinShestakov(false)),
    ("CS",  BremsCompleteScreening(false)),
    ("ABB", BremsAndreevBezrukovBugaev(false)),
    ("SSR", BremsSandrockSoedingreksoRhode(false)),
]

medium = create_ice()
E = 1e6
for (name, p) in params
    # Make a cross section from this parametrization
    handle = make_crosssection_brems_kkp(BremsKelnerKokoulinPetrukhin(false), pdef, medium, 500.0, 0.05, false, false)
    # Actually use the right parametrization
    if name == "KKP"
        handle = make_crosssection_brems_kkp(p, pdef, medium, 500.0, 0.05, false, false)
    elseif name == "PS"
        handle = make_crosssection_brems_ps(p, pdef, medium, 500.0, 0.05, false, false)
    elseif name == "CS"
        handle = make_crosssection_brems_cs(p, pdef, medium, 500.0, 0.05, false, false)
    elseif name == "ABB"
        handle = make_crosssection_brems_abb(p, pdef, medium, 500.0, 0.05, false, false)
    elseif name == "SSR"
        handle = make_crosssection_brems_ssr(p, pdef, medium, 500.0, 0.05, false, false)
    end
    cs = crosssection_at(handle, 0)
    dedx = calculate_dEdx(cs, E)
    println("  $name: dEdx = $(dedx)")
    free_crosssections(handle)
end

println("\nDone.")
