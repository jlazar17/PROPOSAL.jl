# Medium and Component Details Example
# Demonstrates medium properties, component access, and custom media.

using PROPOSAL

println("=== Medium and Component Details ===\n")

# --- Standard media properties ---
println("--- Standard media ---")
for (name, create_fn) in [
    ("Water", create_water),
    ("Ice", create_ice),
    ("Standard Rock", create_standard_rock),
    ("Air", create_air),
    ("Iron", create_iron_medium),
    ("Lead", create_lead_medium),
    ("Liquid Argon", create_liquid_argon),
]
    m = create_fn()
    println("  $name:")
    println("    Mass density: $(get_mass_density(m)) g/cm³")
    println("    Mol density:  $(round(get_mol_density(m), sigdigits=5)) mol/cm³")
    println("    Radiation length: $(round(get_radiation_length(m), sigdigits=5)) g/cm²")
    println("    Num components: $(get_num_components(m))")
    println("    Sum Z/A: $(round(get_ZA(m), sigdigits=5))")
    println()
end

# --- Component details ---
println("--- Components of water ---")
water = create_water()
n_comp = get_components_size(water)
for i in 0:(n_comp - 1)
    comp = get_component(water, i)
    println("  Component $i: $(get_name(comp))")
    println("    Nuclear charge (Z): $(get_nuc_charge(comp))")
    println("    Atomic number (A):  $(get_atomic_num(comp))")
    println("    Atoms in molecule:  $(get_atom_in_molecule(comp))")
    println("    Average nucleon weight: $(round(get_average_nucleon_weight(comp), sigdigits=5)) MeV")
end

# --- Named components ---
println("\n--- Named components ---")
for (name, comp_fn) in [
    ("Hydrogen", ComponentHydrogen),
    ("Carbon", ComponentCarbon),
    ("Nitrogen", ComponentNitrogen),
    ("Oxygen", ComponentOxygen),
    ("Iron", ComponentIron),
    ("Lead", ComponentLead),
]
    comp = comp_fn()
    println("  $name: Z=$(get_nuc_charge(comp)), A=$(get_atomic_num(comp)), hash=$(get_hash(comp))")
end

# --- Component hash lookup ---
println("\n--- Component hash lookup ---")
h2 = ComponentHydrogen()
h = get_hash(h2)
println("  Hydrogen hash: $h")
looked_up = get_component_for_hash(h)
println("  Looked up by hash: $(get_name(looked_up)), Z=$(get_nuc_charge(looked_up))")

# --- Sternheimer parameters ---
println("\n--- Sternheimer parameters (for ionization loss) ---")
ice = create_ice()
println("  Ice:")
println("    I  (mean ionization energy): $(get_I(ice)) eV")
println("    C  (Sternheimer C): $(round(get_C(ice), sigdigits=5))")
println("    a  (Sternheimer a): $(round(get_A(ice), sigdigits=5))")
println("    m  (Sternheimer m): $(round(get_M(ice), sigdigits=5))")
println("    X0 (Sternheimer X0): $(round(get_X0(ice), sigdigits=5))")
println("    X1 (Sternheimer X1): $(round(get_X1(ice), sigdigits=5))")
println("    D0 (Sternheimer D0): $(round(get_D0(ice), sigdigits=5))")

# --- Comparison of water models ---
println("\n--- Water model comparison ---")
for (name, create_fn) in [
    ("Water (standard)", create_water),
    ("PDG2001 Water", create_pdg2001_water),
    ("PDG2020 Water", create_pdg2020_water),
    ("ANTARES Water", create_antares_water),
    ("Cascadia Basin Water", create_cascadia_basin_water),
]
    m = create_fn()
    println("  $name: ρ=$(get_mass_density(m)), X₀=$(round(get_radiation_length(m), sigdigits=6))")
end

println("\nDone.")
