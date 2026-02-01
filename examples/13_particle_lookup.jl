# Particle Lookup and Definitions Example
# Demonstrates particle definitions, builder pattern, and type lookup.

using PROPOSAL

println("=== Particle Definitions ===\n")

# --- Predefined particles ---
println("--- Predefined particle definitions ---")
particles = [
    ("mu-", MuMinusDef()), ("mu+", MuPlusDef()),
    ("e-", EMinusDef()), ("e+", EPlusDef()),
    ("tau-", TauMinusDef()), ("tau+", TauPlusDef()),
    ("gamma", GammaDef()),
    ("pi0", Pi0Def()), ("pi+", PiPlusDef()), ("pi-", PiMinusDef()),
]
for (label, pdef) in particles
    println("  $label: mass=$(get_mass(pdef)) MeV, charge=$(get_charge(pdef)), lifetime=$(get_lifetime(pdef)), type=$(get_particle_type(pdef))")
end

# --- Lookup by particle type code ---
println("\n--- Lookup particle by PDG type code ---")
for code in [13, -13, 11, -11, 15, -15, 22]
    pdef = get_particle_def_for_type(code)
    println("  type=$code: $(get_name(pdef)), mass=$(get_mass(pdef)) MeV")
end

# --- Builder pattern ---
println("\n--- Custom particle via Builder ---")
builder = ParticleDefBuilder()
set_name(builder, "CustomMuon")
set_mass(builder, 105.0)
set_lifetime(builder, 2.2e-6)
set_charge(builder, -1.0)
set_particle_type(builder, 99)
custom = build(builder)
println("  Name: $(get_name(custom))")
println("  Mass: $(get_mass(custom)) MeV")
println("  Charge: $(get_charge(custom))")

# --- Constants ---
println("\n--- Physical constants ---")
println("  Speed of light: $SPEED_OF_LIGHT cm/ns")
println("  Electron mass: $ELECTRON_MASS MeV")
println("  Muon mass: $MUON_MASS MeV")
println("  Tau mass: $TAU_MASS MeV")
println("  Proton mass: $PROTON_MASS MeV")

println("\nDone.")
