# Decay Channels Example
# Demonstrates DecayTable, LeptonicDecayChannel, TwoBodyPhaseSpace,
# and ManyBodyPhaseSpace.

using PROPOSAL

println("=== Decay Channels ===\n")

# --- Leptonic decay channel ---
# mu- -> e- + nu_mu_bar + nu_e (standard muon decay)
println("--- Leptonic decay channel (mu- -> e- + nu_mu_bar + nu_e) ---")
muon = MuMinusDef()
electron = EMinusDef()
nu_mu = NuMuDef()
nu_e = NuEDef()

ldc = LeptonicDecayChannel(nu_mu, electron, nu_e)
ldc_approx = LeptonicDecayChannelApprox(nu_mu, electron, nu_e)
println("  Channel name: $(get_channel_name(ldc))")

# Decay a muon with some kinetic energy
state_rest = ParticleState(PARTICLE_TYPE_MUMINUS, 0, 0, 0, 0, 0, 1, 1000.0)

# Decay products stored in arrays
n_max = 10
types_arr = zeros(Int32, n_max)
energies_arr = zeros(Float64, n_max)
dirs_x = zeros(Float64, n_max)
dirs_y = zeros(Float64, n_max)
dirs_z = zeros(Float64, n_max)

set_random_seed(42)
# Note: LeptonicDecayChannel uses NewtonRaphson which may not converge for all seeds.
# Use LeptonicDecayChannelApprox for robust results.
n_products = decay_channel_decay_to_arrays(ldc_approx, muon, state_rest, types_arr, energies_arr, dirs_x, dirs_y, dirs_z)
println("  Number of decay products: $n_products")
for i in 1:n_products
    println("  Product $i: type=$(types_arr[i]), E=$(round(energies_arr[i], sigdigits=5)) MeV, dir=($(round(dirs_x[i], sigdigits=4)), $(round(dirs_y[i], sigdigits=4)), $(round(dirs_z[i], sigdigits=4)))")
end
total_E = sum(energies_arr[i] for i in 1:n_products)
println("  Total energy: $(round(total_E, sigdigits=6)) MeV (parent energy = 1000 MeV)")

# --- Approximate leptonic decay channel ---
println("\n--- Approximate leptonic decay (faster, less precise) ---")
println("  Channel name: $(get_channel_name(ldc_approx))")

# --- Two-body phase space ---
println("\n--- Two-body phase space (pi0 -> gamma + gamma) ---")
pi0 = Pi0Def()
gamma = GammaDef()
two_body = TwoBodyPhaseSpace(gamma, gamma)
println("  Channel name: $(get_channel_name(two_body))")

state_pi0 = ParticleState(PARTICLE_TYPE_PI0, 0, 0, 0, 0, 0, 1, get_mass(pi0))
set_random_seed(42)
n_products = decay_channel_decay_to_arrays(two_body, pi0, state_pi0, types_arr, energies_arr, dirs_x, dirs_y, dirs_z)
println("  Number of decay products: $n_products")
for i in 1:n_products
    println("  Product $i: type=$(types_arr[i]), E=$(round(energies_arr[i], sigdigits=5)) MeV")
end

# --- Many-body phase space ---
println("\n--- Many-body phase space ---")
mb3 = create_many_body_phase_space_3(electron, nu_mu, nu_e)
println("  Channel name: $(get_channel_name(mb3))")

# --- Decay table ---
println("\n--- Decay table ---")
dt = DecayTable()
add_channel(dt, 1.0, ldc)  # 100% to leptonic channel
println("  Added leptonic channel with BR=1.0")

# Select a channel
set_random_seed(42)
selected = select_channel(dt, random_double())
println("  Selected channel: $(get_channel_name(selected))")

# --- Custom particle with decay table ---
println("\n--- Custom particle with modified decay table ---")
builder = ParticleDefBuilder()
set_particle_def(builder, muon)
set_name(builder, "CustomMuon")
set_particle_type(builder, 9999)  # unique type ID

# Create a new decay table: 50% leptonic, 50% approximate leptonic
dt2 = DecayTable()
add_channel(dt2, 0.5, ldc)
add_channel(dt2, 0.5, ldc_approx)
set_decay_table(builder, dt2)

custom_muon = build(builder)
println("  Custom muon name: $(get_name(custom_muon))")
println("  Custom muon mass: $(get_mass(custom_muon)) MeV")

# --- Stable channel ---
println("\n--- Stable channel ---")
stable = StableChannel()
println("  Stable channel name: $(get_channel_name(stable))")

println("\nDone.")
