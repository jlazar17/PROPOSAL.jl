# Energy Cut Settings Example
# Julia equivalent of Python EnergyCut.ipynb
#
# Energy cuts control the separation between continuous and stochastic energy losses.
# - e_cut: absolute energy threshold in MeV
# - v_cut: relative energy threshold (fraction of particle energy)
# - cont_rand: whether to apply continuous randomization

using PROPOSAL

println("=== Energy Cut Settings ===\n")

# Create different cut configurations
cut_ecut_only = EnergyCutSettings(500.0, 1.0, false)   # only absolute cut at 500 MeV
cut_vcut_only = EnergyCutSettings(Inf, 0.05, false)     # only relative cut at 5%
cut_both = EnergyCutSettings(500.0, 0.05, false)        # both cuts active
cut_stochastic_only = EnergyCutSettings(Inf, 1.0, false) # all losses are continuous

println("e_cut only:       ecut=$(get_ecut(cut_ecut_only)), vcut=$(get_vcut(cut_ecut_only))")
println("v_cut only:       ecut=$(get_ecut(cut_vcut_only)), vcut=$(get_vcut(cut_vcut_only))")
println("both:             ecut=$(get_ecut(cut_both)), vcut=$(get_vcut(cut_both))")
println("stochastic only:  ecut=$(get_ecut(cut_stochastic_only)), vcut=$(get_vcut(cut_stochastic_only))")

# Demonstrate how the effective cut changes with energy
println("\nEffective cut at different energies (with ecut=500, vcut=0.01):")
cut = EnergyCutSettings(500.0, 0.01, false)
for energy in [1e3, 1e4, 1e5, 1e6, 1e9]
    c = PROPOSAL.cut(cut, energy)
    println("  E = $(energy) MeV -> cut = $(c)")
end
