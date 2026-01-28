/**
 * PROPOSAL Julia Bindings
 *
 * Minimal CxxWrap bindings for core PROPOSAL functionality.
 */

#include "jlcxx/jlcxx.hpp"
#include "jlcxx/array.hpp"

#include "PROPOSAL/PROPOSAL.h"

using namespace PROPOSAL;

JLCXX_MODULE define_julia_module(jlcxx::Module& mod) {

    // ========== Math Types ==========

    // Cartesian3D
    mod.add_type<Cartesian3D>("Cartesian3D")
        .constructor<>()
        .constructor<double, double, double>()
        .method("get_x", &Cartesian3D::GetX)
        .method("get_y", &Cartesian3D::GetY)
        .method("get_z", &Cartesian3D::GetZ)
        .method("set_x", &Cartesian3D::SetX)
        .method("set_y", &Cartesian3D::SetY)
        .method("set_z", &Cartesian3D::SetZ)
        .method("magnitude", &Cartesian3D::magnitude);

    // Spherical3D
    mod.add_type<Spherical3D>("Spherical3D")
        .constructor<>()
        .constructor<double, double, double>()
        .method("get_radius", &Spherical3D::GetRadius)
        .method("get_azimuth", &Spherical3D::GetAzimuth)
        .method("get_zenith", &Spherical3D::GetZenith)
        .method("set_radius", &Spherical3D::SetRadius)
        .method("set_azimuth", &Spherical3D::SetAzimuth)
        .method("set_zenith", &Spherical3D::SetZenith);

    // ========== Particle Types ==========

    // ParticleDef struct
    mod.add_type<ParticleDef>("ParticleDef")
        .method("get_name", [](const ParticleDef& p) { return std::string(p.name); })
        .method("get_mass", [](const ParticleDef& p) { return p.mass; })
        .method("get_lifetime", [](const ParticleDef& p) { return p.lifetime; })
        .method("get_charge", [](const ParticleDef& p) { return p.charge; })
        .method("get_particle_type", [](const ParticleDef& p) { return p.particle_type; });

    // Predefined particle definitions
    mod.add_type<MuMinusDef>("MuMinusDef", jlcxx::julia_base_type<ParticleDef>())
        .constructor<>();

    mod.add_type<MuPlusDef>("MuPlusDef", jlcxx::julia_base_type<ParticleDef>())
        .constructor<>();

    mod.add_type<EMinusDef>("EMinusDef", jlcxx::julia_base_type<ParticleDef>())
        .constructor<>();

    mod.add_type<EPlusDef>("EPlusDef", jlcxx::julia_base_type<ParticleDef>())
        .constructor<>();

    mod.add_type<TauMinusDef>("TauMinusDef", jlcxx::julia_base_type<ParticleDef>())
        .constructor<>();

    mod.add_type<TauPlusDef>("TauPlusDef", jlcxx::julia_base_type<ParticleDef>())
        .constructor<>();

    mod.add_type<GammaDef>("GammaDef", jlcxx::julia_base_type<ParticleDef>())
        .constructor<>();

    // ParticleState
    mod.add_type<ParticleState>("ParticleState")
        .constructor<>()
        .method("get_type", [](const ParticleState& p) { return static_cast<int>(p.type); })
        .method("get_position", [](const ParticleState& p) { return Cartesian3D(p.position); })
        .method("get_direction", [](const ParticleState& p) { return Cartesian3D(p.direction); })
        .method("get_energy", [](const ParticleState& p) { return p.energy; })
        .method("get_time", [](const ParticleState& p) { return p.time; })
        .method("get_propagated_distance", [](const ParticleState& p) { return p.propagated_distance; });

    // Factory function for ParticleState
    mod.method("make_particle_state", [](int type, double x, double y, double z,
                                         double dx, double dy, double dz,
                                         double energy, double time, double prop_dist) {
        ParticleState state;
        state.type = type;  // type is int in ParticleState
        state.position = Cartesian3D(x, y, z);
        state.direction = Cartesian3D(dx, dy, dz);
        state.energy = energy;
        state.time = time;
        state.propagated_distance = prop_dist;
        return state;
    });

    // ========== Energy Cut Settings ==========

    mod.add_type<EnergyCutSettings>("EnergyCutSettings")
        .constructor<double, double, bool>()
        .method("get_ecut", &EnergyCutSettings::GetEcut)
        .method("get_vcut", &EnergyCutSettings::GetVcut)
        .method("get_cont_rand", &EnergyCutSettings::GetContRand);

    // ========== Medium ==========

    mod.add_type<Medium>("Medium")
        .method("get_name", &Medium::GetName)
        .method("get_mass_density", &Medium::GetMassDensity)
        .method("get_mol_density", &Medium::GetMolDensity);

    // Medium factory function
    mod.method("create_medium", [](const std::string& name) {
        return std::shared_ptr<Medium>(CreateMedium(name));
    });

    // ========== Secondaries (Track) ==========

    mod.add_type<Secondaries>("Secondaries")
        .method("track_size", [](const Secondaries& s) { return s.GetTrack().size(); })
        .method("get_track_state", [](const Secondaries& s, size_t i) {
            return s.GetTrack().at(i);
        })
        .method("get_initial_state", &Secondaries::GetInitialState)
        .method("get_final_state", &Secondaries::GetFinalState)
        .method("get_track_length", &Secondaries::GetTrackLength)
        .method("get_track_energies_array", [](const Secondaries& s, jlcxx::ArrayRef<double> arr) {
            auto energies = s.GetTrackEnergies();
            for (size_t i = 0; i < energies.size() && i < arr.size(); ++i) {
                arr[i] = energies[i];
            }
            return energies.size();
        })
        .method("get_track_times_array", [](const Secondaries& s, jlcxx::ArrayRef<double> arr) {
            auto times = s.GetTrackTimes();
            for (size_t i = 0; i < times.size() && i < arr.size(); ++i) {
                arr[i] = times[i];
            }
            return times.size();
        })
        .method("get_track_propagated_distances_array", [](const Secondaries& s, jlcxx::ArrayRef<double> arr) {
            auto dists = s.GetTrackPropagatedDistances();
            for (size_t i = 0; i < dists.size() && i < arr.size(); ++i) {
                arr[i] = dists[i];
            }
            return dists.size();
        })
        .method("has_decay", [](const Secondaries& s) -> bool {
            auto types = s.GetTrackTypes();
            for (auto t : types) {
                if (t == InteractionType::Decay) return true;
            }
            return false;
        })
        .method("get_total_continuous_energy_loss", [](const Secondaries& s) -> double {
            auto losses = s.GetContinuousLosses();
            double total = 0.0;
            for (const auto& loss : losses) {
                total += loss.energy;
            }
            return total;
        })
        .method("get_decay_products_to_array", [](const Secondaries& s,
                    jlcxx::ArrayRef<int> types_arr,
                    jlcxx::ArrayRef<double> energies_arr,
                    jlcxx::ArrayRef<double> dx_arr,
                    jlcxx::ArrayRef<double> dy_arr,
                    jlcxx::ArrayRef<double> dz_arr) -> int {
            try {
                auto products = s.GetDecayProducts();
                int n = std::min(static_cast<int>(products.size()),
                                 static_cast<int>(types_arr.size()));
                for (int i = 0; i < n; ++i) {
                    types_arr[i] = static_cast<int>(products[i].type);
                    energies_arr[i] = products[i].energy;
                    auto dir = products[i].direction;
                    dx_arr[i] = dir.GetX();
                    dy_arr[i] = dir.GetY();
                    dz_arr[i] = dir.GetZ();
                }
                return static_cast<int>(products.size());
            } catch (...) {
                return 0;
            }
        });

    // ========== Propagator ==========

    mod.add_type<Propagator>("Propagator")
        .method("propagate", [](Propagator& prop, const ParticleState& initial,
                               double max_distance, double min_energy) {
            return prop.Propagate(initial, max_distance, min_energy);
        });

    // Factory functions for Propagator with each particle type
    mod.method("create_propagator_muminus", [](const std::string& config_path) {
        return Propagator(MuMinusDef(), config_path);
    });
    mod.method("create_propagator_muplus", [](const std::string& config_path) {
        return Propagator(MuPlusDef(), config_path);
    });
    mod.method("create_propagator_eminus", [](const std::string& config_path) {
        return Propagator(EMinusDef(), config_path);
    });
    mod.method("create_propagator_eplus", [](const std::string& config_path) {
        return Propagator(EPlusDef(), config_path);
    });
    mod.method("create_propagator_tauminus", [](const std::string& config_path) {
        return Propagator(TauMinusDef(), config_path);
    });
    mod.method("create_propagator_tauplus", [](const std::string& config_path) {
        return Propagator(TauPlusDef(), config_path);
    });
    mod.method("create_propagator_gamma", [](const std::string& config_path) {
        return Propagator(GammaDef(), config_path);
    });

    // ========== Random Number Generator ==========

    mod.method("set_random_seed", [](int seed) {
        RandomGenerator::Get().SetSeed(seed);
    });

    mod.method("random_double", []() {
        return RandomGenerator::Get().RandomDouble();
    });

    // ========== Version ==========

    mod.method("get_proposal_version", &getPROPOSALVersion);

    // ========== Constants ==========

    mod.set_const("SPEED_OF_LIGHT", SPEED);
    mod.set_const("ELECTRON_MASS", ME);
    mod.set_const("MUON_MASS", MMU);
    mod.set_const("TAU_MASS", MTAU);
    mod.set_const("PROTON_MASS", MP);

    // ========== Particle Type Constants ==========
    mod.set_const("PARTICLE_TYPE_NONE", static_cast<int>(ParticleType::None));
    mod.set_const("PARTICLE_TYPE_EMINUS", static_cast<int>(ParticleType::EMinus));
    mod.set_const("PARTICLE_TYPE_EPLUS", static_cast<int>(ParticleType::EPlus));
    mod.set_const("PARTICLE_TYPE_MUMINUS", static_cast<int>(ParticleType::MuMinus));
    mod.set_const("PARTICLE_TYPE_MUPLUS", static_cast<int>(ParticleType::MuPlus));
    mod.set_const("PARTICLE_TYPE_TAUMINUS", static_cast<int>(ParticleType::TauMinus));
    mod.set_const("PARTICLE_TYPE_TAUPLUS", static_cast<int>(ParticleType::TauPlus));
    mod.set_const("PARTICLE_TYPE_GAMMA", static_cast<int>(ParticleType::Gamma));
}
