module PROPOSAL

using CxxWrap
using Libdl

# Try to load the JLL package if available, otherwise use local development path
const _using_jll = try
    @eval using PROPOSAL_cxxwrap_jll
    true
catch
    false
end

function get_lib_path()
    if _using_jll
        return PROPOSAL_cxxwrap_jll.libPROPOSAL_cxxwrap_path
    else
        # Development fallback - check local build
        if Sys.isapple()
            ext = "dylib"
        elseif Sys.iswindows()
            ext = "dll"
        else
            ext = "so"
        end

        dev_libpath = joinpath(@__DIR__, "..", "build", "lib")
        dev_lib = joinpath(dev_libpath, "libPROPOSAL_cxxwrap.$ext")

        if isfile(dev_lib)
            return dev_lib
        end

        # Check environment variable
        env_path = get(ENV, "PROPOSAL_JL_LIB_PATH", "")
        if !isempty(env_path)
            env_lib = joinpath(env_path, "libPROPOSAL_cxxwrap.$ext")
            if isfile(env_lib)
                return env_lib
            end
        end

        return nothing
    end
end

const _lib_path = get_lib_path()
const _library_available = _lib_path !== nothing

if _library_available
    @wrapmodule(() -> _lib_path, :define_julia_module, Libdl.RTLD_LAZY | Libdl.RTLD_GLOBAL)

    function __init__()
        @initcxx
    end
else
    function __init__()
        @warn """
        PROPOSAL wrapper library not found.

        To use PROPOSAL.jl, either:
        1. Install PROPOSAL_cxxwrap_jll (when available in registry)
        2. Build locally: see deps/binarybuilder/README.md
        3. Set PROPOSAL_JL_LIB_PATH environment variable

        The package structure is ready but bindings are not loaded.
        """
    end
end

# Check if library is available
is_library_available() = _library_available
export is_library_available

# Exports - these will be available when the library is loaded
if _library_available
    # Math types
    export Vector3D, Cartesian3D, Spherical3D, UnitSphericalVector
    export get_x, get_y, get_z, set_x, set_y, set_z, magnitude, normalize, deflect
    export get_radius, get_azimuth, get_zenith, set_radius, set_azimuth, set_zenith
    export cartesian3d_dot

    # ParticleState accessors
    export get_type, get_position, get_direction, get_energy, get_time, get_propagated_distance

    # Particle definitions
    export ParticleDef, ParticleState
    export get_mass, get_lifetime, get_charge, get_particle_type
    export MuMinusDef, MuPlusDef, EMinusDef, EPlusDef
    export TauMinusDef, TauPlusDef, GammaDef
    export Pi0Def, PiMinusDef, PiPlusDef
    export K0Def, KMinusDef, KPlusDef
    export NuEDef, NuEBarDef, NuMuDef, NuMuBarDef, NuTauDef, NuTauBarDef
    export StauMinusDef, StauPlusDef
    export MonopoleDef, SMPMinusDef, SMPPlusDef

    # ParticleDef::Builder
    export ParticleDefBuilder
    export set_name, set_mass, set_low, set_lifetime, set_charge
    export set_decay_table, set_particle_type, set_particle_def, build

    # Component
    export Component
    export get_name, get_nuc_charge, get_atomic_num, get_atom_in_molecule
    export get_log_constant, get_b_prime, get_average_nucleon_weight, get_wood_saxon
    export get_hash, get_component_for_hash

    # Named Components
    export ComponentHydrogen, ComponentCarbon, ComponentNitrogen, ComponentOxygen
    export ComponentSodium, ComponentMagnesium, ComponentSulfur, ComponentChlorine
    export ComponentArgon, ComponentPotassium, ComponentCalcium, ComponentIron
    export ComponentCopper, ComponentLead, ComponentUranium
    export ComponentStandardRock, ComponentFrejusRock

    # Propagation
    export EnergyCutSettings, get_ecut, get_vcut, get_cont_rand
    export Medium, Secondaries, Propagator
    export track_size, get_track_state, get_track_length
    export get_initial_state, get_final_state
    export get_track_energies_array, get_track_times_array, get_track_propagated_distances_array
    export has_decay, get_decay_products_to_array, get_total_continuous_energy_loss
    export cut

    # Medium extended
    export get_mass_density, get_mol_density
    export get_I, get_C, get_A, get_M, get_X0, get_X1, get_D0
    export get_radiation_length, get_MM, get_sum_charge, get_ZA
    export get_num_components, get_sum_nucleons
    export get_component, get_components_size

    # Named Media (factory functions)
    export create_water, create_ice, create_salt, create_calcium_carbonate
    export create_standard_rock, create_frejus_rock
    export create_iron_medium, create_hydrogen_medium, create_lead_medium
    export create_copper_medium, create_uranium_medium
    export create_paraffin, create_air, create_liquid_argon
    export create_antares_water, create_cascadia_basin_water
    export create_pdg2001_water, create_pdg2001_ice, create_pdg2020_water, create_pdg2020_ice

    # Geometry
    export Geometry, Sphere, Cylinder, Box
    export is_inside, is_infront, is_behind, is_entering, is_leaving
    export distance_to_border_first, distance_to_border_second
    export distance_to_closest_approach
    export get_geometry_name, get_hierarchy, set_hierarchy
    export get_inner_radius

    # Density distributions
    export Axis, CartesianAxis, RadialAxis
    export DensityDistribution, DensityHomogeneous, DensityExponential
    export DensityPolynomial, DensitySplines
    export evaluate, integrate, calculate, correct

    # Cross sections
    export CrossSectionBase
    export calculate_dEdx, calculate_dE2dx, calculate_dNdx, calculate_stochastic_loss
    export calculate_cumulative_crosssection, calculate_dNdx_per_target
    export get_interaction_type, get_parametrization_name, get_lower_energy_lim
    export get_crosssection_hash
    export create_crosssections, crosssection_count, crosssection_at, free_crosssections

    # Decay
    export DecayChannel, StableChannel, LeptonicDecayChannelApprox, LeptonicDecayChannel
    export TwoBodyPhaseSpace
    export ManyBodyPhaseSpace, create_many_body_phase_space_2, create_many_body_phase_space_3
    export DecayTable, add_channel, set_stable, set_uniform_sampling, get_channel_name
    export select_channel, decay_channel_decay_to_arrays

    # Secondaries extended
    export StochasticLoss, ContinuousLoss
    export get_parent_particle_energy, get_target_hash
    export get_start_position, get_end_position
    export get_direction_initial, get_direction_final
    export get_time_initial, get_time_final
    export get_state_for_energy, get_state_for_distance
    export get_stochastic_losses_count, get_stochastic_loss_at
    export get_continuous_losses_count, get_continuous_loss_at
    export get_track_positions_count, get_track_position_at
    export get_track_directions_count, get_track_direction_at
    export get_track_types_array, get_target_hashes_array
    export hit_geometry, get_elost
    export get_entry_point, has_entry_point, get_exit_point, has_exit_point
    export get_stochastic_losses_in_geometry_count, get_stochastic_loss_in_geometry_at
    export get_track_in_geometry_count, get_track_in_geometry_at

    # Scattering
    export ScatteringOffset, get_sx, get_sy, get_tx, get_ty
    export MultipleScattering, Highland, HighlandIntegral, Moliere, MoliereInterpol
    export scatter, scattering_angle, scattering_angle_2d, calculate_theta0
    export create_highland_integral
    export StochasticDeflection, required_random_numbers
    export make_stochastic_deflection
    export Scattering
    export scattering_multiple_scatter, scattering_stochastic_deflection
    export scattering_ms_random_numbers, scattering_sd_random_numbers
    export create_scattering_ms_only, create_scattering_with_sd
    export create_scattering_by_types, create_scattering_multiplier

    # Propagation Utility
    export Displacement, Interaction, ContRand, DecayCalc, TimeCalc
    export solve_track_integral, upper_limit_track_integral, get_lower_lim
    export energy_interaction, energy_integral, mean_free_path
    export interaction_rates, interaction_sample_loss
    export variance, energy_randomize
    export energy_decay
    export time_elapsed
    export PropagationUtilityCollection, PropagationUtility
    export set_displacement!, set_interaction!, set_scattering!, set_decay!, set_contrand!, set_time!
    export energy_stochasticloss, energy_distance, length_continuous
    export create_propagation_utility
    export make_displacement, make_interaction, make_contrand
    export make_decay_calc, make_time_calc, make_time_approximate
    export make_propagator_single_sector_sphere, make_propagator_single_sector_cylinder
    export make_propagator_single_sector_box
    export make_propagator_single_sector_sphere_exp, make_propagator_single_sector_cylinder_exp
    export make_propagator_single_sector_box_exp
    export clear_sectors, add_sector_sphere_homogeneous, add_sector_cylinder_homogeneous
    export add_sector_box_homogeneous, add_sector_sphere_exponential
    export make_propagator_from_sectors
    export propagate_with_hierarchy

    # Crosssection Parametrizations
    export KinematicLimits, get_v_min, get_v_max
    export ParametrizationForComponent, ParametrizationForMedium, ParametrizationDirect
    export differential_crosssection, get_kinematic_limits, get_lower_energy_lim_param

    # Shadow effects
    export ShadowEffect, ShadowDuttaRenoSarcevicSeckel, ShadowButkevichMikheyev
    export calculate_shadow_effect

    # Bremsstrahlung parametrizations
    export Bremsstrahlung
    export BremsKelnerKokoulinPetrukhin, BremsPetrukhinShestakov
    export BremsCompleteScreening, BremsAndreevBezrukovBugaev
    export BremsSandrockSoedingreksoRhode, BremsElectronScreening

    # EpairProduction parametrizations
    export EpairProduction, EpairProductionRhoIntegral
    export EpairKelnerKokoulinPetrukhin, EpairSandrockSoedingreksoRhode, EpairForElectronPositron

    # Photonuclear parametrizations
    export Photonuclear, PhotoRealPhotonAssumption, PhotoQ2Integral
    export PhotoZeus, PhotoBezrukovBugaev, PhotoKokoulin, PhotoRhode
    export PhotoAbramowiczLevinLevyMaor91, PhotoAbramowiczLevinLevyMaor97
    export PhotoButkevichMikheyev, PhotoRenoSarcevicSu, PhotoAbtFT, PhotoBlockDurandHa
    export create_photo_allm91_drss, create_photo_allm91_bm
    export create_photo_allm97_drss, create_photo_allm97_bm
    export create_photo_butkevich_drss, create_photo_butkevich_bm
    export create_photo_reno_drss, create_photo_reno_bm
    export create_photo_abtft_drss, create_photo_abtft_bm
    export create_photo_blockdurandha_drss, create_photo_blockdurandha_bm

    # MupairProduction parametrizations
    export MupairProduction, MupairProductionRhoIntegral, MupairKelnerKokoulinPetrukhin

    # WeakInteraction parametrizations
    export WeakInteraction, WeakCooperSarkarMertsch

    # Compton parametrizations
    export Compton, ComptonKleinNishina

    # PhotoPairProduction parametrizations
    export PhotoPairProduction, PhotoPairTsai, PhotoPairKochMotz

    # PhotoMuPairProduction parametrizations
    export PhotoMuPairProduction, PhotoMuPairBurkhardtKelnerKokoulin

    # Ionization parametrizations
    export Ionization, IonizBetheBlochRossi, IonizBergerSeltzerBhabha, IonizBergerSeltzerMoller

    # Annihilation parametrizations
    export Annihilation, AnnihilationHeitler

    # Photoproduction parametrizations
    export Photoproduction
    export PhotoproductionZeus, PhotoproductionBezrukovBugaev, PhotoproductionCaldwell
    export PhotoproductionKokoulin, PhotoproductionRhode
    export PhotoproductionHeck, PhotoproductionHeckC7Shadowing

    # Photoeffect parametrizations
    export Photoeffect, PhotoeffectSauter

    # LPM correction objects
    export BremsLPM, create_brems_lpm
    export EpairLPM, create_epair_lpm
    export PhotoPairLPM, create_photopair_lpm
    export suppression_factor

    # Per-particle standard cross-section factories
    export make_std_crosssection_muminus, make_std_crosssection_muplus
    export make_std_crosssection_eminus, make_std_crosssection_eplus
    export make_std_crosssection_tauminus, make_std_crosssection_tauplus
    export make_std_crosssection_gamma

    # Cross-section from individual parametrizations (returns handle)
    export make_crosssection_brems_kkp, make_crosssection_brems_ps
    export make_crosssection_brems_cs, make_crosssection_brems_abb
    export make_crosssection_brems_ssr, make_crosssection_brems_es
    export make_crosssection_epair_kkp, make_crosssection_epair_ssr, make_crosssection_epair_fep
    export make_crosssection_photo_zeus, make_crosssection_photo_bb
    export make_crosssection_photo_kokoulin, make_crosssection_photo_rhode
    export make_crosssection_photo_allm91, make_crosssection_photo_allm97
    export make_crosssection_photo_bm, make_crosssection_photo_rss
    export make_crosssection_photo_abtft, make_crosssection_photo_bdh
    export make_crosssection_mupair_kkp, make_crosssection_weak_csm
    export make_crosssection_compton_kn
    export make_crosssection_photopair_tsai, make_crosssection_photopair_km
    export make_crosssection_photomupair_bkk
    export make_crosssection_ioniz_bbr, make_crosssection_ioniz_bsb, make_crosssection_ioniz_bsm
    export make_crosssection_annihilation_heitler
    export make_crosssection_photoproduction_zeus, make_crosssection_photoproduction_bb
    export make_crosssection_photoproduction_caldwell, make_crosssection_photoproduction_kokoulin
    export make_crosssection_photoproduction_rhode, make_crosssection_photoproduction_heck
    export make_crosssection_photoproduction_heckc7
    export make_crosssection_photoeffect_sauter

    # InterpolationSettings
    export get_tables_path, set_tables_path
    export get_upper_energy_lim, set_upper_energy_lim
    export get_nodes_dedx, set_nodes_dedx, get_nodes_de2dx, set_nodes_de2dx
    export get_nodes_dndx_e, set_nodes_dndx_e, get_nodes_dndx_v, set_nodes_dndx_v
    export get_nodes_utility, set_nodes_utility
    export get_nodes_rate_interpolant, set_nodes_rate_interpolant

    # PropagationSettings
    export get_max_steps, set_max_steps

    # SecondariesCalculator (handle-based)
    export create_secondaries_calculator, free_secondaries_calculator
    export sec_calc_required_random_numbers, sec_calc_calculate, sec_calc_result_at

    # Math types
    export Polynom, Spline, LinearSpline, CubicSpline
    export get_coefficient, derive, antiderivative

    # Particle lookup
    export get_particle_def_for_type

    # Functions
    export make_particle_state, create_medium, create_propagator
    export set_random_seed, random_double, get_proposal_version
    export create_propagator_muminus, create_propagator_muplus
    export create_propagator_eminus, create_propagator_eplus
    export create_propagator_tauminus, create_propagator_tauplus
    export create_propagator_gamma

    # Logging
    export set_loglevel
    export LOG_TRACE, LOG_DEBUG, LOG_INFO, LOG_WARN, LOG_ERROR, LOG_CRITICAL, LOG_OFF

    # Constants
    export SPEED_OF_LIGHT, ELECTRON_MASS, MUON_MASS, TAU_MASS, PROTON_MASS

    # Particle type constants
    export PARTICLE_TYPE_NONE, PARTICLE_TYPE_EMINUS, PARTICLE_TYPE_EPLUS
    export PARTICLE_TYPE_MUMINUS, PARTICLE_TYPE_MUPLUS
    export PARTICLE_TYPE_TAUMINUS, PARTICLE_TYPE_TAUPLUS, PARTICLE_TYPE_GAMMA
    export PARTICLE_TYPE_PI0, PARTICLE_TYPE_PIPLUS, PARTICLE_TYPE_PIMINUS
    export PARTICLE_TYPE_K0, PARTICLE_TYPE_KPLUS, PARTICLE_TYPE_KMINUS
    export PARTICLE_TYPE_NUE, PARTICLE_TYPE_NUEBAR
    export PARTICLE_TYPE_NUMU, PARTICLE_TYPE_NUMUBAR
    export PARTICLE_TYPE_NUTAU, PARTICLE_TYPE_NUTAUBAR
    export PARTICLE_TYPE_STAUMINUS, PARTICLE_TYPE_STAUPLUS
    export PARTICLE_TYPE_MONOPOLE, PARTICLE_TYPE_SMPPLUS, PARTICLE_TYPE_SMPMINUS
    export PARTICLE_TYPE_HADRON

    # InteractionType constants
    export INTERACTION_TYPE_UNDEFINED, INTERACTION_TYPE_PARTICLE
    export INTERACTION_TYPE_BREMS, INTERACTION_TYPE_IONIZ
    export INTERACTION_TYPE_EPAIR, INTERACTION_TYPE_PHOTONUCLEAR
    export INTERACTION_TYPE_MUPAIR, INTERACTION_TYPE_HADRONS
    export INTERACTION_TYPE_CONTINUOUS_ENERGY_LOSS, INTERACTION_TYPE_WEAKINT
    export INTERACTION_TYPE_COMPTON, INTERACTION_TYPE_DECAY
    export INTERACTION_TYPE_ANNIHILATION, INTERACTION_TYPE_PHOTOPAIR
    export INTERACTION_TYPE_PHOTOPRODUCTION, INTERACTION_TYPE_PHOTOMUPAIR
    export INTERACTION_TYPE_PHOTOEFFECT
end

# High-level Julia API (available when library is loaded)
if _library_available

    # Cartesian3D operator overloads
    Base.:+(a::Cartesian3D, b::Cartesian3D) = cartesian3d_add(a, b)
    Base.:-(a::Cartesian3D, b::Cartesian3D) = cartesian3d_subtract(a, b)
    Base.:-(a::Cartesian3D) = cartesian3d_negate(a)
    Base.:*(a::Cartesian3D, s::Real) = cartesian3d_scale(a, Float64(s))
    Base.:*(s::Real, a::Cartesian3D) = cartesian3d_scale(a, Float64(s))

    # Auto-dereference shared_ptr for utility method calls
    const SharedPtr = CxxWrap.StdLib.SharedPtrAllocated
    mean_free_path(x::SharedPtr{Interaction}, e) = mean_free_path(x[], e)
    energy_interaction(x::SharedPtr{Interaction}, e, r) = energy_interaction(x[], e, r)
    energy_integral(x::SharedPtr{Interaction}, ei, ef) = energy_integral(x[], ei, ef)
    solve_track_integral(x::SharedPtr{Displacement}, u, l) = solve_track_integral(x[], u, l)
    upper_limit_track_integral(x::SharedPtr{Displacement}, e, d) = upper_limit_track_integral(x[], e, d)
    get_lower_lim(x::SharedPtr{Displacement}) = get_lower_lim(x[])
    variance(x::SharedPtr{ContRand}, ei, ef) = variance(x[], ei, ef)
    energy_randomize(x::SharedPtr{ContRand}, ei, ef, r, m) = energy_randomize(x[], ei, ef, r, m)
    energy_decay(x::SharedPtr{DecayCalc}, e, r, d) = energy_decay(x[], e, r, d)
    time_elapsed(x::SharedPtr{TimeCalc}, ei, ef, g, d) = time_elapsed(x[], ei, ef, g, d)

    # Convenience dispatchers for sector-based propagator
    function make_propagator_single_sector(pdef, geo::Sphere, dens::DensityHomogeneous, coll)
        make_propagator_single_sector_sphere(pdef, geo, dens, coll)
    end
    function make_propagator_single_sector(pdef, geo::Cylinder, dens::DensityHomogeneous, coll)
        make_propagator_single_sector_cylinder(pdef, geo, dens, coll)
    end
    function make_propagator_single_sector(pdef, geo::Box, dens::DensityHomogeneous, coll)
        make_propagator_single_sector_box(pdef, geo, dens, coll)
    end
    function make_propagator_single_sector(pdef, geo::Sphere, dens::DensityExponential, coll)
        make_propagator_single_sector_sphere_exp(pdef, geo, dens, coll)
    end
    function make_propagator_single_sector(pdef, geo::Cylinder, dens::DensityExponential, coll)
        make_propagator_single_sector_cylinder_exp(pdef, geo, dens, coll)
    end
    function make_propagator_single_sector(pdef, geo::Box, dens::DensityExponential, coll)
        make_propagator_single_sector_box_exp(pdef, geo, dens, coll)
    end
    export make_propagator_single_sector

    """
        propagate(propagator, state; max_distance=1e20, min_energy=0.0)

    Propagate a particle through the configured medium.
    """
    function propagate(prop::Propagator, state::ParticleState;
                       max_distance::Real=1e20, min_energy::Real=0.0)
        return PROPOSAL.propagate(prop, state, Float64(max_distance), Float64(min_energy))
    end
    export propagate

    """
        ParticleState(particle_type, x, y, z, dx, dy, dz, energy; time=0.0, propagated_distance=0.0)

    Create a new particle state.
    """
    function ParticleState(particle_type::Integer,
                          x::Real, y::Real, z::Real,
                          dx::Real, dy::Real, dz::Real,
                          energy::Real;
                          time::Real=0.0,
                          propagated_distance::Real=0.0)
        return make_particle_state(Int(particle_type),
                                   Float64(x), Float64(y), Float64(z),
                                   Float64(dx), Float64(dy), Float64(dz),
                                   Float64(energy), Float64(time),
                                   Float64(propagated_distance))
    end

    """
        get_components(medium::Medium)

    Return a Vector of Component objects from a Medium.
    """
    function get_components(medium)
        n = get_components_size(medium)
        return [get_component(medium, i) for i in 0:(n-1)]
    end
    export get_components

    """
        get_stochastic_losses(secondaries::Secondaries)

    Return a Vector of StochasticLoss objects.
    """
    function get_stochastic_losses(sec)
        n = get_stochastic_losses_count(sec)
        return [get_stochastic_loss_at(sec, i) for i in 0:(n-1)]
    end
    export get_stochastic_losses

    """
        get_continuous_losses(secondaries::Secondaries)

    Return a Vector of ContinuousLoss objects.
    """
    function get_continuous_losses(sec)
        n = get_continuous_losses_count(sec)
        return [get_continuous_loss_at(sec, i) for i in 0:(n-1)]
    end
    export get_continuous_losses

    """
        get_track_positions(secondaries::Secondaries)

    Return a Vector of Cartesian3D positions along the track.
    """
    function get_track_positions(sec)
        n = get_track_positions_count(sec)
        return [get_track_position_at(sec, i) for i in 0:(n-1)]
    end
    export get_track_positions

    """
        get_track_directions(secondaries::Secondaries)

    Return a Vector of Cartesian3D directions along the track.
    """
    function get_track_directions(sec)
        n = get_track_directions_count(sec)
        return [get_track_direction_at(sec, i) for i in 0:(n-1)]
    end
    export get_track_directions

    """
        distance_to_border(geometry, position, direction)

    Return (first, second) distances to the geometry border.
    """
    # Auto-dereference shared_ptr for Scattering
    scattering_multiple_scatter(x::SharedPtr{Scattering}, args...) = scattering_multiple_scatter(x[], args...)
    scattering_stochastic_deflection(x::SharedPtr{Scattering}, args...) = scattering_stochastic_deflection(x[], args...)
    scattering_ms_random_numbers(x::SharedPtr{Scattering}) = scattering_ms_random_numbers(x[])
    scattering_sd_random_numbers(x::SharedPtr{Scattering}, args...) = scattering_sd_random_numbers(x[], args...)

    # Convenience constructors that accept Julia Vector{Float64}
    Polynom(v::Vector{Float64}) = Polynom(CxxWrap.StdVector(v))
    LinearSpline(x::Vector{Float64}, y::Vector{Float64}) = LinearSpline(CxxWrap.StdVector(x), CxxWrap.StdVector(y))
    CubicSpline(x::Vector{Float64}, y::Vector{Float64}) = CubicSpline(CxxWrap.StdVector(x), CxxWrap.StdVector(y))

    function distance_to_border(geom, pos::Cartesian3D, dir::Cartesian3D)
        return (distance_to_border_first(geom, pos, dir),
                distance_to_border_second(geom, pos, dir))
    end
    export distance_to_border

    """
        calculate_secondaries(calc_handle, loss, component, rnd)

    Calculate secondary particles. Returns a Vector{ParticleState}.
    `calc_handle` is the integer handle from `create_secondaries_calculator`.
    """
    function calculate_secondaries(calc_handle::Integer, loss::StochasticLoss,
                                    comp::Component, rnd::Vector{Float64})
        n = sec_calc_calculate(Int(calc_handle), loss, comp, rnd)
        return [sec_calc_result_at(i) for i in 0:(n-1)]
    end
    export calculate_secondaries
end

end # module
