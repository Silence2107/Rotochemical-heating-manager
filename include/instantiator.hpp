#ifndef INSTANTIATOR_H
#define INSTANTIATOR_H

#include "../include/auxiliaries.h"
#include "../include/constants.h"
#include "../include/cooling.h"
#include "../3rd-party/json/single_include/nlohmann/json.hpp"

#include <fstream>
#include <functional>
#include <vector>
#include <map>
#include <cmath>
#include <string>
#include <sstream>

/// @brief global data powering RHM
namespace instantiator
{
    // (1) EoS setup

    // conversion factors from datafile units to natural units
    double energy_density_conversion,
        pressure_conversion,
        nbar_conversion;

    // filereader
    // SHADOWED -- GETS INSTANTIANED BY THE EXTERNAL INPUTFILE

    // energy density function of baryonic density (natural units)
    std::function<double(double)> energy_density_of_nbar;

    // pressure function of baryonic density (natural units)
    std::function<double(double)> pressure_of_nbar;

    /// @brief baryonic density limits in natural units. _low and _upp represent limits of EoS itself
    /// while _core_limit, _crust_limit represent phase transition boundaries
    double nbar_low,
        nbar_upp,
        nbar_core_limit,
        nbar_crust_limit;
    /// @brief energy density limits in natural units. _low and _upp represent limits of EoS itself <para></para>
    /// while _core_limit represents phase transition boundary
    double edensity_low,
        edensity_core_limit,
        edensity_upp;
    /// @brief pressure limits in natural units. _low and _upp represent limits of EoS itself
    double pressure_low,
        pressure_upp;

    // baryonic density fraction functions of baryonic density (natural units)
    std::map<auxiliaries::phys::Species, std::function<double(double)>> number_densities_of_nbar;

    // fermi momentum functions of baryonic density (natural units)
    std::map<auxiliaries::phys::Species, std::function<double(double)>> k_fermi_of_nbar;

    // effective mass functions of baryonic density (natural units)
    std::map<auxiliaries::phys::Species, std::function<double(double)>> m_stars_of_nbar;

    std::function<double(double)> ion_volume_fr;

    // (2) TOV solver setup

    // Cached EoS interpolator. Only use it if you want to erase cache
    auto eos_interpolator_cached = auxiliaries::math::CachedFunc<std::function<double(double)>,
                                                                 double, const std::vector<double> &, const std::vector<double> &,
                                                                 auxiliaries::math::InterpolationMode, double, bool, bool>(auxiliaries::math::interpolate_cached);
    // Interpolator used for EoS P(rho)
    std::function<double(const std::vector<double> &, const std::vector<double> &, double)> eos_interpolator;

    // nbar(r) cached interpolator. Only use it if you want to erase cache
    auto nbar_interpolator_cached = auxiliaries::math::CachedFunc<std::function<double(double)>,
                                                                  double, const std::vector<double> &, const std::vector<double> &,
                                                                  auxiliaries::math::InterpolationMode, double, bool, bool>(auxiliaries::math::interpolate_cached);
    // Interpolator used for nbar(r)
    std::function<double(const std::vector<double> &, const std::vector<double> &, double)> nbar_interpolator;

    // EoS linspace discretization
    size_t discr_size_EoS;

    // TOV adaption limit
    size_t tov_adapt_limit;

    // TOV solver radius step size in GeV
    double radius_step;

    // TOV solver surface pressure in GeV^4
    double surface_pressure;

    // TOV solver center density in GeV^4
    double center_density;

    // (3) Cooling solver

    // Cooling solver setup

    // desirable relative accuracy of the cooling solvers
    double cooling_newton_step_eps;

    // maximum number of iterations of the cooling solvers
    size_t cooling_newton_max_iter;

    // desirable relative accuracy of the cooling solvers per time step
    double cooling_max_diff_per_t_step;

    // Cooling settings
    double crust_eta;

    // Critical phenomena settings

    bool superfluid_p_1s0,
        superfluid_n_3p2,
        superfluid_n_1s0;

    std::function<double(double)> superfluid_p_temp;
    std::function<double(double)> superfluid_n_temp;
    std::function<double(double)> superconduct_q_gap;

    // Evolution settings
    double t_init,
        t_end,
        base_t_step;
    // estimate for the number of time points (is also used for time step expansion, if enabled)
    double cooling_n_points_estimate;
    // initial temperature profile
    std::function<double(double, double, const std::function<double(double)> &, const std::function<double(double)> &)> initial_t_profile_inf;
    // cooling grid step
    double cooling_radius_step;
    // condition on which to switch to equilibrium cooling
    std::function<bool(double, const std::vector<double> &)> switch_to_equilibrium;

    // time step expansion rate (set to 1.0 for constant time step)
    double exp_rate_estim;

    // (4) Rotochemical heating setup

    // Bij matrix independent entries density dependence on nbar
    // ee, emu, eu, es, mumu, ss
    std::function<double(double)> dne_to_dmue,
        dne_to_dmum,
        dnm_to_dmum,
        dnu_to_dmuu,
        dnu_to_dmus,
        dns_to_dmus;

    // rotational 2 omega omega_dot dependency of time
    std::function<double(double)> omega_sqr_dot;

    /// @brief instantiate the system from json input
    /// @param json_input json inputfile path
    void instantiate_system(const std::string &json_input)
    {
        using namespace instantiator;
        using json = nlohmann::json;

        // bunch of simplifying definitions
        using keys = std::vector<std::string>;
        auto get_interpolation_mode = [](const std::string &mode)
        {
            if (mode == "Linear")
                return auxiliaries::math::InterpolationMode::kLinear;
            else if (mode == "Cubic")
                return auxiliaries::math::InterpolationMode::kCubic;
            else
                RHM_THROW(std::runtime_error, "UI error: Unparsable interpolation mode provided.");
        };

        // read json input
        std::ifstream i(json_input);
        if (!i.is_open())
        {
            RHM_THROW(std::runtime_error, "UI inputfile requested, but the path is invalid.");
        }
        json j = json::parse(i);

        // (1) EoS setup

        // filereader
        auto eos_datafile = j["EoSSetup"]["Datafile"]["Path"];
        if (!(eos_datafile.is_string()))
            RHM_THROW(std::runtime_error, "UI error: Datafile path must be provided as a string.");
        auto eos_datafile_rows = j["EoSSetup"]["Datafile"]["Rows"];
        if (!(eos_datafile_rows.size() == 2))
        {
            if (eos_datafile_rows.is_null())
                eos_datafile_rows = {0, 0};
            else
                RHM_THROW(std::runtime_error, "UI error: Datafile rows must be a pair-array.");
        }
        else if (!(eos_datafile_rows[0].is_number_integer() && eos_datafile_rows[1].is_number_integer()))
            RHM_THROW(std::runtime_error, "UI error: Datafile rows must be provided as integers.");

        auto eos_datafile_cols = j["EoSSetup"]["Datafile"]["Columns"];
        if (!(eos_datafile_cols.size() == 2))
        {
            if (eos_datafile_cols.is_null())
                eos_datafile_cols = {0, 0};
            else
                RHM_THROW(std::runtime_error, "UI error: Datafile cols must be a pair-array.");
        }
        else if (!(eos_datafile_cols[0].is_number_integer() && eos_datafile_cols[1].is_number_integer()))
            RHM_THROW(std::runtime_error, "UI error: Datafile cols must be provided as integers.");

        auxiliaries::math::InterpolationMode eos_datafile_interp_mode;
        auto eos_datafile_interp_read = j["EoSSetup"]["Datafile"]["Interpolation"];
        if (eos_datafile_interp_read.is_null())
            eos_datafile_interp_mode = auxiliaries::math::InterpolationMode::kLinear;
        else if (!(eos_datafile_interp_read.is_string()))
            RHM_THROW(std::runtime_error, "UI error: Datafile interpolation mode may only be a string.");
        else
            eos_datafile_interp_mode = get_interpolation_mode(eos_datafile_interp_read);

        auto nbar_index = j["EoSSetup"]["Quantities"]["BarionicDensity"]["Column"];
        auto nbar_conversion_read = j["EoSSetup"]["Quantities"]["BarionicDensity"]["Units"];
        if (!(nbar_index.is_number_integer()))
            RHM_THROW(std::runtime_error, "UI error: Barionic density column number must be provided as an integer.");
        if (nbar_conversion_read.is_number())
            nbar_conversion = nbar_conversion_read.get<double>();
        else if (nbar_conversion_read.is_string())
        {
            if (nbar_conversion_read == "Gev3")
            {
                nbar_conversion = 1.0;
            }
            else if (nbar_conversion_read == "Fm-3")
            {
                nbar_conversion = 1.0 / constants::conversion::fm3_gev3;
            }
            else
            {
                RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for barionic density.");
            }
        }
        else
            RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for barionic density.");

        // read datafile
        static auto table = auxiliaries::io::read_tabulated_file(eos_datafile, eos_datafile_cols, eos_datafile_rows);

        // data_reader takes input vector and outputs vector of outputs from EoS datafile
        static auto data_reader = auxiliaries::math::CachedFunc<std::vector<auxiliaries::math::CachedFunc<std::function<double(double)>,
                                                                                                          double, const std::vector<double> &, const std::vector<double> &,
                                                                                                          auxiliaries::math::InterpolationMode, double, bool, bool>>,
                                                                double, const std::vector<double> &, size_t>(
            [nbar_index, eos_datafile_interp_mode](std::vector<auxiliaries::math::CachedFunc<std::function<double(double)>,
                                                                                             double, const std::vector<double> &, const std::vector<double> &,
                                                                                             auxiliaries::math::InterpolationMode, double, bool, bool>> &cache,
                                                   const std::vector<double> &input, size_t index)
            {
                if (cache.empty())
                {
                    // fill cache with cached interpolation functions for each column
                    for (size_t i = 0; i < table.size(); ++i)
                    {
                        auto interpolator_cached = auxiliaries::math::CachedFunc<std::function<double(double)>,
                                                                                 double, const std::vector<double> &, const std::vector<double> &,
                                                                                 auxiliaries::math::InterpolationMode, double, bool, bool>(auxiliaries::math::interpolate_cached);
                        cache.push_back(interpolator_cached);
                    }
                }
                // unpack input and convert to datafile units
                double nbar = input[0] / nbar_conversion;
                // return cached interpolation functions, with extrapolation enabled for now
                try
                {
                    return cache[index](table.at(nbar_index), table.at(index), eos_datafile_interp_mode, nbar, true, true);
                }
                catch (std::exception &e)
                {
                    RHM_THROW(std::runtime_error, "Bad EoS request. " + e.what());
                }
            });

        // energy density function of baryonic density (natural units)
        auto energy_density_index = j["EoSSetup"]["Quantities"]["EnergyDensity"]["Column"];
        auto energy_density_conversion_read = j["EoSSetup"]["Quantities"]["EnergyDensity"]["Units"];
        if (!(energy_density_index.is_number_integer()))
            RHM_THROW(std::runtime_error, "UI error: Energy density column number must be provided as an integer.");
        if (energy_density_conversion_read.is_number())
            energy_density_conversion = energy_density_conversion_read.get<double>();
        else if (energy_density_conversion_read.is_string())
        {
            if (energy_density_conversion_read == "Gev4")
            {
                energy_density_conversion = 1.0;
            }
            else if (energy_density_conversion_read == "MevFm-3")
            {
                energy_density_conversion = constants::conversion::mev_over_fm3_gev4;
            }
            else if (energy_density_conversion_read == "GCm-3")
            {
                energy_density_conversion = constants::conversion::g_over_cm3_gev4;
            }
            else
            {
                RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for energy density.");
            }
        }
        else
            RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for energy density.");
        energy_density_of_nbar = [energy_density_index](double nbar)
        { return data_reader({nbar}, energy_density_index) * energy_density_conversion; };

        // pressure function of baryonic density (natural units)
        auto pressure_index = j["EoSSetup"]["Quantities"]["Pressure"]["Column"];
        auto pressure_conversion_read = j["EoSSetup"]["Quantities"]["Pressure"]["Units"];
        if (!(pressure_index.is_number_integer()))
            RHM_THROW(std::runtime_error, "UI error: Pressure column number must be provided as an integer.");
        if (pressure_conversion_read.is_number())
            pressure_conversion = pressure_conversion_read.get<double>();
        else if (pressure_conversion_read.is_string())
        {
            if (pressure_conversion_read == "Gev4")
            {
                pressure_conversion = 1.0;
            }
            else if (pressure_conversion_read == "MevFm-3")
            {
                pressure_conversion = constants::conversion::mev_over_fm3_gev4;
            }
            else if (pressure_conversion_read == "DyneCm-2")
            {
                pressure_conversion = constants::conversion::dyne_over_cm2_gev4;
            }
            else
            {
                RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for pressure.");
            }
        }
        else
            RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for pressure.");
        pressure_of_nbar = [pressure_index](double nbar)
        { return data_reader({nbar}, pressure_index) * pressure_conversion; };

        // baryonic density limits in natural units. _low and _upp represent limits of EoS itself
        // while _core_limit, _crust_limit represent phase transition boundaries
        auto nbar_low_read = j["EoSSetup"]["Quantities"]["BarionicDensity"]["Low"],
             nbar_upp_read = j["EoSSetup"]["Quantities"]["BarionicDensity"]["Upp"],
             nbar_core_limit_read = j["EoSSetup"]["Quantities"]["BarionicDensity"]["CoreLimit"],
             nbar_crust_limit_read = j["EoSSetup"]["Quantities"]["BarionicDensity"]["CrustLimit"];
        if (!(nbar_low_read.is_number() && nbar_upp_read.is_number() && nbar_core_limit_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: Barionic density limits must be provided as numbers.");

        if (!nbar_crust_limit_read.is_number() && !nbar_crust_limit_read.is_null())
            RHM_THROW(std::runtime_error, "UI error: Barionic density crust limit may only be provided as a number.");
        if (nbar_crust_limit_read.is_null())
            nbar_crust_limit_read = nbar_core_limit_read;

        nbar_low = nbar_low_read.get<double>() * nbar_conversion;
        nbar_upp = nbar_upp_read.get<double>() * nbar_conversion;
        nbar_core_limit = nbar_core_limit_read.get<double>() * nbar_conversion;
        nbar_crust_limit = nbar_crust_limit_read.get<double>() * nbar_conversion;

        // energy density limits in natural units. _low and _upp represent limits of EoS itself <para></para>
        // while _core_limit represents phase transition boundary
        auto edensity_low_read = j["EoSSetup"]["Quantities"]["EnergyDensity"]["Low"],
             edensity_core_limit_read = j["EoSSetup"]["Quantities"]["EnergyDensity"]["CoreLimit"],
             edensity_upp_read = j["EoSSetup"]["Quantities"]["EnergyDensity"]["Upp"];

        if (edensity_low_read.is_null() || edensity_low_read == "Deduce")
            edensity_low = energy_density_of_nbar(nbar_low);
        else if (edensity_low_read.is_number())
            edensity_low = edensity_low_read.get<double>() * energy_density_conversion;
        else
            RHM_THROW(std::runtime_error, "UI error: Energy density low limit may only be provided as a number or \"Deduce\".");

        if (edensity_core_limit_read.is_null() || edensity_core_limit_read == "Deduce")
            edensity_core_limit = energy_density_of_nbar(nbar_core_limit);
        else if (edensity_core_limit_read.is_number())
            edensity_core_limit = edensity_core_limit_read.get<double>() * energy_density_conversion;
        else
            RHM_THROW(std::runtime_error, "UI error: Energy density core limit may only be provided as a number or \"Deduce\".");

        if (edensity_upp_read.is_null() || edensity_upp_read == "Deduce")
            edensity_upp = energy_density_of_nbar(nbar_upp);
        else if (edensity_upp_read.is_number())
            edensity_upp = edensity_upp_read.get<double>() * energy_density_conversion;
        else
            RHM_THROW(std::runtime_error, "UI error: Energy density upp limit may only be provided as a number or \"Deduce\".");

        // pressure limits in natural units. _low and _upp represent limits of EoS itself
        auto pressure_low_read = j["EoSSetup"]["Quantities"]["Pressure"]["Low"],
             pressure_upp_read = j["EoSSetup"]["Quantities"]["Pressure"]["Upp"];

        if (pressure_low_read.is_null() || pressure_low_read == "Deduce")
            pressure_low = pressure_of_nbar(nbar_low);
        else if (pressure_low_read.is_number())
            pressure_low = pressure_low_read.get<double>() * pressure_conversion;
        else
            RHM_THROW(std::runtime_error, "UI error: Pressure low limit may only be provided as a number or \"Deduce\".");

        if (pressure_upp_read.is_null() || pressure_upp_read == "Deduce")
            pressure_upp = pressure_of_nbar(nbar_upp);
        else if (pressure_upp_read.is_number())
            pressure_upp = pressure_upp_read.get<double>() * pressure_conversion;
        else
            RHM_THROW(std::runtime_error, "UI error: Pressure upp limit may only be provided as a number or \"Deduce\".");

        // (2) TOV solver setup

        if (j["TOVSolver"].is_null())
            RHM_THROW(std::runtime_error, "UI error: TOV solver setup is essential for any simulation and must be provided.");

        auxiliaries::math::InterpolationMode eos_interp_mode;
        auto eos_interp_mode_read = j["TOVSolver"]["EoSInterpolation"];
        if (eos_interp_mode_read.is_null())
            eos_interp_mode = auxiliaries::math::InterpolationMode::kLinear;
        else if (!(eos_interp_mode_read.is_string()))
            RHM_THROW(std::runtime_error, "UI error: EoS interpolation mode must be a string.");
        else
            eos_interp_mode = get_interpolation_mode(eos_interp_mode_read);

        // Interpolator used for EoS P(rho)
        eos_interpolator = [eos_interp_mode](const std::vector<double> &input, const std::vector<double> &output, double val)
        {
            return eos_interpolator_cached(input, output, eos_interp_mode, val, false, true);
        };

        auxiliaries::math::InterpolationMode nbar_interp_mode;
        auto nbar_interp_mode_read = j["TOVSolver"]["BarionicDensityInterpolation"];
        if (nbar_interp_mode_read.is_null())
            nbar_interp_mode = auxiliaries::math::InterpolationMode::kLinear;
        else if (!(nbar_interp_mode_read.is_string()))
            RHM_THROW(std::runtime_error, "UI error: nbar interpolation mode must be a string.");
        else
            nbar_interp_mode = get_interpolation_mode(nbar_interp_mode_read);

        // Interpolator used for nbar(r)
        nbar_interpolator = [nbar_interp_mode](const std::vector<double> &input, const std::vector<double> &output, double val)
        {
            return nbar_interpolator_cached(input, output, nbar_interp_mode, val, false, true);
        };

        // EoS linspace discretization
        auto discr_size_EoS_read = j["TOVSolver"]["EoSDiscretization"];
        if (discr_size_EoS_read.is_null())
            discr_size_EoS = 1000;
        else if (!(discr_size_EoS_read.is_number_integer()))
            RHM_THROW(std::runtime_error, "UI error: EoS discretization must be provided as an integer.");
        else
            discr_size_EoS = discr_size_EoS_read.get<size_t>();

        // TOV adaption limit
        auto tov_adapt_limit_read = j["TOVSolver"]["AdaptionLimit"];
        if (tov_adapt_limit_read.is_null())
            tov_adapt_limit = 20;
        else if (!(tov_adapt_limit_read.is_number_integer()))
            RHM_THROW(std::runtime_error, "UI error: TOV adaption limit must be provided as an integer.");
        else
            tov_adapt_limit = tov_adapt_limit_read.get<size_t>();

        // TOV solver radius step size in GeV
        auto radius_step_read = j["TOVSolver"]["RadiusStep"];
        auto tov_length_conversion_read = j["TOVSolver"]["LengthUnits"];
        double tov_length_conversion;
        if (tov_length_conversion_read.is_number())
            tov_length_conversion = tov_length_conversion_read.get<double>();
        else if (tov_length_conversion_read.is_string())
        {
            if (tov_length_conversion_read == "Gev-1")
            {
                tov_length_conversion = 1.0;
            }
            else if (tov_length_conversion_read == "Km")
            {
                tov_length_conversion = constants::conversion::km_gev;
            }
            else if (tov_length_conversion_read == "M")
            {
                tov_length_conversion = 1E-3 * constants::conversion::km_gev;
            }
            else if (tov_length_conversion_read == "Cm")
            {
                tov_length_conversion = 1E-5 * constants::conversion::km_gev;
            }
            else
            {
                RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for TOV length.");
            }
        }
        else
            RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for TOV length.");
        if (!(radius_step_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: TOV solver radius step must be provided as a number.");
        else
            radius_step = radius_step_read.get<double>() * tov_length_conversion;

        // TOV solver surface pressure & center density
        auto tov_surface_pressure_provided_as_read = j["TOVSolver"]["SurfacePressure"]["ProvidedAs"];
        auto tov_surface_pressure_read = j["TOVSolver"]["SurfacePressure"]["Value"];
        if (!(tov_surface_pressure_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: TOV surface pressure must be provided as a number.");

        if (tov_surface_pressure_provided_as_read == "LinspacedMinToMax")
        {
            surface_pressure = tov_surface_pressure_read.get<double>() * (pressure_upp - pressure_low) + pressure_low;
        }
        else if (tov_surface_pressure_provided_as_read == "Same")
        {
            surface_pressure = tov_surface_pressure_read.get<double>() * pressure_conversion;
        }
        else
        {
            RHM_THROW(std::runtime_error, "UI error: TOV surface pressure may only be provided in \"LinspacedMinToMax\" or \"Same\" modes.");
        }

        auto center_density_read = j["TOVSolver"]["CenterDensity"]["Value"];
        auto tov_density_provided_as_read = j["TOVSolver"]["CenterDensity"]["ProvidedAs"];
        if (!(center_density_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: TOV solver center density must be provided as a number.");

        if (tov_density_provided_as_read == "LinspacedMinToMax")
        {
            center_density = center_density_read.get<double>() * (edensity_upp - edensity_low) + edensity_low;
        }
        else if (tov_density_provided_as_read == "Same")
        {
            center_density = center_density_read.get<double>() * energy_density_conversion;
        }
        else
        {
            RHM_THROW(std::runtime_error, "UI error: TOV center density may only be provided in \"LinspacedMinToMax\" or \"Same\" modes.");
        }

        // (->1) EoS Setup

        if (j["CoolingSolver"].is_null())
            return; // user presents no interest in cooling

        // provided particles
        auto particles_read = j["EoSSetup"]["Particles"];
        if (!particles_read.is_array())
            RHM_THROW(std::runtime_error, "UI error: Particle types must be provided as an array.");

        std::vector<auxiliaries::phys::Species> particles;
        for (auto particle_name : particles_read)
        {
            using namespace constants::species;
            std::vector<auxiliaries::phys::Species> known_particles =
                {neutron, proton, electron, muon, tau, uquark, dquark, squark};
            if (particle_name.is_string())
            {
                // identify particle with the list of known ones
                bool identified = false;
                for (const auto &known_particle : known_particles)
                {
                    if (particle_name == known_particle.name())
                    {
                        particles.push_back(known_particle);
                        identified = true;
                        break;
                    }
                }
                if (!identified)
                    RHM_THROW(std::runtime_error, "UI error: " + particle_name.get<std::string>().c_str() + " is not a recognized particle.");
            }
            else
                RHM_THROW(std::runtime_error, "UI error: Particle type must be provided as a string.");
        }

        // number density functions of baryonic density (natural units) &&
        // fermi momentum functions of baryonic density (natural units)

        for (const auto &particle : particles)
        {
            auto particle_name = particle.name();
            auto particle_density_read = j["EoSSetup"]["Quantities"]["NumberDensities"][particle_name];

            if (particle_density_read.is_null())
                RHM_THROW(std::runtime_error, "UI error: Particle density info is not provided for " + particle.name() + ".");

            auto particle_density_column_read = particle_density_read["Column"];
            if (!(particle_density_column_read.is_number_integer()))
                RHM_THROW(std::runtime_error, "UI error: " + particle.name() + " density column number must be provided as an integer.");

            auto particle_density_conversion_read = particle_density_read["Units"];
            double particle_density_conversion;
            if (particle_density_conversion_read.is_null())
                particle_density_conversion = nbar_conversion;
            else if (particle_density_conversion_read.is_number())
            {
                particle_density_conversion = particle_density_conversion_read.get<double>();
            }
            else if (particle_density_conversion_read.is_string())
            {
                // for "Density"
                if (particle_density_conversion_read == "Gev3")
                {
                    particle_density_conversion = 1.0;
                }
                else if (particle_density_conversion_read == "Fm-3")
                {
                    particle_density_conversion = 1.0 / constants::conversion::fm3_gev3;
                }
                // for "DensityFraction"
                else if (particle_density_conversion_read == "DimLess")
                {
                    particle_density_conversion = 1.0;
                }
                // for "KFermi"
                else if (particle_density_conversion_read == "Gev")
                {
                    particle_density_conversion = 1.0;
                }
                else if (particle_density_conversion_read == "Fm-1")
                {
                    particle_density_conversion = pow(1.0 / constants::conversion::fm3_gev3, 1.0 / 3.0);
                }
                else
                {
                    RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for " + particle_name + " density.");
                }
            }
            else
                RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for " + particle_name + " density.");

            auto particle_density_provided_as_read = particle_density_read["ProvidedAs"];
            if (particle_density_provided_as_read.is_null())
                particle_density_provided_as_read = "Density";

            if (particle_density_provided_as_read == "Density")
            {
                number_densities_of_nbar.insert(
                    {particle, [particle_density_column_read, particle_density_conversion](double nbar)
                     {
                         return data_reader({nbar}, particle_density_column_read) * particle_density_conversion;
                     }});
            }
            else if (particle_density_provided_as_read == "DensityFraction")
            {
                number_densities_of_nbar.insert(
                    {particle, [particle_density_column_read, particle_density_conversion](double nbar)
                     {
                         return data_reader({nbar}, particle_density_column_read) * particle_density_conversion * nbar;
                     }});
            }
            else if (particle_density_provided_as_read == "KFermi")
            {
                number_densities_of_nbar.insert(
                    {particle, [particle_density_column_read, particle_density_conversion](double nbar)
                     {
                         using constants::scientific::Pi;
                         return pow(data_reader({nbar}, particle_density_column_read) * particle_density_conversion, 3.0) / (3.0 * Pi * Pi);
                     }});
            }
            else
                RHM_THROW(std::runtime_error, "UI error: Particle density may only be provided in \"Density\", \"DensityFraction\" or \"KFermi\" modes.");
            // assemble fermi momentum functions for fermions
            if (particle.classify() != auxiliaries::phys::Species::ParticleClassification::kMeson)
                k_fermi_of_nbar.insert(
                    {particle, [particle](double nbar)
                     {
                         using constants::scientific::Pi;
                         return pow(3.0 * Pi * Pi * number_densities_of_nbar[particle](nbar), 1.0 / 3.0);
                     }});
        }

        // effective mass functions of baryonic density (natural units)

        for (const auto &particle : particles)
        {
            auto particle_name = particle.name();
            auto particle_mst_read = j["EoSSetup"]["Quantities"]["EffectiveMasses"][particle_name];

            if (particle.classify() == auxiliaries::phys::Species::ParticleClassification::kMeson)
            {
                if (!particle_mst_read.is_null())
                    RHM_THROW(std::runtime_error, "UI error: Effective mass is not expected for " + particle.name() + ".");
                continue; // mesons have no effective mass
            }

            if (particle_mst_read.is_null())
                RHM_THROW(std::runtime_error, "UI error: Effective mass info is not provided for " + particle.name() + ".");
            auto particle_mst_column_read = particle_mst_read["Column"];
            auto particle_mst_provided_as_read = particle_mst_read["ProvidedAs"];
            if (!(particle_mst_column_read.is_number_integer()) && particle_mst_provided_as_read != "FermiEnergy")
                RHM_THROW(std::runtime_error, "UI error: " + particle.name() + " effective mass column number must be provided as an integer.");

            auto particle_mst_conversion_read = particle_mst_read["Units"];
            double particle_mst_conversion;
            if (particle_mst_conversion_read.is_null())
            {
                if (particle_mst_provided_as_read != "FermiEnergy")
                    RHM_THROW(std::runtime_error, "UI error: Effective mass units must be provided for " + particle.name() + ".");
            }
            else if (particle_mst_conversion_read.is_number())
            {
                particle_mst_conversion = particle_mst_conversion_read.get<double>();
            }
            else if (particle_mst_conversion_read.is_string())
            {
                if (particle_mst_conversion_read == "Gev")
                {
                    particle_mst_conversion = 1.0;
                }
                else if (particle_mst_conversion_read == "Mev")
                {
                    particle_mst_conversion = 1.0 / constants::conversion::gev_over_mev;
                }
                else if (particle_mst_conversion_read == "NucleonMass")
                {
                    particle_mst_conversion = constants::species::neutron.mass();
                }
                else
                {
                    RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for " + particle.name() + " effective mass.");
                }
            }
            else
                RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for " + particle.name() + " effective mass.");
            if (particle_mst_provided_as_read == "FermiEnergy")
            {
                m_stars_of_nbar.insert(
                    {particle, [particle](double nbar)
                     {
                         return sqrt(pow(k_fermi_of_nbar[particle](nbar), 2.0) + pow(particle.mass(), 2.0));
                     }});
            }
            else if (particle_mst_provided_as_read == "EffectiveMass")
            {
                m_stars_of_nbar.insert(
                    {particle, [particle_mst_column_read, particle_mst_conversion](double nbar)
                     {
                         return data_reader({nbar}, particle_mst_column_read) * particle_mst_conversion;
                     }});
            }
            else
                RHM_THROW(std::runtime_error, "UI error: Particle effective mass may only be provided in \"FermiEnergy\" or \"EffectiveMass\" modes.");
        }

        // ion volume fraction function of baryonic density (natural units)

        auto ion_volume_fr_read = j["EoSSetup"]["Quantities"]["IonVolumeFraction"];

        auto ion_volume_provided_as_read = ion_volume_fr_read["ProvidedAs"];
        if (ion_volume_provided_as_read.is_null() || ion_volume_provided_as_read == "Absent")
            ion_volume_fr = [](double nbar)
            { return 0.0; };
        else if (ion_volume_provided_as_read == "ExcludedVolume")
            ion_volume_fr = [](double nbar)
            {
                if (nbar >= nbar_crust_limit)
                    return 0.0;
                using namespace constants::conversion;
                using namespace constants::scientific;
                auto eta_ion = 4.0 / 3 * Pi * pow(1.1, 3.0) * fm3_gev3 * energy_density_of_nbar(nbar) / constants::species::neutron.mass();
                return std::min(1.0, eta_ion);
            };
        else if (ion_volume_provided_as_read == "IonVolumeFraction")
        {
            auto ion_volume_column_read = ion_volume_fr_read["Column"];
            if (!(ion_volume_column_read.is_number_integer()))
                RHM_THROW(std::runtime_error, "UI error: Ion volume fraction column number must be provided as an integer.");
            ion_volume_fr = [ion_volume_column_read](double nbar)
            {
                return data_reader({nbar}, ion_volume_column_read);
            };
        }
        else
            RHM_THROW(std::runtime_error, "UI error: Ion volume fraction may only be provided in \"Absent\", \"ExcludedVolume\" or \"IonVolumeFraction\" modes.");

        // (3) Cooling solver

        // Cooling solver setup

        // desirable relative accuracy of the cooling solvers
        auto cooling_newton_step_eps_read = j["CoolingSolver"]["NewtonTolerance"];
        if (cooling_newton_step_eps_read.is_null())
            cooling_newton_step_eps = 1E-5;
        else if (!(cooling_newton_step_eps_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: Cooling solver Newton relative tolerance must be provided as a number.");
        else
            cooling_newton_step_eps = cooling_newton_step_eps_read.get<double>();

        // maximum number of iterations of the cooling solvers
        auto cooling_newton_max_iter_read = j["CoolingSolver"]["NewtonMaxIter"];
        if (cooling_newton_max_iter_read.is_null())
            cooling_newton_max_iter = 50;
        else if (!(cooling_newton_max_iter_read.is_number_integer()))
            RHM_THROW(std::runtime_error, "UI error: Cooling solver Newton max iter must be provided as an integer.");
        else
            cooling_newton_max_iter = cooling_newton_max_iter_read.get<size_t>();

        // desirable relative accuracy of the cooling solvers per time step
        auto cooling_newton_step_eps_per_step_read = j["CoolingSolver"]["StepTolerance"];
        if (cooling_newton_step_eps_per_step_read.is_null())
            cooling_newton_step_eps = 0.05;
        else if (!(cooling_newton_step_eps_per_step_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: Cooling solver relative tolerance per step must be provided as a number.");
        else
            cooling_newton_step_eps = cooling_newton_step_eps_per_step_read.get<double>();

        // initial temperature profile
        auto initial_t_profile_read = j["CoolingSolver"]["TemperatureProfile"];
        auto cooling_temp_conversion_read = j["CoolingSolver"]["TemperatureUnits"];
        double cooling_temp_conversion;
        if (cooling_temp_conversion_read.is_number())
            cooling_temp_conversion = cooling_temp_conversion_read.get<double>();
        else if (cooling_temp_conversion_read.is_string())
        {
            if (cooling_temp_conversion_read == "Gev")
            {
                cooling_temp_conversion = 1.0;
            }
            else if (cooling_temp_conversion_read == "Mev")
            {
                cooling_temp_conversion = 1.0 / constants::conversion::gev_over_mev;
            }
            else if (cooling_temp_conversion_read == "K")
            {
                cooling_temp_conversion = 1.0 / constants::conversion::gev_over_k;
            }
            else
            {
                RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for temperature.");
            }
        }
        else
            RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for temperature.");

        auto initial_t_profile_provided_as_read = initial_t_profile_read["ProvidedAs"];
        std::function<double(double, double, const std::function<double(double)> &)> redshift_factor;
        if (initial_t_profile_provided_as_read == "Redshifted")
        {
            redshift_factor = [](double r, double r_ns, const std::function<double(double)> &exp_phi)
            {
                return 1;
            };
        }
        else if (initial_t_profile_provided_as_read == "SurfaceRedshifted")
        {
            redshift_factor = [](double r, double r_ns, const std::function<double(double)> &exp_phi)
            {
                return exp_phi(r_ns);
            };
        }
        else if (initial_t_profile_provided_as_read == "Local")
        {
            redshift_factor = [](double r, double r_ns, const std::function<double(double)> &exp_phi)
            {
                return exp_phi(r);
            };
        }
        else
            RHM_THROW(std::runtime_error, "UI error: Initial temperature profile must be provided as \"Redshifted\", \"SurfaceRedshifted\" or \"Local\" .");

        auto initial_t_profile_mode_read = initial_t_profile_read["Mode"];
        auto initial_t_profile_params_read = initial_t_profile_read["Parameters"];
        if (!initial_t_profile_params_read.is_array())
            RHM_THROW(std::runtime_error, "UI error: Initial temperature profile parameters must be provided as an array.");

        if (initial_t_profile_mode_read == "Flat")
        {
            // expected parameters: [T]
            if (initial_t_profile_params_read.size() != 1)
                RHM_THROW(std::runtime_error, "UI error: Initial temperature profile parameters must be of size 1 for the \"Flat\" mode.");
            auto temp_init_read = initial_t_profile_params_read[0];
            double temp_init;
            if (!(temp_init_read.is_number()))
                RHM_THROW(std::runtime_error, "UI error: Initial temperature must be provided as a number.");
            else
                temp_init = temp_init_read.get<double>() * cooling_temp_conversion;
            initial_t_profile_inf = [temp_init, redshift_factor](double r, double r_ns, const std::function<double(double)> &exp_phi, const std::function<double(double)> &nbar_of_r)
            {
                return temp_init * redshift_factor(r, r_ns, exp_phi);
            };
        }
        else if (initial_t_profile_mode_read == "CoreJump")
        {
            // expected parameters: [T_core, T_crust]
            if (initial_t_profile_params_read.size() != 2)
                RHM_THROW(std::runtime_error, "UI error: Initial temperature profile parameters must be of size 2 for the \"CoreJump\" mode.");
            auto temp_core_read = initial_t_profile_params_read[0];
            auto temp_crust_read = initial_t_profile_params_read[1];
            double temp_core, temp_crust;
            if (!(temp_core_read.is_number()))
                RHM_THROW(std::runtime_error, "UI error: Initial temperature in the core must be provided as a number.");
            else
                temp_core = temp_core_read.get<double>() * cooling_temp_conversion;
            if (!(temp_crust_read.is_number()))
                RHM_THROW(std::runtime_error, "UI error: Initial temperature in the crust must be provided as a number.");
            else
                temp_crust = temp_crust_read.get<double>() * cooling_temp_conversion;

            initial_t_profile_inf = [temp_core, temp_crust, redshift_factor](double r, double r_ns, const std::function<double(double)> &exp_phi, const std::function<double(double)> &nbar_of_r)
            {
                return (nbar_of_r(r) > nbar_core_limit ? temp_core : temp_crust) * redshift_factor(r, r_ns, exp_phi);
            };
        }
        else
            RHM_THROW(std::runtime_error, "UI error: Initial temperature profile must be provided in \"Flat\" or \"CoreJump\" mode .");
        // Evolution settings
        auto time_conversion_read = j["CoolingSolver"]["TimeUnits"];
        double time_conversion;
        if (time_conversion_read.is_number())
            time_conversion = time_conversion_read.get<double>();
        else if (time_conversion_read.is_string())
        {
            if (time_conversion_read == "Gev-1")
            {
                time_conversion = 1.0;
            }
            else if (time_conversion_read == "Yr")
            {
                time_conversion = 1E-6 * constants::conversion::myr_over_s * constants::conversion::gev_s;
            }
            else if (time_conversion_read == "Myr")
            {
                time_conversion = constants::conversion::myr_over_s * constants::conversion::gev_s;
            }
            else if (time_conversion_read == "S")
            {
                time_conversion = constants::conversion::gev_s;
            }
            else
            {
                RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for time.");
            }
        }
        else
            RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for time.");
        auto time_init_read = j["CoolingSolver"]["TimeInit"];
        if (time_init_read.is_null())
            t_init = 0.0 * time_conversion;
        else if (!(time_init_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: Initial time may only be provided as a number.");
        else
            t_init = time_init_read.get<double>() * time_conversion;

        auto time_end_read = j["CoolingSolver"]["TimeEnd"];
        if (!(time_end_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: Final time must be provided as a number.");
        else
            t_end = time_end_read.get<double>() * time_conversion;

        auto base_time_step_read = j["CoolingSolver"]["TimeBaseStep"];
        if (!(base_time_step_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: Base time step must be provided as a number.");
        else
            base_t_step = base_time_step_read.get<double>() * time_conversion;

        // estimate for the number of time points (is also used for time step expansion, if enabled)
        auto n_points_estimate_read = j["CoolingSolver"]["NumberPointsEstimate"];
        if (!(n_points_estimate_read.is_number_integer()))
            RHM_THROW(std::runtime_error, "UI error: Estimate of number of cooling time stamps must be provided as an integer.");
        else
            cooling_n_points_estimate = n_points_estimate_read.get<size_t>();

        // time step expansion factor
        auto time_step_expansion_factor_read = j["CoolingSolver"]["ExpansionRate"];
        if (time_step_expansion_factor_read.is_null() || time_step_expansion_factor_read == "Deduce")
            exp_rate_estim = pow((t_end - t_init) / base_t_step, 1.0 / cooling_n_points_estimate) *
                             pow((pow((t_end - t_init) / base_t_step, 1.0 / cooling_n_points_estimate) - 1), 1.0 / cooling_n_points_estimate);
        else if (!(time_step_expansion_factor_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: Time step expansion factor may only be provided as a number or \"Deduce\".");
        else
            exp_rate_estim = time_step_expansion_factor_read.get<double>();

        // cooling grid setup
        auto cooling_radius_step_read = j["CoolingSolver"]["RadiusStep"];
        auto cooling_length_conversion_read = j["CoolingSolver"]["LengthUnits"];
        double cooling_length_conversion;
        if (cooling_length_conversion_read.is_number())
            cooling_length_conversion = cooling_length_conversion_read.get<double>();
        else if (cooling_length_conversion_read.is_string())
        {
            if (cooling_length_conversion_read == "Gev-1")
            {
                cooling_length_conversion = 1.0;
            }
            else if (cooling_length_conversion_read == "Km")
            {
                cooling_length_conversion = constants::conversion::km_gev;
            }
            else if (cooling_length_conversion_read == "M")
            {
                cooling_length_conversion = 1E-3 * constants::conversion::km_gev;
            }
            else if (cooling_length_conversion_read == "Cm")
            {
                cooling_length_conversion = 1E-5 * constants::conversion::km_gev;
            }
            else
            {
                RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for cooling length.");
            }
        }
        else
            RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for cooling length.");
        if (!(cooling_radius_step_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: Cooling solver radius step must be provided as a number.");
        else
            cooling_radius_step = cooling_radius_step_read.get<double>() * cooling_length_conversion;

        // condition on which to switch to equilibrium cooling
        auto cooling_enable_equilibrium_mode_read = j["CoolingSolver"]["EnableEquilibrium"]["Mode"];
        if (!(cooling_enable_equilibrium_mode_read.is_string()))
            switch_to_equilibrium = [](double, const std::vector<double> &)
            { return false; };
        else
        {
            if (cooling_enable_equilibrium_mode_read == "Immediately")
                switch_to_equilibrium = [](double, const std::vector<double> &)
                { return true; };
            else if (cooling_enable_equilibrium_mode_read == "Never")
                switch_to_equilibrium = [](double, const std::vector<double> &)
                { return false; };
            else if (cooling_enable_equilibrium_mode_read == "Conditional")
            {
                auto cooling_enable_equilibrium_condition1_read = j["CoolingSolver"]["EnableEquilibrium"]["Conditions"]["UponReachingTime"];
                auto cooling_enable_equilibrium_condition2_read = j["CoolingSolver"]["EnableEquilibrium"]["Conditions"]["UponProfileFlattening"];

                // let's make it less efficient but more readable
                switch_to_equilibrium = [cooling_enable_equilibrium_condition1_read,
                                         cooling_enable_equilibrium_condition2_read,
                                         time_conversion](double t_curr, const std::vector<double> &t_profile)
                {
                    if (!cooling_enable_equilibrium_condition1_read.is_null())
                    {
                        if (!(cooling_enable_equilibrium_condition1_read.is_number()))
                            RHM_THROW(std::runtime_error, "UI error: Time for switching to equilibrium may only be provided as a number.");
                        else if (t_curr < cooling_enable_equilibrium_condition1_read.get<double>() * time_conversion)
                            return false;
                    }
                    if (!cooling_enable_equilibrium_condition2_read.is_null())
                    {
                        if (!(cooling_enable_equilibrium_condition2_read.is_number()))
                            RHM_THROW(std::runtime_error, "UI error: Profile flattening ratio for switching to equilibrium may only be provided as a number.");
                        else if (std::abs(t_profile.end()[-2] - t_profile.front()) / t_profile.end()[-2] > cooling_enable_equilibrium_condition2_read.get<double>())
                            return false;
                    }
                    return true;
                };
            }
            else
                RHM_THROW(std::runtime_error, "UI error: Equilibrium condition may only be provided in \"Immediately\", \"Never\" or \"Conditional\" modes.");
        }

        // Cooling settings
        auto crust_eta_read = j["EoSSetup"]["Misc"]["CrustalEta"];
        if (crust_eta_read.is_null() || !(crust_eta_read.is_number()))
            RHM_THROW(std::runtime_error, "UI error: Light element share (crustal #eta) must be provided as a number.");
        else
            crust_eta = crust_eta_read.get<double>();

        // Critical phenomena settings
        auto superfluid_p_1s0_read = j["EoSSetup"]["Misc"]["ProtonSuperfluidity1S0"];
        auto superfluid_n_3p2_read = j["EoSSetup"]["Misc"]["NeutronSuperfluidity3P2"];
        auto superfluid_n_1s0_read = j["EoSSetup"]["Misc"]["NeutronSuperfluidity1S0"];

        superfluid_p_1s0 = !superfluid_p_1s0_read.is_null();
        superfluid_n_3p2 = !superfluid_n_3p2_read.is_null();
        superfluid_n_1s0 = !superfluid_n_1s0_read.is_null();

        auto select_crit_temp_model = [](const std::string &model_name)
        {
            /*kAO,
            kCCDK,
            kA,
            kB,
            kC,
            kA2,
            kHadronToQGP*/
            using namespace auxiliaries::phys;
            if (model_name == "AO")
                return CriticalTemperatureModel::kAO;
            else if (model_name == "CCDK")
                return CriticalTemperatureModel::kCCDK;
            else if (model_name == "A")
                return CriticalTemperatureModel::kA;
            else if (model_name == "B")
                return CriticalTemperatureModel::kB;
            else if (model_name == "C")
                return CriticalTemperatureModel::kC;
            else if (model_name == "A2")
                return CriticalTemperatureModel::kA2;
            else if (model_name == "HadronToQGP")
                return CriticalTemperatureModel::kHadronToQGP;
            else
                RHM_THROW(std::runtime_error, "UI error: Critical temperature model must be provided as a string among \"AO\", \"CCDK\", \"A\", \"B\", \"C\", \"A2\" or \"HadronToQGP\".");
        };

        if (!superfluid_p_1s0_read.is_string() && superfluid_p_1s0)
            RHM_THROW(std::runtime_error, "UI error: Proton superfluidity may only be provided as a string (namely model name).");
        if (!superfluid_n_3p2_read.is_string() && superfluid_n_3p2)
            RHM_THROW(std::runtime_error, "UI error: Neutron superfluidity may only be provided as a string (namely model name).");
        if (!superfluid_n_1s0_read.is_string() && superfluid_n_1s0)
            RHM_THROW(std::runtime_error, "UI error: Neutron superfluidity may only be provided as a string (namely model name).");

        superfluid_p_temp = [select_crit_temp_model, superfluid_p_1s0_read](double k_fermi)
        {
            if (superfluid_p_1s0)
            {
                using namespace auxiliaries::phys;
                return critical_temperature(k_fermi, select_crit_temp_model(superfluid_p_1s0_read.get<std::string>()));
            }
            return 0.0;
        };
        superfluid_n_temp = [select_crit_temp_model, superfluid_n_3p2_read, superfluid_n_1s0_read](double k_fermi)
        {
            using namespace auxiliaries::phys;
            using constants::species::neutron;
            if (superfluid_n_3p2 && superfluid_n_1s0)
            {
                if (k_fermi <= k_fermi_of_nbar[neutron](nbar_core_limit))
                    return critical_temperature(k_fermi, select_crit_temp_model(superfluid_n_1s0_read.get<std::string>()));
                else
                    return critical_temperature(k_fermi, select_crit_temp_model(superfluid_n_3p2_read.get<std::string>()));
            }
            else if (superfluid_n_3p2)
                return critical_temperature(k_fermi, select_crit_temp_model(superfluid_n_3p2_read.get<std::string>()));
            else if (superfluid_n_1s0)
                return critical_temperature(k_fermi, select_crit_temp_model(superfluid_n_1s0_read.get<std::string>()));
            else
                return 0.0;
        };

        // superconducting gap
        auto superconduct_gap_conversion_read = j["EoSSetup"]["Quantities"]["QuarkSuperconductingGap"]["Units"];
        auto superconduct_gap_column_read = j["EoSSetup"]["Quantities"]["QuarkSuperconductingGap"]["Column"];
        if (!superconduct_gap_column_read.is_null() && !superconduct_gap_column_read.is_number_integer())
            RHM_THROW(std::runtime_error, "UI error: Quark superconducting gap column number may only be provided as an integer.");

        double superconduct_gap_conversion;
        if (superconduct_gap_conversion_read.is_number())
        {
            superconduct_gap_conversion = superconduct_gap_conversion_read.get<double>();
        }
        else if (superconduct_gap_conversion_read.is_string())
        {
            if (superconduct_gap_conversion_read == "Gev")
            {
                superconduct_gap_conversion = 1.0;
            }
            else if (superconduct_gap_conversion_read == "Mev")
            {
                superconduct_gap_conversion = 1.0 / constants::conversion::gev_over_mev;
            }
            else if (superconduct_gap_conversion_read == "Fm-1")
            {
                superconduct_gap_conversion = pow(1.0 / constants::conversion::fm3_gev3, 1.0 / 3.0);
            }
            else
            {
                RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for quark superconducting gap.");
            }
        }
        else if (!superconduct_gap_column_read.is_null())
            RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for quark superconducting gap.");

        if (!superconduct_gap_column_read.is_null())
        {
            superconduct_q_gap = [superconduct_gap_column_read, superconduct_gap_conversion](double nbar)
            {
                return data_reader({nbar}, superconduct_gap_column_read) * superconduct_gap_conversion;
            };
        }
        else
        {
            superconduct_q_gap = [](double nbar)
            {
                return 0.0;
            };
        }

        // (4) Rotochemical heating setup

        if (j["RHSolver"].is_null())
            return; // user does not want rotochemical heating

        auto bij_density_read = j["EoSSetup"]["Quantities"]["DensityChemPotentialDerivatives"];

        // we expect npemuds (H{npem} + Q{udse}) matter at most.
        using namespace constants::species;

        // we first see which entries are provided by the user

        auto bee_density_read = bij_density_read["Electron"]["Electron"],
             bem_density_read = bij_density_read["Electron"]["Muon"],
             bmm_density_read = bij_density_read["Muon"]["Muon"],
             buu_density_read = bij_density_read["Uquark"]["Uquark"],
             bus_density_read = bij_density_read["Uquark"]["Squark"],
             bss_density_read = bij_density_read["Squark"]["Squark"];
        bem_density_read = (bem_density_read.is_null() ? bij_density_read["Muon"]["Electron"] : bem_density_read);
        bus_density_read = (bus_density_read.is_null() ? bij_density_read["Squark"]["Uquark"] : bus_density_read);

        // if all npemuds matter plays role in rotochemical heating,
        // the number of rows is reduced to 4 (via conservation laws),
        // namely, (dNe, dNm, dNu, dNs), with 6 independent entries within the matrix.

        // all checks reflecting whether the user provided any good info will be done in the main program

        auto bij_vars = std::vector<std::reference_wrapper<decltype(dne_to_dmue)>>({dne_to_dmue, dne_to_dmum, dnm_to_dmum, dnu_to_dmuu, dnu_to_dmus, dns_to_dmus});
        auto bij_reads = std::vector<decltype(bee_density_read)>({bee_density_read, bem_density_read, bmm_density_read, buu_density_read, bus_density_read, bss_density_read});

        for (size_t entry_num = 0; entry_num < bij_vars.size(); ++entry_num)
        {
            auto read = bij_reads[entry_num];
            auto &entry = bij_vars[entry_num];
            if (!read.is_null())
            {
                auto column_read = read["Column"];
                auto conversion_read = read["Units"];
                double conversion;
                if (conversion_read.is_null())
                    RHM_THROW(std::runtime_error, "UI error: Units must be provided for each \"DensityChemPotentialDerivatives\" entry.");
                else if (conversion_read.is_number())
                    conversion = conversion_read.get<double>();
                else if (conversion_read.is_string())
                {
                    if (conversion_read == "Gev2")
                    {
                        conversion = 1.0;
                    }
                    else if (conversion_read == "Mev-1Fm-3")
                    {
                        conversion = constants::conversion::gev_over_mev / constants::conversion::fm3_gev3;
                    }
                    else
                    {
                        RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for an \"DensityChemPotentialDerivatives\" entry.");
                    }
                }
                else
                    RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for an \"DensityChemPotentialDerivatives\" entry.");
                if (!column_read.is_number_integer())
                    RHM_THROW(std::runtime_error, "UI error: Column number must be provided as an integer for an \"DensityChemPotentialDerivatives\" entry.");
                entry.get() = [column_read, conversion](double nbar)
                {
                    return data_reader({nbar}, column_read) * conversion;
                };
            }
            else
                entry.get() = [](double)
                {
                    return 0.0;
                };
        }
        auto roto_time_units = j["RHSolver"]["TimeUnits"];
        double roto_time_conversion;
        if (roto_time_units.is_null())
            RHM_THROW(std::runtime_error, "UI error: Time units must be provided for rotochemical heating.");
        else if (roto_time_units.is_number())
            roto_time_conversion = roto_time_units.get<double>();
        else if (roto_time_units.is_string())
        {
            if (roto_time_units == "Gev-1")
            {
                roto_time_conversion = 1.0;
            }
            else if (roto_time_units == "Yr")
            {
                roto_time_conversion = 1E-6 * constants::conversion::myr_over_s * constants::conversion::gev_s;
            }
            else if (roto_time_units == "S")
            {
                roto_time_conversion = constants::conversion::gev_s;
            }
            else if (roto_time_units == "Ms")
            {
                roto_time_conversion = 1E-3 * constants::conversion::gev_s;
            }
            else
            {
                RHM_THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for rotochemical heating time.");
            }
        }
        else
            RHM_THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for rotochemical heating time.");
        // rotational 2 omega omega_dot dependency of time
        auto roto_omega_sqr_dot_read = j["RHSolver"]["RotationalOmegaSquareDot"];
        if (!roto_omega_sqr_dot_read.is_array())
            RHM_THROW(std::runtime_error, "UI error: Rotational law and related settings must be provided as an array.");
        else
        {
            auto roto_omega_sqr_dot_provided_as_read = roto_omega_sqr_dot_read[0];
            if (roto_omega_sqr_dot_provided_as_read == "BeyondMagneticDipole")
            {
                auto braking_index_read = roto_omega_sqr_dot_read[1];
                auto p0_read = roto_omega_sqr_dot_read[2],
                     p0_dot_read = roto_omega_sqr_dot_read[3];
                if (!(braking_index_read.is_number()))
                    RHM_THROW(std::runtime_error, "UI error: Braking index for rotational law must be provided as a number.");
                if (!(p0_read.is_number()))
                    RHM_THROW(std::runtime_error, "UI error: Initial rotational period must be provided as a number.");
                if (!(p0_dot_read.is_number()))
                    RHM_THROW(std::runtime_error, "UI error: Initial rotational period derivative must be provided as a number.");

                auto braking_index = braking_index_read.get<double>(),
                     p0 = roto_time_conversion * p0_read.get<double>(),
                     p0_dot = p0_dot_read.get<double>();

                omega_sqr_dot = [braking_index, p0, p0_dot](double t)
                {
                    using constants::scientific::Pi;
                    if (braking_index == 1)
                    {
                        return -8 * Pi * Pi * p0_dot / pow(p0, 3.0) * std::exp(-2 * p0_dot * t / p0);
                    }
                    return -8 * Pi * Pi * p0_dot / pow(p0, 3.0) * pow(1 + (braking_index - 1) * p0_dot * t / p0, (braking_index + 1) / (1 - braking_index));
                };
            }
            else
                RHM_THROW(std::runtime_error, "UI error: Rotational law and related settings must be provided in \"BeyondMagneticDipole\" mode.");
        }
    }
}

#endif