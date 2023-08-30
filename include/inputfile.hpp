#ifndef INPUTFILE_H
#define INPUTFILE_H

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

/// @brief inputfile powering the rotochemical manager
namespace inputfile
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
    std::map<auxiliaries::phys::Species, std::function<double(double)>> bar_densities_of_nbar;

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

    // TOV solver radius step size in GeV
    double radius_step;

    // TOV solver density step size in GeV^4
    double density_step;

    // TOV solver center density in GeV^4
    double center_density;

    // (3) Cooling solver

    // Cooling solver setup

    // desirable relative accuracy of the cooling solvers
    double cooling_newton_step_eps;

    // maximum number of iterations of the cooling solvers
    size_t cooling_newton_max_iter;

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
    std::function<double(double, double)> initial_t_profile_inf;
    // cooling grid step
    double cooling_radius_step;
    // condition on which to switch to equilibrium cooling
    std::function<bool(double, const std::vector<double> &)> switch_to_equilibrium;

    // time step expansion rate (set to 1.0 for constant time step)
    double exp_rate_estim;

    /// @brief instantiate the inputfile from json input
    /// @param json_input json inputfile path
    void instantiate_system(const std::string &json_input)
    {
        using namespace inputfile;
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
                THROW(std::runtime_error, "UI error: Unparsable interpolation mode provided.");
        };

        // read json input
        std::ifstream i(json_input);
        if (!i.is_open())
        {
            THROW(std::runtime_error, "UI inputfile requested, but the path is invalid.");
        }
        json j = json::parse(i);

        // (1) EoS setup

        // filereader
        auto eos_datafile = j["EoSSetup"]["Datafile"]["Path"];
        auto eos_datafile_rows = j["EoSSetup"]["Datafile"]["Rows"];
        if (!(eos_datafile_rows.size() == 2))
        {
            if (eos_datafile_rows.is_null())
                eos_datafile_rows = {0, 0};
            else
                THROW(std::runtime_error, "UI error: Datafile rows must be a pair-array.");
        }

        auto eos_datafile_cols = j["EoSSetup"]["Datafile"]["Columns"];
        if (!(eos_datafile_cols.size() == 2))
        {
            if (eos_datafile_cols.is_null())
                eos_datafile_cols = {0, 0};
            else
                THROW(std::runtime_error, "UI error: Datafile cols must be a pair-array.");
        }

        auxiliaries::math::InterpolationMode eos_datafile_interp_mode;
        auto eos_datafile_interp_read = j["EoSSetup"]["Datafile"]["Interpolation"];
        if (eos_datafile_interp_read.is_null())
            eos_datafile_interp_mode = auxiliaries::math::InterpolationMode::kLinear;
        else if (!(eos_datafile_interp_read.is_string()))
            THROW(std::runtime_error, "UI error: Datafile interpolation mode must be a string.");
        else
            eos_datafile_interp_mode = get_interpolation_mode(eos_datafile_interp_read);

        auto nbar_index = j["EoSSetup"]["Quantities"]["BarionicDensity"]["Column"];
        auto nbar_conversion_read = j["EoSSetup"]["Quantities"]["BarionicDensity"]["Units"];
        if (!(nbar_index.is_number_integer()))
            THROW(std::runtime_error, "UI error: Barionic density column number must be provided as an integer.");
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
                THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for barionic density.");
            }
        }
        else
            THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for barionic density.");

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
                return cache[index](table[nbar_index], table[index], eos_datafile_interp_mode, nbar, true, true);
            });

        // energy density function of baryonic density (natural units)
        auto energy_density_index = j["EoSSetup"]["Quantities"]["EnergyDensity"]["Column"];
        auto energy_density_conversion_read = j["EoSSetup"]["Quantities"]["EnergyDensity"]["Units"];
        if (!(energy_density_index.is_number_integer()))
            THROW(std::runtime_error, "UI error: Energy density column number must be provided as an integer.");
        if (energy_density_conversion_read.is_number())
            energy_density_conversion = energy_density_conversion_read.get<double>();
        else if (energy_density_conversion_read.is_string())
        {
            if (energy_density_conversion_read == "Gev4")
            {
                energy_density_conversion = 1.0;
            }
            else if (energy_density_conversion_read == "MevOverFm3")
            {
                energy_density_conversion = constants::conversion::mev_over_fm3_gev4;
            }
            else if (energy_density_conversion_read == "GOverCm3")
            {
                energy_density_conversion = constants::conversion::g_over_cm3_gev4;
            }
            else
            {
                THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for energy density.");
            }
        }
        else
            THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for energy density.");
        energy_density_of_nbar = [energy_density_index](double nbar)
        { return data_reader({nbar}, energy_density_index) * energy_density_conversion; };

        // pressure function of baryonic density (natural units)
        auto pressure_index = j["EoSSetup"]["Quantities"]["Pressure"]["Column"];
        auto pressure_conversion_read = j["EoSSetup"]["Quantities"]["Pressure"]["Units"];
        if (!(pressure_index.is_number_integer()))
            THROW(std::runtime_error, "UI error: Pressure column number must be provided as an integer.");
        if (pressure_conversion_read.is_number())
            pressure_conversion = pressure_conversion_read.get<double>();
        else if (pressure_conversion_read.is_string())
        {
            if (pressure_conversion_read == "Gev4")
            {
                pressure_conversion = 1.0;
            }
            else if (pressure_conversion_read == "MevOverFm3")
            {
                pressure_conversion = constants::conversion::mev_over_fm3_gev4;
            }
            else if (pressure_conversion_read == "DyneOverCm2")
            {
                pressure_conversion = constants::conversion::dyne_over_cm2_gev4;
            }
            else
            {
                THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for pressure.");
            }
        }
        else
            THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for pressure.");
        pressure_of_nbar = [pressure_index](double nbar)
        { return data_reader({nbar}, pressure_index) * pressure_conversion; };

        // baryonic density limits in natural units. _low and _upp represent limits of EoS itself
        // while _core_limit, _crust_limit represent phase transition boundaries
        auto nbar_low_read = j["EoSSetup"]["Quantities"]["BarionicDensity"]["Low"],
             nbar_upp_read = j["EoSSetup"]["Quantities"]["BarionicDensity"]["Upp"],
             nbar_core_limit_read = j["EoSSetup"]["Quantities"]["BarionicDensity"]["CoreLimit"],
             nbar_crust_limit_read = j["EoSSetup"]["Quantities"]["BarionicDensity"]["CrustLimit"];
        if (!(nbar_low_read.is_number() && nbar_upp_read.is_number() && nbar_core_limit_read.is_number()))
            THROW(std::runtime_error, "UI error: Barionic density limits must be provided as numbers.");

        if (!nbar_crust_limit_read.is_number() && !nbar_crust_limit_read.is_null())
            THROW(std::runtime_error, "UI error: Barionic density crust limit may only be provided as a number.");
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
            THROW(std::runtime_error, "UI error: Energy density low limit may only be provided as a number or \"Deduce\".");

        if (edensity_core_limit_read.is_null() || edensity_core_limit_read == "Deduce")
            edensity_core_limit = energy_density_of_nbar(nbar_core_limit);
        else if (edensity_core_limit_read.is_number())
            edensity_core_limit = edensity_core_limit_read.get<double>() * energy_density_conversion;
        else
            THROW(std::runtime_error, "UI error: Energy density core limit may only be provided as a number or \"Deduce\".");

        if (edensity_upp_read.is_null() || edensity_upp_read == "Deduce")
            edensity_upp = energy_density_of_nbar(nbar_upp);
        else if (edensity_upp_read.is_number())
            edensity_upp = edensity_upp_read.get<double>() * energy_density_conversion;
        else
            THROW(std::runtime_error, "UI error: Energy density upp limit may only be provided as a number or \"Deduce\".");

        // pressure limits in natural units. _low and _upp represent limits of EoS itself
        auto pressure_low_read = j["EoSSetup"]["Quantities"]["Pressure"]["Low"],
             pressure_upp_read = j["EoSSetup"]["Quantities"]["Pressure"]["Upp"];

        if (pressure_low_read.is_null() || pressure_low_read == "Deduce")
            pressure_low = pressure_of_nbar(nbar_low);
        else if (pressure_low_read.is_number())
            pressure_low = pressure_low_read.get<double>() * pressure_conversion;
        else
            THROW(std::runtime_error, "UI error: Pressure low limit may only be provided as a number or \"Deduce\".");

        if (pressure_upp_read.is_null() || pressure_upp_read == "Deduce")
            pressure_upp = pressure_of_nbar(nbar_upp);
        else if (pressure_upp_read.is_number())
            pressure_upp = pressure_upp_read.get<double>() * pressure_conversion;
        else
            THROW(std::runtime_error, "UI error: Pressure upp limit may only be provided as a number or \"Deduce\".");

        // TOV solver setup

        auxiliaries::math::InterpolationMode eos_interp_mode;
        auto eos_interp_mode_read = j["TOVSolver"]["EoSInterpolation"];
        if (eos_interp_mode_read.is_null())
            eos_interp_mode = auxiliaries::math::InterpolationMode::kLinear;
        else if (!(eos_interp_mode_read.is_string()))
            THROW(std::runtime_error, "UI error: EoS interpolation mode must be a string.");
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
            THROW(std::runtime_error, "UI error: nbar interpolation mode must be a string.");
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
            THROW(std::runtime_error, "UI error: EoS discretization must be provided as an integer.");
        else
            discr_size_EoS = discr_size_EoS_read.get<size_t>();

        // TOV solver radius step size in GeV
        auto radius_step_read = j["TOVSolver"]["RadiusStep"]; // reads in km
        if (radius_step_read.is_null())
            radius_step = 0.01 * constants::conversion::km_gev; // default value is 10 meters
        else if (!(radius_step_read.is_number()))
            THROW(std::runtime_error, "UI error: TOV solver radius step must be provided as a number.");
        else
            radius_step = radius_step_read.get<double>() * constants::conversion::km_gev;

        // TOV solver density step size in GeV^4
        auto density_step_read = j["TOVSolver"]["DensityStep"]; // reads in fraction of maxima
        if (density_step_read.is_null())
            density_step = 1E-8 * edensity_upp; // default value is 1E-8 of maxima
        else if (!(density_step_read.is_number()))
            THROW(std::runtime_error, "UI error: TOV solver density step must be provided as a number.");
        else
            density_step = density_step_read.get<double>() * edensity_upp;

        // TOV solver center density in GeV^4
        auto center_density_read = j["TOVSolver"]["CenterDensity"]; // reads in fraction of maxima
        if (!(center_density_read.is_number()))
            THROW(std::runtime_error, "UI error: TOV solver center density must be provided as a number.");
        else
            center_density = center_density_read.get<double>() * edensity_upp;

        // provided particles
        auto particles_read = j["EoSSetup"]["Particles"];
        if (particles_read.size() < 1)
            return; // no particles provided -> user expresses no interest in cooling
        if (!particles_read.is_array())
            THROW(std::runtime_error, "UI error: Particle types must be provided as an array.");

        std::vector<auxiliaries::phys::Species> particles;
        for (auto particle_name : particles_read)
        {
            using namespace constants::species;
            std::vector<auxiliaries::phys::Species> known_particles =
                {neutron, proton, electron, muon, tau, uquark, dquark, squark};
            if (particle_name.is_string())
            {
                // identify particle with the list of known ones
                for (const auto &known_particle : known_particles)
                {
                    if (particle_name == known_particle.name())
                    {
                        particles.push_back(known_particle);
                        break;
                    }
                }
            }
            else
                THROW(std::runtime_error, "UI error: Particle type must be provided as a string.");
        }

        // baryonic density functions of baryonic density (natural units) &&
        // fermi momentum functions of baryonic density (natural units)
        auto particle_density_conversion_read = j["EoSSetup"]["Quantities"]["BarionicDensities"]["Units"];
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
                THROW(std::runtime_error, "UI error: Unexpected conversion unit provided for particle density.");
            }
        }
        else
            THROW(std::runtime_error, "UI error: Unparsable conversion unit provided for particle density.");
        for (const auto &particle : particles)
        {
            auto particle_name = particle.name();
            auto particle_density_read = j["EoSSetup"]["Quantities"]["BarionicDensities"][particle_name];
            
            // I for now assume that only mesons are not fermions, which is of course very brave. Consider TODO
            if (particle_density_read.is_null() && !(particle.classify() == auxiliaries::phys::Species::ParticleClassification::kMeson))
                THROW(std::runtime_error, "UI error: Particle density info is not provided for " + particle.name() + ".");
            
            auto particle_density_column_read = particle_density_read["Column"];
            if (!(particle_density_column_read.is_number_integer()))
                THROW(std::runtime_error, "UI error: Particle density column number must be provided as an integer.");
            
            auto particle_density_provided_as_read = particle_density_read["ProvidedAs"];
            if (particle_density_provided_as_read.is_null())
                particle_density_provided_as_read = "Density";

            bar_densities_of_nbar.insert(
                {particle, [particle_density_column_read, particle_density_provided_as_read, particle_density_conversion](double nbar)
                 {
                     if (particle_density_provided_as_read == "Density")
                     {
                         return data_reader({nbar}, particle_density_column_read) * particle_density_conversion;
                     }
                     else if (particle_density_provided_as_read == "DensityFraction")
                     {
                         return data_reader({nbar}, particle_density_column_read) * particle_density_conversion * nbar;
                     }
                     else if (particle_density_provided_as_read == "KFermi")
                     {
                         using constants::scientific::Pi;
                         return pow(data_reader({nbar}, particle_density_column_read) * particle_density_conversion, 3.0) / (3.0 * Pi * Pi);
                     }
                     else
                        THROW(std::runtime_error, "UI error: Particle density may only be provided in \"Density\", \"DensityFraction\" or \"KFermi\" modes.");
                 }});
            k_fermi_of_nbar.insert(
                {particle, [particle](double nbar) 
                 {
                    using constants::scientific::Pi;
                    return pow(3.0 * Pi * Pi * bar_densities_of_nbar[particle](nbar), 1.0 / 3.0);
                 }});
        }
    }
}

#endif