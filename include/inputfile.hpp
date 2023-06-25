#ifndef INPUTFILE_H
#define INPUTFILE_H

#include "../include/auxiliaries.h"
#include "../include/constants.h"
#include "../include/cooling.h"

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
    double energy_density_conversion = constants::conversion::g_over_cm3_gev4,
           pressure_conversion = constants::conversion::dyne_over_cm2_gev4,
           nbar_conversion = 1.0 / constants::conversion::fm3_gev3;

    // read datafile
    auto table = auxiliaries::read_tabulated_file("/mnt/d/VSProjects/Rotochemical-heating-manager/data/APR_EOS_Acc_Fe.dat", {0, 0}, {7, 236});

    // data_reader takes input vector and outputs vector of outputs from EoS datafile
    auto data_reader = auxiliaries::CachedFunc<std::vector<auxiliaries::CachedFunc<std::function<double(double)>,
                                                                                   double, const std::vector<double> &, const std::vector<double> &,
                                                                                   auxiliaries::InterpolationMode, double, bool, bool>>,
                                               double, const std::vector<double> &, size_t>(
        [](std::vector<auxiliaries::CachedFunc<std::function<double(double)>,
                                               double, const std::vector<double> &, const std::vector<double> &,
                                               auxiliaries::InterpolationMode, double, bool, bool>> &cache,
           const std::vector<double> &input, size_t index)
        {
            if (cache.empty())
            {
                // fill cache with cached interpolation functions for each column
                for (size_t i = 0; i < table.size(); ++i)
                {
                    auto interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                                       double, const std::vector<double> &, const std::vector<double> &,
                                                                       auxiliaries::InterpolationMode, double, bool, bool>(auxiliaries::interpolate_cached);
                    cache.push_back(interpolator_cached);
                }
            }
            // unpack input and convert to datafile units
            double nbar = input[0] / nbar_conversion;
            // return cached interpolation functions, with extrapolation enabled for now
            return cache[index](table[2], table[index], auxiliaries::InterpolationMode::kLinear, nbar, true, true);
        });

    // energy density function of baryonic density (natural units)
    std::function<double(double)> energy_density_of_nbar = [](double nbar)
    { return data_reader({nbar}, 0) * energy_density_conversion; };

    // pressure function of baryonic density (natural units)
    std::function<double(double)> pressure_of_nbar = [](double nbar)
    { return data_reader({nbar}, 1) * pressure_conversion; };

    /// @brief baryonic density limits in natural units. _low and _upp represent limits of EoS itself
    /// while _core_limit, _crust_limit represent phase transition boundaries
    double nbar_low = 6.023E-13 * nbar_conversion,
           nbar_upp = 1.89 * nbar_conversion,
           nbar_core_limit = 9E-2 * nbar_conversion,
           nbar_crust_limit = 2.096E-2 * nbar_conversion;
    /// @brief energy density limits in natural units. _low and _upp represent limits of EoS itself <para></para>
    /// while _core_limit represents phase transition boundary
    double edensity_low = energy_density_of_nbar(nbar_low),
           edensity_core_limit = energy_density_of_nbar(nbar_core_limit),
           edensity_upp = energy_density_of_nbar(nbar_upp);
    /// @brief pressure limits in natural units. _low and _upp represent limits of EoS itself
    double pressure_low = pressure_of_nbar(nbar_low),
           pressure_upp = pressure_of_nbar(nbar_upp);

    // baryonic density fraction functions of baryonic density (natural units)
    std::map<auxiliaries::Species, std::function<double(double)>> Y_i_functions_of_nbar =
        {
            {constants::scientific::electron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar}, 3) : 0.0; }},
            {constants::scientific::muon, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar}, 4) : 0.0; }},
            {constants::scientific::neutron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar}, 5) : 1.0; }},
            {constants::scientific::proton, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar}, 6) : 0.0; }}};

    // fermi momentum functions of baryonic density (natural units)
    std::map<auxiliaries::Species, std::function<double(double)>> k_fermi_of_nbar =
        {
            {constants::scientific::electron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 3) * nbar, 1.0 / 3) : 0.0; }},
            {constants::scientific::muon, [](double nbar)
             { return (nbar >= nbar_core_limit) ? pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 4) * nbar, 1.0 / 3) : 0.0; }},
            {constants::scientific::neutron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 5) * nbar, 1.0 / 3) : pow(3 * constants::scientific::Pi * constants::scientific::Pi * nbar, 1.0 / 3); }},
            {constants::scientific::proton, [](double nbar)
             { return (nbar >= nbar_core_limit) ? pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 6) * nbar, 1.0 / 3) : 0.0; }}};

    // effective mass functions of baryonic density (natural units)
    std::map<auxiliaries::Species, std::function<double(double)>> m_stars_of_nbar =
        {
            {constants::scientific::electron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? sqrt(constants::scientific::M_e * constants::scientific::M_e + pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 3) * nbar, 2.0 / 3)) : constants::scientific::M_e; }},
            {constants::scientific::muon, [](double nbar)
             { return (nbar >= nbar_core_limit) ? sqrt(constants::scientific::M_mu * constants::scientific::M_mu + pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 4) * nbar, 2.0 / 3)) : constants::scientific::M_mu; }},
            {constants::scientific::neutron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar}, 12) * constants::scientific::M_N : constants::scientific::M_N; }},
            {constants::scientific::proton, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar}, 11) * constants::scientific::M_N : constants::scientific::M_N; }}};

    std::function<double(double)> ion_volume_fr = [](double nbar)
    {
        if (nbar >= nbar_crust_limit)
            return 0.0;
        using namespace constants::conversion;
        using namespace constants::scientific;
        auto eta_ion = 4.0 / 3 * Pi * pow(1.1, 3.0) * fm3_gev3 * energy_density_of_nbar(nbar) / M_N;
        return std::min(1.0, eta_ion);
    };

    // (2) TOV solver setup

    // Cached EoS interpolator. Only use it if you want to erase cache
    auto eos_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                           double, const std::vector<double> &, const std::vector<double> &,
                                                           auxiliaries::InterpolationMode, double, bool, bool>(auxiliaries::interpolate_cached);
    // Interpolator used for EoS P(rho)
    auto eos_interpolator = [](const std::vector<double> &input, const std::vector<double> &output, double val)
    {
        return eos_interpolator_cached(input, output, auxiliaries::InterpolationMode::kLinear, val, false, true);
    };

    // nbar(r) cached interpolator. Only use it if you want to erase cache
    auto nbar_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                            double, const std::vector<double> &, const std::vector<double> &,
                                                            auxiliaries::InterpolationMode, double, bool, bool>(auxiliaries::interpolate_cached);
    // Interpolator used for nbar(r)
    auto nbar_interpolator = [](const std::vector<double> &input, const std::vector<double> &output, double val)
    {
        return nbar_interpolator_cached(input, output, auxiliaries::InterpolationMode::kLinear, val, false, true);
    };

    // EoS linspace discretization
    size_t discr_size_EoS = 1000;

    // TOV solver radius step size in GeV
    double radius_step = 0.001 * 5E19;

    // TOV solver density step size in GeV^4
    double density_step = 1E-8 * edensity_upp;

    // TOV solver center density in GeV^4
    double center_density = 108.3 / 500.0 * edensity_upp;

    // (3) Cooling solver

    // Cooling solver setup
    auto cooling_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                               double, const std::vector<double> &, const std::vector<double> &,
                                                               auxiliaries::InterpolationMode, double, bool, bool>(auxiliaries::interpolate_cached);
    auto cooling_interpolator = [](const std::vector<double> &x, const std::vector<double> &y, double val)
    {
        return cooling_interpolator_cached(x, y, auxiliaries::InterpolationMode::kCubic, val, false, true);
    };

    // Cooling settings
    double crust_eta = 2.26E-18;

    // Critical phenomena settings

    bool superfluid_p_1s0 = true,
         superfluid_n_3p2 = true,
         superfluid_n_1s0 = true;

    std::function<double(double)> superfluid_p_temp = [](double k_fermi)
    {
        if (superfluid_p_1s0)
        {
            using namespace cooling::predefined::auxiliary;
            return critical_temperature(k_fermi, CriticalTemperatureModel::kCCDK);
        }
        return 0.0;
    };
    std::function<double(double)> superfluid_n_temp = [](double k_fermi)
    {
        if (superfluid_n_3p2 || superfluid_n_1s0)
        {
            using namespace cooling::predefined::auxiliary;
            return critical_temperature(k_fermi, CriticalTemperatureModel::kCCDK);
        }
        return 0.0;
    };
    std::function<double(double)> superconduct_q_gap = [](double nbar)
    {
        return 0.0;
    };

    // Evolution settings
    double t_init = 0.0 * constants::conversion::myr_over_s * constants::conversion::gev_s,
           t_end = 6.1E-0 * constants::conversion::myr_over_s * constants::conversion::gev_s,
           base_t_step = 1.0E-18 * constants::conversion::myr_over_s * constants::conversion::gev_s;
    // estimate for the number of time points (used for time step expansion, if enabled)
    double cooling_n_points_estimate = 1000;
    double T_init_local = 5E9 / constants::conversion::gev_over_k;

    // time step expansion rate (set to 1.0 for constant time step)
    double exp_rate_estim = pow((t_end - t_init) / base_t_step, 1.0 / cooling_n_points_estimate) *
                            pow((pow((t_end - t_init) / base_t_step, 1.0 / cooling_n_points_estimate) - 1), 1.0 / cooling_n_points_estimate);
}

#endif