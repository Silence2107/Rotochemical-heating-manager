#ifndef INSTANTIATOR_H
#define INSTANTIATOR_H

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

/// @brief global data powering RHM
namespace instantiator
{
    // (0) System setup
    auxiliaries::io::LogLevel log_level = auxiliaries::io::kError;

    // (1) EoS setup

    // conversion factors from datafile units to natural units
    double energy_density_conversion = constants::conversion::g_over_cm3_gev4,
           pressure_conversion = constants::conversion::dyne_over_cm2_gev4,
           nbar_conversion = 1.0 / constants::conversion::fm3_gev3;

    // read datafile
    auto table = auxiliaries::io::read_tabulated_file("presupplied/APR4/APR_EOS_Acc_Fe_RHMstandard.dat", {0, 0}, {7, 237});

    // data_reader takes input vector and outputs vector of outputs from EoS datafile
    auto data_reader = auxiliaries::math::CachedFunc<std::vector<auxiliaries::math::CachedInterpolatorWrap>,
                                                     double, const std::vector<double> &, size_t>(
        [](std::vector<auxiliaries::math::CachedInterpolatorWrap> &cache,
           const std::vector<double> &input, size_t index)
        {
            if (cache.empty())
            {
                // fill cache with cached interpolation functions for each column
                for (size_t i = 0; i < table.size(); ++i)
                {
                    auto interpolator_cached = auxiliaries::math::CachedInterpolatorWrap(auxiliaries::math::interpolate_cached);
                    cache.push_back(interpolator_cached);
                }
            }
            // unpack input and convert to datafile units
            double nbar = input[0] / nbar_conversion;
            // return cached interpolation functions, with extrapolation enabled for now
            return cache[index](table[2], table[index], auxiliaries::math::InterpolationMode::kLinear, nbar, true, true);
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
           nbar_sf_shift = 9E-2 * nbar_conversion;
    /// @brief energy density limits in natural units. _low and _upp represent limits of EoS itself <para></para>
    /// while _core_limit represents phase transition boundary
    double edensity_low = energy_density_of_nbar(nbar_low),
           edensity_upp = energy_density_of_nbar(nbar_upp);
    /// @brief pressure limits in natural units. _low and _upp represent limits of EoS itself
    double pressure_low = pressure_of_nbar(nbar_low),
           pressure_upp = pressure_of_nbar(nbar_upp);

    // baryonic density fraction functions of baryonic density (natural units)
    std::map<auxiliaries::phys::Species, std::function<double(double)>> number_densities_of_nbar =
        {
            {constants::species::electron, [](double nbar)
             { return data_reader({nbar}, 3) * nbar; }},
            {constants::species::muon, [](double nbar)
             { return data_reader({nbar}, 4) * nbar; }},
            {constants::species::neutron, [](double nbar)
             { return data_reader({nbar}, 5) * nbar; }},
            {constants::species::proton, [](double nbar)
             { return data_reader({nbar}, 6) * nbar; }}};

    // fermi momentum functions of baryonic density (natural units)
    std::map<auxiliaries::phys::Species, std::function<double(double)>> k_fermi_of_nbar =
        {
            {constants::species::electron, [](double nbar)
             { return pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 3) * nbar, 1.0 / 3); }},
            {constants::species::muon, [](double nbar)
             { return pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 4) * nbar, 1.0 / 3); }},
            {constants::species::neutron, [](double nbar)
             { return pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 5) * nbar, 1.0 / 3); }},
            {constants::species::proton, [](double nbar)
             { return pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 6) * nbar, 1.0 / 3); }}};

    // effective mass functions of baryonic density (natural units)
    std::map<auxiliaries::phys::Species, std::function<double(double)>> m_stars_of_nbar =
        {
            {constants::species::electron, [](double nbar)
             { return sqrt(constants::species::electron.mass() * constants::species::electron.mass() + pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 3) * nbar, 2.0 / 3)); }},
            {constants::species::muon, [](double nbar)
             { return sqrt(constants::species::muon.mass() * constants::species::muon.mass() + pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar}, 4) * nbar, 2.0 / 3)); }},
            {constants::species::neutron, [](double nbar)
             { return data_reader({nbar}, 12) * constants::species::neutron.mass(); }},
            {constants::species::proton, [](double nbar)
             { return data_reader({nbar}, 11) * constants::species::proton.mass(); }}};

    std::function<double(double)> ion_volume_fr = [](double nbar)
    {
        using namespace constants::conversion;
        using namespace constants::scientific;
        auto eta_ion = 4.0 / 3 * Pi * pow(1.1, 3.0) * fm3_gev3 * energy_density_of_nbar(nbar) / constants::species::neutron.mass();
        return eta_ion > 1.0 ? 0.0 : eta_ion;
    };

    // (2) TOV solver setup

    // Cached EoS interpolator. Only use it if you want to erase cache
    auto eos_interpolator_cached = auxiliaries::math::CachedInterpolatorWrap(auxiliaries::math::interpolate_cached);
    // Interpolator used for EoS P/rho
    auto eos_interpolator = [](const std::vector<double> &input, const std::vector<double> &output, double val)
    {
        return eos_interpolator_cached(input, output, auxiliaries::math::InterpolationMode::kLinear, val, false, true);
    };

    // interpolation mode for radial functions (nbar(r), m(r), ...)
    auxiliaries::math::InterpolationMode radial_interp_mode = auxiliaries::math::InterpolationMode::kLinear;

    // nbar(r) cached interpolator. Only use it if you want to erase cache
    auto nbar_interpolator_cached = auxiliaries::math::CachedInterpolatorWrap(auxiliaries::math::interpolate_cached);
    // Interpolator used for nbar(r)
    auto nbar_interpolator = [](const std::vector<double> &input, const std::vector<double> &output, double val)
    {
        return nbar_interpolator_cached(input, output, auxiliaries::math::InterpolationMode::kLinear, val, false, true);
    };

    // EoS linspace discretization
    size_t discr_size_EoS = 1000;

    // TOV adaption limit
    size_t tov_adapt_limit = 20;

    // TOV solver radius step size in GeV
    double radius_step = 0.01 * constants::conversion::km_gev;

    // TOV solver surface pressure in GeV^4
    double surface_pressure = 1E-8 * pressure_upp;

    // TOV solver center pressure in GeV^4
    double center_pressure = 0.2166 * pressure_upp;

    // (3) Cooling solver

    // Cooling solver setup

    // desirable relative accuracy of the cooling solvers
    double cooling_newton_step_eps = 1e-5;

    // maximum number of iterations of the cooling solvers
    size_t cooling_newton_max_iter = 50;

    // desirable relative accuracy of the cooling solvers per time step
    double cooling_max_diff_per_t_step = 0.05;

    // Cooling settings
    double crust_eta = 2.26E-18;

    // Critical phenomena settings

    bool superfluid_p_1s0 = false,
         superfluid_n_3p2 = false,
         superfluid_n_1s0 = false;

    std::function<double(double)> superfluid_p_temp = [](double nbar)
    {
        if (superfluid_p_1s0)
        {
            using namespace auxiliaries::phys;
            using constants::species::proton;
            return critical_temperature(k_fermi_of_nbar.at(proton)(nbar), CriticalTemperatureModel::kCCDK_PS);
        }
        return 0.0;
    };
    std::function<double(double)> superfluid_n_temp = [](double nbar)
    {
        using namespace auxiliaries::phys;
        using constants::species::neutron;
        double k_fermi = k_fermi_of_nbar.at(neutron)(nbar);
        if (superfluid_n_3p2 && superfluid_n_1s0)
        {
            if (nbar <= nbar_sf_shift)
                return critical_temperature(k_fermi, CriticalTemperatureModel::kSFB_NS);
            else
                return critical_temperature(k_fermi, CriticalTemperatureModel::kAO_NT);
        }
        else if (superfluid_n_3p2)
            return critical_temperature(k_fermi, CriticalTemperatureModel::kAO_NT);
        else if (superfluid_n_1s0)
            return critical_temperature(k_fermi, CriticalTemperatureModel::kSFB_NS);
        else
            return 0.0;
    };
    std::function<double(double)> superconduct_q_gap = [](double nbar)
    {
        return 0.0;
    };

    // Evolution settings
    double t_init = 0.0 * constants::conversion::myr_over_s * constants::conversion::gev_s,
           t_end = 1E3 * constants::conversion::myr_over_s * constants::conversion::gev_s,
           base_t_step = 1.0E-18 * constants::conversion::myr_over_s * constants::conversion::gev_s;
    // estimate for the number of time points (is also used for time step expansion, if enabled)
    size_t cooling_n_points_estimate = 1000;
    // initial temperature profile
    auto initial_t_profile_inf = [](double r, double r_ns, const std::function<double(double)> &exp_phi, const std::function<double(double)> &nbar_of_r)
    {
        using namespace constants::conversion;
        return 5E9 * exp_phi(r_ns) / gev_over_k;
    };
    // cooling grid step
    double cooling_radius_step = 10 * radius_step;
    // condition on which to switch to equilibrium cooling
    auto switch_to_equilibrium = [](double t_curr, const std::vector<double> &t_profile)
    {
        // return false; // pure non-equilibrium
        // return true; // pure equilibrium
        // switch at "rather plain profile" and "relevantly late"
        return 1E-5 * constants::conversion::myr_over_s * constants::conversion::gev_s < t_curr &&
               std::abs(t_profile.end()[-2] - t_profile.front()) / t_profile.end()[-2] < 0.01;
    };

    // time step expansion rate (set to 1.0 for constant time step)
    double exp_rate_estim = pow((t_end - t_init) / base_t_step, 1.0 / cooling_n_points_estimate) *
                            pow((pow((t_end - t_init) / base_t_step, 1.0 / cooling_n_points_estimate) - 1), 1.0 / cooling_n_points_estimate);

    // (4) Rotochemical heating setup

    // Bij density matrix dependence on nbar
    std::map<std::pair<auxiliaries::phys::Species, auxiliaries::phys::Species>, std::function<double(double)>> dni_to_dmuj;

    // rotational 2 omega omega_dot dependency of time
    std::function<double(double)> omega_sqr_dot = [](double t)
    {
        size_t braking_index = 3;
        double p0 = 1 * 1E-3 * constants::conversion::gev_s,
               p0_dot = 1.0E-20;
        using constants::scientific::Pi;
        if (braking_index == 1)
        {
            return -8 * Pi * Pi * p0_dot / pow(p0, 3.0) * std::exp(-2 * p0_dot * t / p0);
        }
        return -8 * Pi * Pi * p0_dot / pow(p0, 3.0) * pow(1 + (braking_index - 1) * p0_dot * t / p0, (braking_index + 1) / (1 - braking_index));
    };
}

#endif