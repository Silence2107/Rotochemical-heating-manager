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
    /// @brief energy density limits in datafile units. _low and _upp represent limits of EoS itself <para></para>
    /// while _core_limit represents phase transition boundary
    double edensity_low = 1.000E3,
           edensity_core_limit = 1.5197E14,
           edensity_upp = 7.4456E15;
    /// @brief pressure limits in datafile units. _low and _upp represent limits of EoS itself
    double pressure_low = 1.003E17,
           pressure_upp = 8.3308E36;
    /// @brief baryonic density limits in datafile units. _low and _upp represent limits of EoS itself
    /// while _core_limit, _crust_limit represent phase transition boundaries
    double nbar_low = 6.023E-13,
           nbar_upp = 1.89,
           nbar_core_limit = 9E-2,
           nbar_crust_limit = 2.096E-2;

    // datafile with EoS data
    std::ifstream fstr("/mnt/d/VSProjects/Rotochemical-heating-manager/data/APR_EOS_Acc_Fe.dat");

    // Cached file reader. Only use it if you want to erase cache
    auto reader_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>,
                                                 std::vector<double>, const std::vector<double> &, std::ifstream &,
                                                 const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &>(
        [](std::vector<std::vector<double>> &cache, const std::vector<double> &input, std::ifstream &fstr, const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator)
        {
            /// (barionic density &gt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, electron fraction, muon -//-, neutron -//-, proton -//-, lambda -//-, sigma- -//-, sigma0 -//-, sigma+ -//-, m star proton -//-, m star neutron -//-, m star lambda -//-, m star sigma- -//-, m star sigma0 -//-, m star sigma+ -//-)
            ///	(barionic density &lt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, Acell, Aion, Z, [empty])
            if (cache.empty())
            {
                std::string nextline;
                size_t line_number = 0;
                cache.resize(17); // 17 columns (expected max amount)
                for (size_t cols = 0; cols < cache.size(); ++cols)
                    cache[cols].resize(229); // 229 rows (expected amount)
                for (; std::getline(fstr, nextline) && line_number < 6; ++line_number)
                    ; // skip 6 lines
                while (std::getline(fstr >> std::ws, nextline) && line_number < 235)
                {
                    nextline = auxiliaries::retrieve_cleared_line(nextline); // clear line
                    std::stringstream strstr(nextline);
                    std::string word;
                    for (int i = 0; std::getline(strstr, word, ' '); ++i)
                        cache[i][line_number - 6] = std::stod(word);
                    ++line_number;
                } // read all important data into cache
            }
            std::vector<double> output;
            double nbar = input[0];                 // barionic density (input[0])
            if (nbar > nbar_upp || nbar < nbar_low) // we do not have data beyond these values
                throw std::runtime_error("Data request out of range; Encountered in eos_reader::reader_cached");
            output.resize(cache.size()); // output size
            if (nbar > nbar_core_limit)
            {
                for (int i = 0; i < cache.size(); ++i)
                    output[i] = interpolator(cache[2], cache[i], nbar);
            }
            else if (nbar < nbar_crust_limit)
            {
                for (int i = 0; i < cache.size(); ++i)
                    output[i] = interpolator(cache[2], cache[i], nbar);
            }
            else
            {
                // extract data between i.e. at phase transition;
                // for densities, pressure and baryonic density, we introduce slight slope for monotony; other entries get copypasted depending on nbar
                if (nbar > (nbar_core_limit + nbar_crust_limit) / 2.0)
                {
                    output = std::vector<double>({1.5197E+14, 9.2819E+32, 9.0000E-02, 3.1606E-02, 0.0000E+00, 9.6839E-01, 3.1606E-02, 0, 0, 0, 0, 7.3203E-01, 8.8003E-01, 0, 0, 0, 0});
                    // I choose slopes by hand : split 10%/80%/10%
                    output[0] -= 2.0 * (nbar_core_limit - nbar) / (nbar_core_limit - nbar_crust_limit) * (1.5197E+14 - 3.493E+13) / 10.0;
                    output[1] -= 2.0 * (nbar_core_limit - nbar) / (nbar_core_limit - nbar_crust_limit) * (9.2819E+32 - 7.311E+31) / 10.0;
                    output[2] -= 2.0 * (nbar_core_limit - nbar) / (nbar_core_limit - nbar_crust_limit) * (9.0000E-02 - 2.096E-02) / 10.0;
                }
                else
                {
                    output = std::vector<double>({3.493E+13, 7.311E+31, 2.096E-02, 1127., 124., 26., 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0});
                    // I choose slopes by hand : split 10%/80%/10%
                    output[0] += 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (1.5197E+14 - 3.493E+13) / 10.0;
                    output[1] += 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (9.2819E+32 - 7.311E+31) / 10.0;
                    output[2] += 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (9.0000E-02 - 2.096E-02) / 10.0;
                }
            }
            return output;
        });

    // data_reader takes input vector and outputs vector of outputs from EoS datafile
    auto data_reader = [](const std::vector<double> &input)
    {
        // linearly interpolated EOS
        return auxiliaries::eos_data(
            input, [&](const std::vector<double> &input, std::ifstream &fstr)
            { return reader_cached(input, fstr, [&](const std::vector<double> &input, const std::vector<double> &output, double val)
                                   { return auxiliaries::interpolate(input, output, auxiliaries::InterpolationMode::kLinear, val, false); }); },
            fstr);
    };

    // energy density function of baryonic density (units are given by datafile)
    std::function<double(double)> energy_density_of_nbar = [](double nbar)
    { return data_reader({nbar})[0]; };

    // pressure function of baryonic density (units are given by datafile)
    std::function<double(double)> pressure_of_nbar = [](double nbar)
    { return data_reader({nbar})[1]; };

    // baryonic density fraction functions of baryonic density (units are given by datafile)
    std::map<auxiliaries::Species, std::function<double(double)>> Y_i_functions_of_nbar =
        {
            {constants::scientific::electron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar})[3] : 0.0; }},
            {constants::scientific::muon, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar})[4] : 0.0; }},
            {constants::scientific::neutron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar})[5] : 1.0; }},
            {constants::scientific::proton, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar})[6] : 0.0; }}};

    // fermi momentum functions of baryonic density (GeV units)
    std::map<auxiliaries::Species, std::function<double(double)>> k_fermi_of_nbar =
        {
            {constants::scientific::electron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar})[3] * nbar * nbar_conversion, 1.0 / 3) : 0.0; }},
            {constants::scientific::muon, [](double nbar)
             { return (nbar >= nbar_core_limit) ? pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar})[4] * nbar * nbar_conversion, 1.0 / 3) : 0.0; }},
            {constants::scientific::neutron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar})[5] * nbar * nbar_conversion, 1.0 / 3) : pow(3 * constants::scientific::Pi * constants::scientific::Pi * nbar * nbar_conversion, 1.0 / 3); }},
            {constants::scientific::proton, [](double nbar)
             { return (nbar >= nbar_core_limit) ? pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar})[6] * nbar * nbar_conversion, 1.0 / 3) : 0.0; }}};

    // effective mass functions of baryonic density (GeV units)
    std::map<auxiliaries::Species, std::function<double(double)>> m_stars_of_nbar =
        {
            {constants::scientific::electron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? sqrt(constants::scientific::M_e * constants::scientific::M_e + pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar})[3] * nbar * nbar_conversion, 2.0 / 3)) : constants::scientific::M_e; }},
            {constants::scientific::muon, [](double nbar)
             { return (nbar >= nbar_core_limit) ? sqrt(constants::scientific::M_mu * constants::scientific::M_mu + pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar})[4] * nbar * nbar_conversion, 2.0 / 3)) : constants::scientific::M_mu; }},
            {constants::scientific::neutron, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar})[12] * constants::scientific::M_N : constants::scientific::M_N; }},
            {constants::scientific::proton, [](double nbar)
             { return (nbar >= nbar_core_limit) ? data_reader({nbar})[11] * constants::scientific::M_N : constants::scientific::M_N; }}};

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
                                                           auxiliaries::InterpolationMode, double, bool>(auxiliaries::interpolate_cached);
    // Interpolator used for EoS P(rho)
    auto eos_interpolator = [](const std::vector<double> &input, const std::vector<double> &output, double val)
    {
        return eos_interpolator_cached(input, output, auxiliaries::InterpolationMode::kLinear, val, false);
    };

    // nbar(r) cached interpolator. Only use it if you want to erase cache
    auto nbar_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                            double, const std::vector<double> &, const std::vector<double> &,
                                                            auxiliaries::InterpolationMode, double, bool>(auxiliaries::interpolate_cached);
    // Interpolator used for nbar(r)
    auto nbar_interpolator = [](const std::vector<double> &input, const std::vector<double> &output, double val)
    {
        return nbar_interpolator_cached(input, output, auxiliaries::InterpolationMode::kLinear, val, false);
    };

    // EoS linspace discretization
    size_t discr_size_EoS = 1000;

    // TOV solver radius step size in GeV
    double radius_step = 0.001 * 5E19;

    // TOV solver density step size in GeV^4
    double density_step = 1E-8 * edensity_upp * energy_density_conversion;

    // TOV solver center density in GeV^4
    double center_density = 108.3 / 500.0 * edensity_upp * energy_density_conversion;

    // (3) Cooling solver

    // Cooling solver setup
    auto cooling_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                               double, const std::vector<double> &, const std::vector<double> &,
                                                               auxiliaries::InterpolationMode, double, bool>(auxiliaries::interpolate_cached);
    auto cooling_interpolator = [](const std::vector<double> &x, const std::vector<double> &y, double val)
    {
        return cooling_interpolator_cached(x, y, auxiliaries::InterpolationMode::kCubic, val, false);
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