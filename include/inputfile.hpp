#ifndef INPUTFILE_H
#define INPUTFILE_H

#include "../include/auxiliaries.h"
#include "../include/constants.h"
#include "../include/eos_reader.h"

#include <fstream>
#include <functional>
#include <vector>
#include <cmath>

/// @brief inputfile powering the rotochemical manager
namespace inputfile
{
    // (1) definition of data_reader -- takes input vector and outputs vector of outputs from EoS datafile

    // datafile with EoS data
    std::ifstream fstr("/mnt/d/VSProjects/Rotochemical-heating-manager/data/APR_EOS_Acc_Fe.dat");

    // APR4 cached EoS. Only use it if you want to erase cache
    auto apr_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>,
                                              std::vector<double>, const std::vector<double> &, std::ifstream &,
                                              const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &>(eos_reader::predefined::apr4_cached);

    // data_reader takes input vector and outputs vector of outputs from EoS datafile
    auto data_reader = [](const std::vector<double> &input)
    {
        // linearly interpolated APR4
        return eos_reader::eos_data(
            input, [&](const std::vector<double> &input, std::ifstream &fstr)
            { return apr_cached(input, fstr, [&](const std::vector<double> &input, const std::vector<double> &output, double val)
                                { return auxiliaries::interpolate(input, output, auxiliaries::InterpolationMode::kLinear, val, false); }); },
            fstr);
    };

    // (2) EoS additional setup

    // conversion factors from datafile units to natural units
    double energy_density_conversion = constants::conversion::g_over_cm3_gev4,
           pressure_conversion = constants::conversion::dyne_over_cm2_gev4,
           nbar_conversion = 1.0 / constants::conversion::fm3_gev3;
    /// @brief energy density (in g cm^-3) limits in APR4. _low and _upp represent limits of EoS itself <para></para>
    /// while _core_limit represents phase transition boundary
    double edensity_low = 1.000E3,
           edensity_core_limit = 1.5197E14,
           edensity_upp = 7.4456E15;
    /// @brief pressure (in dyne cm^-2) limits in APR4. _low and _upp represent limits of EoS itself
    double pressure_low = 1.003E17,
           pressure_upp = 8.3308E36;
    /// @brief baryonic density (in fm^-3) limits in APR4. _low and _upp represent limits of EoS itself
    /// while _core_limit and _crust_limit represent phase transition boundaries
    double nbar_low = 6.023E-13,
           nbar_upp = 1.89,
           nbar_core_limit = 9E-2,
           nbar_crust_limit = 2.096E-2; // eos barionic density limits

    // energy density function of baryonic density (units are given by datafile)
    std::function<double(double)> energy_density_of_nbar = [](double nbar)
    { return data_reader({nbar})[0]; };

    // pressure function of baryonic density (units are given by datafile)
    std::function<double(double)> pressure_of_nbar = [](double nbar)
    { return data_reader({nbar})[1]; };

    // baryonic density fraction functions of baryonic density (units are given by datafile)
    std::vector<std::function<double(double)>> Y_i_functions_of_nbar =
        {
            [](double nbar)
            { return (nbar >= nbar_core_limit) ? data_reader({nbar})[3] : 0.0; }, // electron fraction
            [](double nbar)
            { return (nbar >= nbar_core_limit) ? data_reader({nbar})[4] : 0.0; }, // muon fraction
            [](double nbar)
            { return (nbar >= nbar_core_limit) ? data_reader({nbar})[5] : 1.0; }, // neutron fraction
            [](double nbar)
            { return (nbar >= nbar_core_limit) ? data_reader({nbar})[6] : 0.0; } // proton fraction

    };

    // effective mass functions of baryonic density (GeV units)
    std::vector<std::function<double(double)>> m_stars_of_nbar =
        {
            [](double nbar)
            { return (nbar >= nbar_core_limit) ? data_reader({nbar})[12] : constants::scientific::M_N; }, // neutron
            [](double nbar)
            { return (nbar >= nbar_core_limit) ? data_reader({nbar})[11] : constants::scientific::M_N; } // proton

    };

    // fermi momentum functions of baryonic density (GeV units)
    std::vector<std::function<double(double)>> k_fermi_of_nbar =
        {
            [](double nbar)
            { return (nbar >= nbar_core_limit) ? pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar})[5] * nbar * nbar_conversion, 1.0 / 3) : pow(3 * constants::scientific::Pi * constants::scientific::Pi * nbar * nbar_conversion, 1.0 / 3); }, // neutron
            [](double nbar)
            { return (nbar >= nbar_core_limit) ? pow(3 * constants::scientific::Pi * constants::scientific::Pi * data_reader({nbar})[6] * nbar * nbar_conversion, 1.0 / 3) : 0.0; } // proton

    };

    // (3) TOV solver setup

    // APR4 cached EoS interpolator. Only use it if you want to erase cache
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
    double density_step = 0.0001 * edensity_upp * energy_density_conversion;

    // TOV solver center density in GeV^4
    double center_density = 0.2 * edensity_upp * energy_density_conversion;
}

#endif