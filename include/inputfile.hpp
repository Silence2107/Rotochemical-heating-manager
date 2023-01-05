#ifndef INPUTFILE_H
#define INPUTFILE_H

#include "../include/auxiliaries.h"
#include "../include/constants.h"
#include "../include/eos_reader.h"

#include <fstream>
#include <functional>
#include <vector>

/// @brief inputfile powering the rotochemical manager
namespace inputfile
{
    // (1) definition of data_reader -- takes input vector and outputs vector of outputs from EoS datafile

    // datafile with EoS data
    std::ifstream fstr("/mnt/d/VSProjects/Rotochemical-heating-manager/data/IST_NS.TXT");

    // ISTNS cached EoS. Only use it if you want to erase cache
    auto ist_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>,
                                              std::vector<double>, const std::vector<double> &, std::ifstream &,
                                              const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &>(eos_reader::predefined::ist_for_ns_cached);

    // data_reader takes input vector and outputs vector of outputs from EoS datafile
    auto data_reader = [](const std::vector<double> &input)
    {
        // linearly interpolated ISTNS
        return eos_reader::eos_data(
            input, [&](const std::vector<double> &input, std::ifstream &fstr)
            { return ist_cached(input, fstr, [&](const std::vector<double> &input, const std::vector<double> &output, double val)
                                { return auxiliaries::interpolate(input, output, auxiliaries::InterpolationMode::kLinear, val, false); }); },
            fstr);
    };

    // (2) EoS additional setup

    // energy density function of baryonic density (units are given by datafile)
    std::function<double(double)> energy_density_of_nbar = [](double nbar)
    { return data_reader({nbar})[1]; };

    // pressure function of baryonic density (units are given by datafile)
    std::function<double(double)> pressure_of_nbar = [](double nbar)
    { return data_reader({nbar})[0]; };

    // baryonic density fraction functions of baryonic density (units are given by datafile)
    std::vector<std::function<double(double)>> Y_i_functions_of_nbar =
        {
            [](double nbar)
            { return data_reader({nbar})[5] / nbar; }, // electron fraction
            [](double nbar)
            { return data_reader({nbar})[6] / nbar; }, // neutron fraction
            [](double nbar)
            { return data_reader({nbar})[7] / nbar; } // proton fraction

    };

    // conversion factors from datafile units to natural units
    double energy_density_conversion = constants::conversion::mev_over_fm3_gev4,
           pressure_conversion = constants::conversion::mev_over_fm3_gev4,
           nbar_conversion = 1.0 / constants::conversion::fm3_gev3;
    /// @brief energy density (in MeV fm^-3) limits in IST. _low and _upp represent limits of EoS itself
    /// while _core_limit represent core boundary
    double edensity_low = 2.2134491254971723E-8,
           edensity_upp = 15892.136580408434,
           edensity_core_limit = 94.214131003471735;
    /// @brief pressure (in MeV fm^-3) limits in IST. _low and _upp represent limits of EoS
    double pressure_low = 5.1065210580102853E-14,
           pressure_upp = 228788.58970172083;
    /// @brief baryonic density (in fm^-3) limits in IST. _low and _upp represent limits of EoS itself
    /// while _core_limit and _crust_limit represent to-from crust boundaries
    double nbar_low = 2.3739996827636742E-11,
           nbar_upp = 2.3189838273277710,
           nbar_crust_limit = 9.9798029952044190E-2,
           nbar_core_limit = 9.9999913289570197E-2;

    // (3) TOV solver setup

    // ISTNS cached EoS interpolator. Only use it if you want to erase cache
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