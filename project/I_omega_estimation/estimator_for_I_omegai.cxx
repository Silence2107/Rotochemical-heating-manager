
#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "../include/eos_reader.h"  // allows to read EoS datafiles
#include "../include/tov_solver.h"  // contains TOV solver
#include "../include/constants.h"   // contains constants
#include "../include/auxiliaries.h" // contains auxiliary functionality

int main(int argc, char **argv)
{
    if (argc != 2)
    {
        std::cerr << "Usage: " << argv[0] << " <center density fraction>" << std::endl;
        return 1;
    }
    // Setup that calculates I_{\omega_i} (based on provided central_density/max_eos_density argv[1]), which is the following quantity:
    // I_{\omega i} = \int_{core} dY_i/dP * dP/d\omega^2 * n * dV, with estimation dP/d\omega^2 \approx -P/\omega_K^2
    // where \omega_K is the Keplerian frequency and i being particle species

    // SETUP ------------------------------------------------------------------------

    //IST NS SETUP
    // (1) definition of data_reader -- takes input vector and outputs vector of outputs from EoS datafile
    std::ifstream fstr("../../data/IST_NS.TXT");
    auto ist_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>,
                                              std::vector<double>, const std::vector<double> &, std::ifstream &,
                                              const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &>(eos_reader::predefined::ist_for_ns_cached);
    auto data_reader = [&fstr, &ist_cached](const std::vector<double> &input)
    {
        // linearly interpolated IST NS
        return eos_reader::eos_data(
            input, [&](const std::vector<double> &input, std::ifstream &fstr)
            { return ist_cached(input, fstr, [&](const std::vector<double> &input, const std::vector<double> &output, double val)
                                { return auxiliaries::interpolate(input, output, auxiliaries::InterpolationMode::kLinear, val, false); }); },
            fstr);
    };

    // (2) EoS additional setup
    size_t output_index_energy_density = 1, // index of energy density in EoS output
        output_index_pressure = 0;          // index of pressure in EoS output
    std::vector<std::function<double(double)>> Y_i_functions_of_nbar =
        {
            [&data_reader](double nbar)
            { return data_reader({nbar})[5] / nbar; }, // electron fraction
            [&data_reader](double nbar)
            { return data_reader({nbar})[6] / nbar; }, // neutron fraction
            [&data_reader](double nbar)
            { return data_reader({nbar})[7] / nbar; }  // proton fraction

        };                                                                       // Species fractions as a function of baryonic density
    double energy_density_conversion = constants::conversion::mev_over_fm3_gev4, // conversion factor from EoS energy density to natural units
        pressure_conversion = constants::conversion::mev_over_fm3_gev4,          // conversion factor from EoS pressure to natural units
        nbar_conversion = 1.0 / constants::conversion::fm3_gev3;                 // conversion factor from EoS nbar to natural units
    double rho_min = constants::ist_ns::edensity_low,
           rho_core = constants::ist_ns::edensity_core_limit,
           rho_max = constants::ist_ns::edensity_upp; // eos density limits
    double p_min = constants::ist_ns::pressure_low,
           p_max = constants::ist_ns::pressure_upp; // eos pressure limits
    double nbar_min = constants::ist_ns::nbar_low,
           nbar_max = constants::ist_ns::nbar_upp; // eos barionic density limits

    // (3) TOV solver setup
    auto eos_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                           double, const std::vector<double> &, const std::vector<double> &,
                                                           auxiliaries::InterpolationMode, double, bool>(auxiliaries::interpolate_cached);
    auto eos_interpolator = [&eos_interpolator_cached](const std::vector<double> &input, const std::vector<double> &output, double val)
    {
        // P(rho) cached cubic interpolator, checks enabled
        return eos_interpolator_cached(input, output, auxiliaries::InterpolationMode::kLinear, val, false);
    }; // defines EoS interpolator
    auto nbar_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                            double, const std::vector<double> &, const std::vector<double> &,
                                                            auxiliaries::InterpolationMode, double, bool>(auxiliaries::interpolate_cached);
    auto nbar_interpolator = [&nbar_interpolator_cached](const std::vector<double> &input, const std::vector<double> &output, double val)
    {
        // nbar(r) cached cubic interpolator, checks enabled
        return nbar_interpolator_cached(input, output, auxiliaries::InterpolationMode::kLinear, val, false);
    };                                                                  // defines nbar interpolator
    size_t discr_size_EoS = 1000;                                       // linspace size for EoS interpolation
    double radius_step = 0.001 * 5E19;                                  // TOV solver radius step size in GeV
    double density_step = 0.0001 * rho_max * energy_density_conversion; // TOV solver density step size in GeV^4
    double center_density = std::stod(argv[1]) * rho_max * energy_density_conversion;  // TOV solver center density in GeV^4
    

    /* // APR4 SETUP
    // (1) definition of data_reader -- takes input vector and outputs vector of outputs from EoS datafile
    std::ifstream fstr("../data/APR_EOS_Acc_Fe.dat");
    auto apr_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>,
                                              std::vector<double>, const std::vector<double> &, std::ifstream &,
                                              const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &>(eos_reader::predefined::apr4_cached);
    auto data_reader = [&fstr, &apr_cached](const std::vector<double> &input)
    {
        // linearly interpolated IST NS
        return eos_reader::eos_data(
            input, [&](const std::vector<double> &input, std::ifstream &fstr)
            { return apr_cached(input, fstr, [&](const std::vector<double> &input, const std::vector<double> &output, double val)
                                { return auxiliaries::interpolate(input, output, auxiliaries::InterpolationMode::kLinear, val, false); }); },
            fstr);
    };

    // (2) EoS additional setup
    size_t output_index_energy_density = 0, // index of energy density in EoS output
        output_index_pressure = 1;          // index of pressure in EoS output
    std::vector<std::function<double(double)>> Y_i_functions_of_nbar =
        {
            [&data_reader](double nbar)
            { return data_reader({nbar})[3]; }, // electron fraction
            [&data_reader](double nbar)
            { return data_reader({nbar})[4]; }, // muon fraction
            [&data_reader](double nbar)
            { return data_reader({nbar})[5]; }, // neutron fraction
            [&data_reader](double nbar)
            { return data_reader({nbar})[6]; } // proton fraction

        };                                                                     // Species fractions as a function of baryonic density
    double energy_density_conversion = constants::conversion::g_over_cm3_gev4, // conversion factor from EoS energy density to natural units
        pressure_conversion = constants::conversion::dyne_over_cm2_gev4,       // conversion factor from EoS pressure to natural units
        nbar_conversion = 1.0 / constants::conversion::fm3_gev3;               // conversion factor from EoS nbar to natural units
    double rho_min = constants::apr4::edensity_low,
           rho_core = constants::apr4::edensity_core_limit,
           rho_max = constants::apr4::edensity_upp; // eos density limits
    double p_min = constants::apr4::pressure_low,
           p_max = constants::apr4::pressure_upp; // eos pressure limits
    double nbar_min = constants::apr4::nbar_low,
           nbar_max = constants::apr4::nbar_upp; // eos barionic density limits

    // (3) TOV solver setup
    auto eos_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                           double, const std::vector<double> &, const std::vector<double> &,
                                                           auxiliaries::InterpolationMode, double, bool>(auxiliaries::interpolate_cached);
    auto eos_interpolator = [&eos_interpolator_cached](const std::vector<double> &input, const std::vector<double> &output, double val)
    {
        // P(rho) cached cubic interpolator, checks enabled
        return eos_interpolator_cached(input, output, auxiliaries::InterpolationMode::kLinear, val, false);
    }; // defines EoS interpolator
    auto nbar_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                            double, const std::vector<double> &, const std::vector<double> &,
                                                            auxiliaries::InterpolationMode, double, bool>(auxiliaries::interpolate_cached);
    auto nbar_interpolator = [&nbar_interpolator_cached](const std::vector<double> &input, const std::vector<double> &output, double val)
    {
        // nbar(r) cached cubic interpolator, checks enabled
        return nbar_interpolator_cached(input, output, auxiliaries::InterpolationMode::kLinear, val, false);
    };                                                                                // defines nbar interpolator
    size_t discr_size_EoS = 1000;                                                     // linspace size for EoS interpolation
    double radius_step = 0.001 * 5E19;                                                // TOV solver radius step size in GeV
    double density_step = 0.0001 * rho_max * energy_density_conversion;               // TOV solver density step size in GeV^4
    double center_density = std::stod(argv[1]) * rho_max * energy_density_conversion; // TOV solver center density in GeV^4
    */

    // RUN --------------------------------------------------------------------------

    // EoS definition

    auto eos_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>, double, double>(
        [&](std::vector<std::vector<double>> &cache, double rho)
        {
            if (rho < 0 || rho > rho_max * energy_density_conversion)
                throw std::runtime_error("Data request out of range; Encountered in main::eos_cached");
            if (rho <= rho_min * energy_density_conversion)
                return 0.0;
            if (cache.empty() || cache[0].size() != discr_size_EoS)
            {                                                                                        // then fill/refill cache
                cache = std::vector<std::vector<double>>(2, std::vector<double>(discr_size_EoS, 0)); // initialize 2xdiscr_size_EoS matrix
                std::vector<double> x(discr_size_EoS, 0);
                for (int i = 1; i < discr_size_EoS - 1; ++i)
                { // cache EoS for further efficiency
                    x[i] = i * (nbar_max - nbar_min) / discr_size_EoS + nbar_min;
                    cache[0][i] = energy_density_conversion * data_reader(std::vector<double>({x[i]}))[output_index_energy_density];
                    cache[1][i] = pressure_conversion * data_reader(std::vector<double>({x[i]}))[output_index_pressure];
                }
                x[0] = nbar_min;
                x[x.size() - 1] = nbar_max;
                cache[0][0] = energy_density_conversion * rho_min;
                cache[0][cache[0].size() - 1] = energy_density_conversion * rho_max;
                cache[1][0] = pressure_conversion * p_min;
                cache[1][cache[1].size() - 1] = pressure_conversion * p_max;
                eos_interpolator_cached.erase(); // clean up cached interpolator
            }
            return eos_interpolator(cache[0], cache[1], rho);
        });

    // TOV solver

    auto tov_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                              const std::function<double(double)> &, double, double, double, double>(tov_solver::tov_solution);
    auto tov = [&tov_cached, &eos_cached, radius_step, density_step, center_density](double r)
    {
        // TOV solution cached
        return tov_cached(eos_cached, r, center_density, radius_step, density_step);
    };

    // we also need nbar(r) so that to map data_reader output to radius
    auto nbar = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>, double>(
        [&](std::vector<std::vector<double>> &cache, double r)
        {
            // cache contains {r, n_B(r)} arrays; recaching is not supported at the moment, call ::erase instead
            // returns [0] -> nbar at given point and [1] -> radius at which crust begins
            // warning : initial_density should not be less (or even of order of) that density at core_limit; exception would be thrown otherwise

            if (cache.empty())
            {
                nbar_interpolator_cached.erase(); // clean up cached interpolator
                double R_ns = tov(0.0)[4];
                cache = std::vector<std::vector<double>>(2, std::vector<double>());
                for (double r_current = 0; r_current < R_ns; r_current += radius_step)
                    cache[0].push_back(r_current);
                for (size_t i = 0; i < cache[0].size(); ++i)
                {
                    double r_current = cache[0][i];
                    // now we somehow have to find corresponding n_B
                    // let's stick to densities
                    double density_at_r = tov(r_current)[1];
                    // density_prec = TMath::Abs(tov(eos, r_current + radius_step, initial_density, radius_step, density_step)[1] - density_at_r); // we need these for bisection; maybe better use density_step
                    if (density_at_r <= energy_density_conversion * rho_core)
                    { // until we reach core limit
                        // cache[1][i] = 0.0;
                        // continue;
                        break; // finish writing
                    }
                    double nbar_low = nbar_min, nbar_upp = nbar_max; // we need these for bisection search; in fm-3 units for now
                    double nbar_mid = (nbar_low + nbar_upp) / 2.0;
                    // while (conversion::g_over_cm3_gev4 * TMath::Abs(istdat(std::vector<double>({nbar_low}))[0] - istdat(std::vector<double>({nbar_upp}))[0]) >= density_prec)
                    while (fabs(nbar_upp - nbar_low) > nbar_min)
                    {
                        // while we are too far from appropriate precision for nbar estimate
                        // recalculate via bisection method
                        nbar_mid = (nbar_low + nbar_upp) / 2.0;
                        if (energy_density_conversion * data_reader(std::vector<double>({nbar_mid}))[output_index_energy_density] > density_at_r)
                            nbar_upp = nbar_mid;
                        else
                            nbar_low = nbar_mid;
                    }
                    cache[1].push_back(nbar_mid);
                }
                cache[0].resize(cache[1].size()); // truncate radii array so that to fit to nbar
            }
            return std::vector<double>({nbar_interpolator(cache[0], cache[1], r), cache[0].back()});
        });

    // ready to calculate

    double r_ns = tov(0.0)[4];
    double r_crust = nbar(0.0)[1];
    double m_ns = tov(r_ns)[0];
    std::cout << m_ns * constants::conversion::gev_over_msol << " ";
    for (size_t species = 0; species < Y_i_functions_of_nbar.size(); ++species)
    {
        double I_i = 0.0;
        auto Y_i = Y_i_functions_of_nbar[species];
        double omega_k_sqr = pow(2.0 / 3, 3.0) * constants::scientific::G * m_ns / (r_ns * r_ns * r_ns);
        auto integrand = [&](double r)
        {
            return -4.0 * M_PI * r * r * nbar_conversion * nbar(r)[0] * 1.0 / omega_k_sqr * (Y_i(nbar(r + radius_step)[0]) - Y_i(nbar(r)[0])) / (tov(r + radius_step)[3] / tov(r)[3] - 1);
        };
        for (double r = 0; r < r_crust - radius_step; r += radius_step)
            I_i += integrand(r);
        std::cout << I_i * radius_step / (constants::conversion::gev_s * constants::conversion::gev_s) << " ";
    }
}