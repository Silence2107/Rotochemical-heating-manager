
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
#include "../include/inputfile.hpp" // inputfile

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

    using namespace inputfile;

    // INPUT:
    // central density
    auto center_density = std::stod(argv[1]) * edensity_upp * energy_density_conversion;

    // EoS definition

    auto eos_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>, double, double>(
        [&](std::vector<std::vector<double>> &cache, double rho)
        {
            if (rho < 0 || rho > edensity_upp * energy_density_conversion)
                throw std::runtime_error("Data request out of range; Encountered in main::eos_cached");
            if (rho <= edensity_low * energy_density_conversion)
                return 0.0;
            if (cache.empty() || cache[0].size() != discr_size_EoS)
            {                                                                                        // then fill/refill cache
                cache = std::vector<std::vector<double>>(2, std::vector<double>(discr_size_EoS, 0)); // initialize 2xdiscr_size_EoS matrix
                std::vector<double> x(discr_size_EoS, 0);
                for (int i = 1; i < discr_size_EoS - 1; ++i)
                { // cache EoS for further efficiency
                    x[i] = i * (nbar_upp - nbar_low) / discr_size_EoS + nbar_low;
                    cache[0][i] = energy_density_conversion * energy_density_of_nbar(x[i]);
                    cache[1][i] = pressure_conversion * pressure_of_nbar(x[i]);
                }
                x[0] = nbar_low;
                x[x.size() - 1] = nbar_upp;
                cache[0][0] = energy_density_conversion * edensity_low;
                cache[0][cache[0].size() - 1] = energy_density_conversion * edensity_upp;
                cache[1][0] = pressure_conversion * pressure_low;
                cache[1][cache[1].size() - 1] = pressure_conversion * pressure_upp;
                eos_interpolator_cached.erase(); // clean up cached interpolator
            }
            return eos_interpolator(cache[0], cache[1], rho);
        });

    // TOV solver

    auto tov_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                              const std::function<double(double)> &, double, double, double, double>(tov_solver::tov_solution);
    auto tov = [&tov_cached, &eos_cached, &center_density](double r)
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
                    if (density_at_r <= energy_density_conversion * edensity_core_limit)
                    { // until we reach core limit
                        // cache[1][i] = 0.0;
                        // continue;
                        break; // finish writing
                    }
                    double nbar_left = nbar_low, nbar_right = nbar_upp; // we need these for bisection search; in fm-3 units for now
                    double nbar_mid = (nbar_left + nbar_right) / 2.0;
                    while (fabs(nbar_right - nbar_left) > nbar_low)
                    {
                        // while we are too far from appropriate precision for nbar estimate
                        // recalculate via bisection method
                        nbar_mid = (nbar_left + nbar_right) / 2.0;
                        if (energy_density_conversion * energy_density_of_nbar(nbar_mid) > density_at_r)
                            nbar_right = nbar_mid;
                        else
                            nbar_left = nbar_mid;
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
    for (auto species = Y_i_functions_of_nbar.begin(); species != Y_i_functions_of_nbar.end(); ++species)
    {
        double I_i = 0.0;
        auto Y_i = species->second;
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