
#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>
#include <cmath>

#include "../../3rd-party/argparse/argparse.hpp"

#include "../../include/tov_solver.h"     // contains TOV solver
#include "../../include/constants.h"      // contains constants
#include "../../include/auxiliaries.h"    // contains auxiliary functionality
#include "../../include/instantiator.hpp" // instantiator

int main(int argc, char **argv)
{
    // Setup that calculates I_{\omega_i} (based on provided central_density/max_eos_density argv[1]), which is the following quantity:
    // I_{\omega i} = \int_{core} dY_i/dP * dP/d\omega^2 * dN, with estimation dP/d\omega^2 \approx -P/\omega_K^2
    // where \omega_K is the Keplerian frequency and i being particle species
    argparse::ArgumentParser parser("estimator_for_I_omegai", "Estimates I_omega quantities (rotochemical heating related) based on EoS", "Argparse powered by SiLeader");

    parser.addArgument({"--inputfile"}, "json input file path (required)");
    parser.addArgument({"--center_pressure"}, "center pressure linspaced fraction (optional, default: read from inputfile)");

    auto args = parser.parseArgs(argc, argv);

    using namespace instantiator;
    instantiator::instantiate_system(args.get<std::string>("inputfile"), {"TOV", "COOL"});

    double center_pressure = instantiator::center_pressure;
    if (args.has("center_pressure"))
        center_pressure = std::stod(args.get<std::string>("center_pressure")) * (instantiator::pressure_upp - instantiator::pressure_low) + instantiator::pressure_low;

    // EoS definition

    auto eos_inv_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, double, double>(
        [&](std::vector<std::vector<double>> &cache, double p)
        {
            if (p < pressure_low || p > pressure_upp)
                RHM_THROW(std::runtime_error, "Data request out of range.");
            if (cache.empty() || cache[0].size() != discr_size_EoS)
            {                                                                                        // then fill/refill cache
                cache = std::vector<std::vector<double>>(2, std::vector<double>(discr_size_EoS, 0)); // initialize 2xdiscr_size_EoS matrix
                std::vector<double> x(discr_size_EoS, 0);
                for (size_t i = 0; i < discr_size_EoS; ++i)
                { // cache EoS for further efficiency
                    x[i] = nbar_low * pow(nbar_upp / nbar_low, i / (discr_size_EoS - 1.0));
                    cache[0][i] = pressure_of_nbar(x[i]);
                    cache[1][i] = energy_density_of_nbar(x[i]);
                }
                eos_interpolator_cached.erase(); // clean up cached interpolator
            }
            return eos_interpolator(cache[0], cache[1], p);
        });

    // TOV solver

    auto tov_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                                    const std::function<double(double)> &, double, double, double, double, size_t>(tov_solver::tov_solution);
    auto tov = [&tov_cached, &eos_inv_cached, center_pressure](double r)
    {
        // TOV solution cached
        return tov_cached(eos_inv_cached, r, center_pressure, radius_step, surface_pressure, tov_adapt_limit);
    };

    auto nbar = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, double, double>(
        [&](std::vector<std::vector<double>> &cache, double r)
        {
            // cache and evaluate nbar(r) for given r

            if (cache.empty())
            {
                nbar_interpolator_cached.erase(); // clean up cached interpolator
                double R_ns = tov(0.0)[4];
                cache = std::vector<std::vector<double>>(2, std::vector<double>());
                for (double r_current = 0; r_current < R_ns; r_current += radius_step)
                    cache[0].push_back(r_current);
                cache[0].push_back(R_ns);
                for (size_t i = 0; i < cache[0].size(); ++i)
                {
                    double r_current = cache[0][i];
                    // now we somehow have to find corresponding n_B
                    // let's stick to densities
                    double density_at_r = tov(r_current)[1];

                    double nbar_left = nbar_low, nbar_right = nbar_upp; // we need these for bisection search;
                    double nbar_mid = (nbar_left + nbar_right) / 2.0;
                    while (fabs(nbar_right - nbar_left) > nbar_low)
                    {
                        // while we are too far from appropriate precision for nbar estimate
                        // recalculate via bisection method
                        nbar_mid = (nbar_left + nbar_right) / 2.0;
                        double left_val = energy_density_of_nbar(nbar_left) - density_at_r,
                               right_val = energy_density_of_nbar(nbar_right) - density_at_r,
                               mid_val = energy_density_of_nbar(nbar_mid) - density_at_r;
                        if (left_val * mid_val <= 0)
                            nbar_right = nbar_mid;
                        else if (right_val * mid_val <= 0)
                            nbar_left = nbar_mid;
                        else
                            RHM_THROW(std::runtime_error, "Bisection method failed. Investigate manually or report to the team.");
                    }
                    cache[1].push_back(nbar_mid);
                }
            }
            return nbar_interpolator(cache[0], cache[1], r);
        });

    // run

    double r_ns = tov(0.0)[4];
    double m_ns = tov(r_ns)[0];

    std::cout << "M / M_sol " << m_ns * constants::conversion::gev_over_msol << "\n";

    auto exp_lambda = [&tov](double r)
    {
        return pow(1 - 2 * constants::scientific::G * tov(r)[0] / r, -0.5);
    };
    for (auto species = number_densities_of_nbar.begin(); species != number_densities_of_nbar.end(); ++species)
    {
        auto n_i = species->second;
        double omega_k_sqr = pow(2.0 / 3, 3.0) * constants::scientific::G * m_ns / (r_ns * r_ns * r_ns);
        auto integrand = [&](double r)
        {
            return -nbar(r) * 1.0 / omega_k_sqr * (n_i(nbar(r + radius_step)) / nbar(r + radius_step) - n_i(nbar(r)) / nbar(r)) / (tov(r + radius_step)[3] / tov(r)[3] - 1);
        };
        double I_i = auxiliaries::math::integrate_volume<>(std::function<double(double)>(integrand), 0.0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)();
        std::cout << species->first.name() << " " << I_i / (constants::conversion::gev_s * constants::conversion::gev_s) << "\n";
    }
}