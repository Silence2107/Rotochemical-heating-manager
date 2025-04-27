
#include <iostream>
#include <fstream>
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

    auto eos_inv_cached = edensity_of_pressure;

    // TOV solver

    auto tov_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                                    const std::function<double(double)> &, double, double, double, double, size_t>(tov_solver::tov_solution);
    auto tov = [&tov_cached, &eos_inv_cached, center_pressure](double r)
    {
        // TOV solution cached
        return tov_cached(eos_inv_cached, r, center_pressure, radius_step, surface_pressure, tov_adapt_limit);
    };

    auto nbar = [&](double r)
    {
        return nbar_of_pressure(tov(r)[3]);
    };

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