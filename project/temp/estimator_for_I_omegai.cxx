
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

    auto tov_df = tov_solver::tov_solution(eos_inv_cached, center_pressure, radius_step, surface_pressure, pressure_low, tov_adapt_limit);

    std::vector<double> df_nbar(tov_df[0].size()), df_exp_phi(tov_df[0].size()), df_exp_lambda(tov_df[0].size()), df_pressure(tov_df[0].size());

    double r_ns = tov_df[0].back();
    double m_ns = tov_df[1].back();
    for (size_t i = 0; i < tov_df[0].size(); ++i)
    {
        df_pressure[i] = tov_df[2][i];
        df_nbar[i] = nbar_of_pressure(tov_df[2][i]);
        df_exp_phi[i] = std::exp(tov_df[3][i]);
        if (i == 0)
            df_exp_lambda[i] = 1.0;
        else
            df_exp_lambda[i] = std::pow(1 - 2 * constants::scientific::G * tov_df[1][i] / tov_df[0][i], -0.5);
    }
    auto nbar = auxiliaries::math::Interpolator(tov_df[0], df_nbar, radial_interp_mode);
    auto exp_phi = auxiliaries::math::Interpolator(tov_df[0], df_exp_phi, radial_interp_mode);
    auto exp_lambda = auxiliaries::math::Interpolator(tov_df[0], df_exp_lambda, radial_interp_mode);
    auto pressure_of_r = auxiliaries::math::Interpolator(tov_df[0], df_pressure, radial_interp_mode);

    std::cout << "M / M_sol " << m_ns * constants::conversion::gev_over_msol << "\n";

    for (auto species = number_densities_of_nbar.begin(); species != number_densities_of_nbar.end(); ++species)
    {
        auto n_i = species->second;
        double omega_k_sqr = pow(2.0 / 3, 3.0) * constants::scientific::G * m_ns / (r_ns * r_ns * r_ns);
        auto integrand = [&](double r)
        {
            return -nbar(r) * 1.0 / omega_k_sqr * (n_i(nbar(r + radius_step)) / nbar(r + radius_step) - n_i(nbar(r)) / nbar(r)) / (pressure_of_r(r + radius_step) / pressure_of_r(r) - 1);
        };
        double I_i = auxiliaries::math::integrate_volume<>(std::function<double(double)>(integrand), 0.0, r_ns - radius_step, exp_lambda, auxiliaries::math::IntegrationMode::kRectangular)();
        std::cout << species->first.name() << " " << I_i / (constants::conversion::gev_s * constants::conversion::gev_s) << "\n";
    }
}