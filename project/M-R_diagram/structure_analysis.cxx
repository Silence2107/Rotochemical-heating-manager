
#include "../../include/auxiliaries.h"
#include "../../include/cooling.h"
#include "../../include/constants.h"
#include "../../include/tov_solver.h"
#include "../../include/instantiator.hpp"

#include "../../3rd-party/argparse/argparse.hpp"

#include <vector>
#include <functional>
#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

#if RHM_HAS_ROOT
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>
#endif

int main(int argc, char **argv)
{
    std::string program_name = "structure_analysis";
    argparse::ArgumentParser parser(program_name, "Extract structural information based on EoS. By default, prints the current mass and radius of the star.", "Argparse powered by SiLeader");

    parser.addArgument({"--inputfile"}, "json input file path (required)");
    parser.addArgument({"--max_mass"}, "whether to search for max mass point on M-R curve (optional, value-free, default: disabled)", argparse::ArgumentType::StoreTrue);
    parser.addArgument({"--deconfinement"}, "whether to search for deconfinement emergence and quark onset mass (optional, value-free, default: disabled)", argparse::ArgumentType::StoreTrue);
    auto args = parser.parseArgs(argc, argv);
    bool search_max_mass = args.has("max_mass");
    bool search_deconfinement = args.has("deconfinement");

    using namespace instantiator;
    if(search_deconfinement)
        instantiator::instantiate_system(args.get<std::string>("inputfile"), {"TOV", "COOL"});
    else
        instantiator::instantiate_system(args.get<std::string>("inputfile"), {"TOV"});

    auxiliaries::io::Logger logger(program_name);

    logger.log([]()
               { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
               [&]()
               { return "Instantiation complete"; });
    // RUN --------------------------------------------------------------------------

    // EoS definition

    auto eos_inv_cached = edensity_of_pressure;

    // TOV solver

    auto tov_df = tov_solver::tov_solution(eos_inv_cached, center_pressure, radius_step, surface_pressure, pressure_low, tov_adapt_limit);

    std::vector<double> df_nbar(tov_df[0].size()), df_exp_phi(tov_df[0].size()), df_exp_lambda(tov_df[0].size());

    double r_ns = tov_df[0].back();
    double m_ns = tov_df[1].back();
    for (size_t i = 0; i < tov_df[0].size(); ++i)
    {
        df_nbar[i] = nbar_of_pressure(tov_df[2][i]);
        df_exp_phi[i] = std::exp(tov_df[3][i]);
        if (i == 0)
            df_exp_lambda[i] = 1.0;
        else
            df_exp_lambda[i] = std::pow(1 - 2 * constants::scientific::G * tov_df[1][i] / tov_df[0][i], -0.5);
    }
    auto nbar = auxiliaries::math::Interpolator(tov_df[0], df_nbar, radial_interp_mode);
    // auto exp_phi = auxiliaries::math::Interpolator(tov_df[0], df_exp_phi, radial_interp_mode);
    // auto exp_lambda = auxiliaries::math::Interpolator(tov_df[0], df_exp_lambda, radial_interp_mode);

    size_t indent = 20;
    std::cout << std::left << std::setw(indent) << "Mass[Ms]" << std::setw(indent) << "Radius[km]";
    if(search_max_mass)
        std::cout << std::setw(indent) << "Max_mass[Ms]" << std::setw(indent) << "Max_mass_radius[km]";
    if(search_deconfinement)
        std::cout << std::setw(indent) << "Quark_onset[Ms]" << std::setw(indent) << "Deconf_emergence[km]";
    std::cout << '\n';
    // print current mass and radius
    std::cout << std::left << std::setw(indent) << m_ns * constants::conversion::gev_over_msol;
    std::cout << std::setw(indent) << r_ns / constants::conversion::km_gev;

    if (search_max_mass)
    {
        // produce total M-R curve
        auto get_m_r_at_pressure = [&](double pressure)
        {
            // TOV solver
            auto tov_df = tov_solver::tov_solution(eos_inv_cached, pressure, radius_step, surface_pressure, pressure_low, tov_adapt_limit);

            double r_ns = tov_df[0].back();
            double m_ns = tov_df[1].back();
            // think twice here if you need to clean up any global cache
            // We memorized P(rho), but cleaning it is doing extra unnecessary work

            return std::vector<double>({r_ns, m_ns});
        };

        // size_t indent = 20;
        size_t selection_size = 1000;
        double left_fraction = 0.001,
               right_fraction = 0.999;
        std::vector<double> x, y, z;
        // assemble data for different center pressures
        for (size_t count = 0; count < selection_size; ++count)
        {
            using namespace constants::conversion;
            double frac = left_fraction * pow((right_fraction / left_fraction), count / (selection_size - 1.0));
            double pressure = frac * (pressure_upp - pressure_low) + pressure_low;
            auto point = get_m_r_at_pressure(pressure);
            x.push_back(point[0] / km_gev);
            y.push_back(point[1] * gev_over_msol);
            z.push_back(pressure / pressure_conversion);
        }
        // retrieve index of max mass
        auto max_mass_it = std::max_element(y.begin(), y.end());
        // print max mass, corresponding radius
        std::cout << std::setw(indent) << *max_mass_it;
        std::cout << std::setw(indent) << x[std::distance(y.begin(), max_mass_it)];
        logger.log([]()
                   { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
                   [&]()
                   { return "Maximum mass resolved"; });
    }
    if (search_deconfinement)
    {
        // calculate the mass of quark appearance
        if (number_densities_of_nbar.find(constants::species::uquark) != number_densities_of_nbar.end())
        {
            auto udens = number_densities_of_nbar[constants::species::uquark];
            double nb_left = nbar_low,
                   nb_right = nbar_upp;
            double q_onset_density = (nb_left + nb_right) / 2;
            if (udens(nbar_upp) == 0.0)
                q_onset_density = nbar_upp;
            else if (udens(nbar_low) != 0.0)
                q_onset_density = nbar_low;
            else
                while (nb_right - nb_left > nbar_low)
                {
                    q_onset_density = (nb_left + nb_right) / 2;
                    if (udens(q_onset_density) == 0.0)
                    {
                        nb_left = q_onset_density;
                        continue;
                    }
                    else
                    {
                        nb_right = q_onset_density;
                        continue;
                    }
                }
            double q_onset_pressure = pressure_of_nbar(q_onset_density);

            auto tov_df_at_transition = tov_solver::tov_solution(eos_inv_cached, q_onset_pressure, radius_step, surface_pressure, pressure_low, tov_adapt_limit);
            // double r_ns_at_transition = tov_df_at_transition[0].back();
            double m_ns_at_transition = tov_df_at_transition[1].back();
            // print quark onset mass
            std::cout << std::setw(indent) << m_ns_at_transition * constants::conversion::gev_over_msol;

            // deconfinement emergence
            double r_left = 0.0, r_right = r_ns, r_deconf = r_ns / 2;
            auto deconf_equation = [&](double r)
            {
                return nbar(r) - q_onset_density;
            };
            while (r_right - r_left > radius_step / 2)
            {
                r_deconf = (r_left + r_right) / 2;
                if (deconf_equation(r_left) * deconf_equation(r_deconf) <= 0)
                    r_right = r_deconf;
                else if (deconf_equation(r_right) * deconf_equation(r_deconf) <= 0)
                    r_left = r_deconf;
                else
                {
                    r_deconf = NAN;
                    break;
                }
            }
            // print deconfinement radius
            std::cout << std::setw(indent) << r_deconf/constants::conversion::km_gev;
        }
        else
        {
            std::cout << std::setw(indent) << NAN;
            std::cout << std::setw(indent) << NAN;
        }
        logger.log([]()
                   { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
                   [&]()
                   { return "Deconfinement properties resolved"; });
    }
}