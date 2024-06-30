
#include "../../include/auxiliaries.h"
#include "../../include/cooling.h"
#include "../../include/constants.h"
#include "../../include/tov_solver.h"
#include "../../include/instantiator.hpp"

#include "../../3rd-party/argparse/argparse.hpp"

#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
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
    argparse::ArgumentParser parser("summary", "Extract general information based on EoS", "Argparse powered by SiLeader");

#if RHM_REQUIRES_INPUTFILE
    parser.addArgument({"--inputfile"}, "json input file path (required)");
#endif
    parser.addArgument({"--no_mr_analysis"}, "whether to search for max mass, deconfinement emergence etc (optional, value-free, default: print)", argparse::ArgumentType::StoreTrue);

    auto args = parser.parseArgs(argc, argv);
    bool do_mr_analysis = !args.has("no_mr_analysis");

    using namespace instantiator;
#if RHM_REQUIRES_INPUTFILE
    instantiator::instantiate_system(args.get<std::string>("inputfile"), {"TOV"});
#endif

    // RUN --------------------------------------------------------------------------

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

    auto tov_cached = auxiliaries::math::CachedFunc<std::vector<std::function<double(double)>>, std::vector<double>,
                                                    const std::function<double(double)> &, double, double, double,
                                                    double, size_t, auxiliaries::math::InterpolationMode>(tov_solver::tov_solution);
    auto tov = [&tov_cached, &eos_inv_cached](double r)
    {
        // TOV solution cached
        return tov_cached(eos_inv_cached, r, center_pressure, radius_step, surface_pressure, tov_adapt_limit, radial_interp_mode);
    };

    auto nbar = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, double, double>(
        [&](std::vector<std::vector<double>> &cache, double r) mutable
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

    nbar(0.0); // calculate r_core
    double r_ns = tov(0.0)[4];
    double m_ns = tov(r_ns)[0];
    /*
    auto exp_phi = [&tov](double r)
    {
        return std::exp(tov(r)[2]);
    };

    auto exp_lambda = [&tov](double r)
    {
        return pow(1 - 2 * constants::scientific::G * tov(r)[0] / r, -0.5);
    };*/

    size_t indent = 20;
    std::cout << std::left << std::setw(indent) << "Mass (M sol) : " << m_ns * constants::conversion::gev_over_msol << '\n';
    std::cout << std::left << std::setw(indent) << "Radius (km) : " << r_ns / constants::conversion::km_gev << '\n';

    if (!do_mr_analysis)
        return 0;
    // analyze full M-R
    auto get_m_r_at_pressure = [&](double pressure)
    {
        // TOV solver
        auto tov_cached = auxiliaries::math::CachedFunc<std::vector<std::function<double(double)>>, std::vector<double>,
                                                        const std::function<double(double)> &, double, double, double,
                                                        double, size_t, auxiliaries::math::InterpolationMode>(tov_solver::tov_solution);
        auto tov = [&tov_cached, &eos_inv_cached, pressure](double r)
        {
            // TOV solution cached
            return tov_cached(eos_inv_cached, r, pressure, radius_step, surface_pressure, tov_adapt_limit, radial_interp_mode);
        };

        double r_ns = tov(0.0)[4];
        double m_ns = tov(r_ns)[0];
        // think twice here if you need to clean up any global cache
        // We memorized P(rho), but cleaning it is doing extra unnecessary work

        return std::vector<double>({r_ns, m_ns});
    };

    //size_t indent = 20;
    size_t selection_size = 1000;
    double left_fraction = 0.001,
           right_fraction = 0.999;
    std::vector<double> x, y, z;
    // assemble data for different center pressures
    for (size_t count = 0; count < selection_size; ++count)
    {
        using namespace constants::conversion;
        double frac = left_fraction + count * (right_fraction - left_fraction) / (selection_size - 1);
        double pressure = frac * (pressure_upp - pressure_low) + pressure_low;
        auto point = get_m_r_at_pressure(pressure);
        x.push_back(point[0] / km_gev);
        y.push_back(point[1] * gev_over_msol);
        z.push_back(pressure / pressure_conversion);
    }
    // retrieve index of max mass
    auto max_mass_it = std::max_element(y.begin(), y.end());
    // print max mass, corresponding radius
    std::cout << std::left << std::setw(indent) << "Max mass (M sol) : " << *max_mass_it << '\n';
    std::cout << std::left << std::setw(indent) << "Max mass radius (km) : " << x[std::distance(y.begin(), max_mass_it)] << '\n';

    /*
    // calculate the mass of core appearance
    tov_cached.erase();
    double transition_pressure = pressure_of_nbar(nbar_core_limit);

    auto tov_at_transition = [&tov_cached, &eos_inv_cached, transition_pressure](double r)
    {
        // TOV solution cached
        return tov_cached(eos_inv_cached, r, transition_pressure, radius_step, surface_pressure, tov_adapt_limit, radial_interp_mode);
    };

    // print core emergence radius and mass
    double r_ns_at_transition = tov_at_transition(0.0)[4];
    double m_ns_at_transition = tov_at_transition(r_ns_at_transition)[0];
    std::cout << std::left << std::setw(indent) << "Core emergence mass (M sol) : " << m_ns_at_transition * constants::conversion::gev_over_msol << '\n';
    std::cout << std::left << std::setw(indent) << "Core emergence radius (km) : " << r_ns_at_transition / constants::conversion::km_gev << '\n';*/
}