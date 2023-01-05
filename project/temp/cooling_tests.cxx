
#include "../include/auxiliaries.h"
#include "../include/cooling.h"
#include "../include/eos_reader.h"
#include "../include/constants.h"
#include "../include/tov_solver.h"
#include "../include/inputfile.hpp"

#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>

int main()
{
    using namespace inputfile;

    // Cooling solver setup
    auto cooling_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                               double, const std::vector<double> &, const std::vector<double> &,
                                                               auxiliaries::InterpolationMode, double, bool>(auxiliaries::interpolate_cached);
    auto cooling_interpolator = [&cooling_interpolator_cached](const std::vector<double> &x, const std::vector<double> &y, double val)
    {
        return cooling_interpolator_cached(x, y, auxiliaries::InterpolationMode::kCubic, val, false);
    };

    auto cooling_solver_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>, double, double, const std::function<double(double, double)> &, double, double,
                                                         const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &>(cooling::solver::stationary_cooling_cached);

    // RUN --------------------------------------------------------------------------

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
    auto tov = [&tov_cached, &eos_cached](double r)
    {
        // TOV solution cached
        return tov_cached(eos_cached, r, center_density, radius_step, density_step);
    };

    // ready to run

    double r_ns = tov(0.0)[4];
    double m_ns = tov(r_ns)[0];
    double eta = 1E-18;

    auto photon_luminosity = cooling::predefined::photonic::surface_luminosity(r_ns, m_ns, eta);

    // plot luminocity for temperatures between 1 and 100 MeV
    std::vector<double> x(100, 0);
    std::vector<double> y(100, 0);
    for (int i = 0; i < 100; ++i)
    {
        x[i] = 1 + i * 0.99;
        y[i] = photon_luminosity(0, x[i] / 1000.0) * constants::conversion::gev_s / constants::conversion::erg_over_gev;
    }

    TCanvas *c1 = new TCanvas("c1", "c1");
    auto gr = new TGraph(100, x.data(), y.data());
    gr->Draw("AL");
    // title offset
    gr->GetYaxis()->SetTitleOffset(1.5);

    gr->GetXaxis()->SetTitle("T [MeV]");
    gr->GetYaxis()->SetTitle("L [erg/s]");
    c1->SaveAs("photon_luminosity.pdf");
}