
#include "../../include/auxiliaries.h"
#include "../../include/cooling.h"
#include "../../include/constants.h"
#include "../../include/tov_solver.h"
#include "../../include/inputfile.hpp"

#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <iostream>

#include <sstream>
#include <fstream>

#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TFile.h>
#include <TStyle.h>

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        std::cout << "Usage: " << argv[0] << " <pdf_path=znpe.pdf> <rootfile_path=None>" << std::endl;
    }
    std::string pdf_path = (argc > 1) ? argv[1] : "znpe.pdf";
    bool rootfile_creation = (argc > 2);
    using namespace inputfile;

    // RUN --------------------------------------------------------------------------

    auto eos_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, double, double>(
        [&](std::vector<std::vector<double>> &cache, double rho)
        {
            if (rho < 0 || rho > edensity_upp)
                throw std::runtime_error("Data request out of range; Encountered in main::eos_cached");
            if (rho <= edensity_low)
                return 0.0;
            if (cache.empty() || cache[0].size() != discr_size_EoS)
            {                                                                                        // then fill/refill cache
                cache = std::vector<std::vector<double>>(2, std::vector<double>(discr_size_EoS, 0)); // initialize 2xdiscr_size_EoS matrix
                std::vector<double> x(discr_size_EoS, 0);
                for (int i = 1; i < discr_size_EoS - 1; ++i)
                { // cache EoS for further efficiency
                    x[i] = i * (nbar_upp - nbar_low) / discr_size_EoS + nbar_low;
                    cache[0][i] = energy_density_of_nbar(x[i]);
                    cache[1][i] = pressure_of_nbar(x[i]);
                }
                x[0] = nbar_low;
                x[x.size() - 1] = nbar_upp;
                cache[0][0] = edensity_low;
                cache[0][cache[0].size() - 1] = edensity_upp;
                cache[1][0] = pressure_low;
                cache[1][cache[1].size() - 1] = pressure_upp;
                eos_interpolator_cached.erase(); // clean up cached interpolator
            }
            return eos_interpolator(cache[0], cache[1], rho);
        });

    // EoS definition
    auto get_znpe = [&](double edensity)
    {
        // TOV solver

        auto tov_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                                const std::function<double(double)> &, double, double, double, double>(tov_solver::tov_solution);
        auto tov = [&tov_cached, &eos_cached, edensity](double r)
        {
            // TOV solution cached
            return tov_cached(eos_cached, r, edensity, radius_step, density_step);
        };

        double r_crust;
        bool crust_found = false;

        auto nbar = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, double, double>(
            [&](std::vector<std::vector<double>> &cache, double r)
            {
                // cache contains {r, n_B(r)} arrays; recaching is not supported at the moment, call ::erase instead
                // return nbar(r) for given r

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
                        if (!crust_found && density_at_r < edensity_core_limit)
                        {
                            crust_found = true;
                            r_crust = r_current;
                        }
                        double nbar_left = nbar_low, nbar_right = nbar_upp; // we need these for bisection search;
                        double nbar_mid = (nbar_left + nbar_right) / 2.0;
                        while (fabs(nbar_right - nbar_left) > nbar_low)
                        {
                            // while we are too far from appropriate precision for nbar estimate
                            // recalculate via bisection method
                            nbar_mid = (nbar_left + nbar_right) / 2.0;
                            if (energy_density_of_nbar(nbar_mid) > density_at_r)
                                nbar_right = nbar_mid;
                            else
                                nbar_left = nbar_mid;
                        }
                        cache[1].push_back(nbar_mid);
                    }
                    cache[0].push_back(R_ns);
                    cache[1].push_back(0.0);
                }
                return nbar_interpolator(cache[0], cache[1], r);
            });

        double r_ns = tov(0.0)[4];
        double m_ns = tov(r_ns)[0];
        nbar(0.0); // invoke in order to find r_crust
        if (!crust_found)
            r_crust = r_ns;

        auto exp_phi = [&tov](double r)
        {
            return std::exp(tov(r)[2]);
        };

        auto exp_lambda = [&tov](double r)
        {
            return pow(1 - 2 * constants::scientific::G * tov(r)[0] / r, -0.5);
        };

        // Znpe rotochemical
        std::function<double(double)> integrand = [&](double r)
        {
            using namespace constants::conversion;
            return data_reader({nbar(r)}, 43) / exp_phi(r) / fm3_gev3 * gev_over_mev;
        };
        double z_npe = 1/auxiliaries::math::integrate_volume<>(integrand, 0.0, r_crust, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)();
        
        // clear
        nbar_interpolator_cached.erase();
        return std::vector<double>({m_ns, z_npe});
    };

    // tabulate
    size_t n = 500;
    size_t offset = 15;
    std::vector<double> x, y;
    // assemble data for different center densities
    for(size_t count = 0.08 * n; count < 0.2 * n; ++count)
    {
        using namespace constants::conversion;
        double edensity = (count * edensity_upp) / n;
        auto point = get_znpe(edensity);
        x.push_back(point[0] * gev_over_msol);
        y.push_back(point[1] / constants::conversion::erg_over_gev);
        std::cout << "At " << count << " out of " << n << " cycles. M = " << x.back() << " Ms, Ze = " << y.back() << " ergs \n";
    }

    // draw
    TCanvas *c1 = new TCanvas("c1", "c1");
    auto gr = new TGraph(x.size(), x.data(), y.data());
    if (rootfile_creation)
    {
        std::string rootfile_path = argv[2];
        TFile *f = new TFile(rootfile_path.c_str(), "RECREATE");
        gr->Write();
        f->Close();
    }
    gr->SetLineColor(kBlue);
    gr->Draw("AC");
    
    gPad->SetTicks();
    gPad->SetTopMargin(0.05);
    gPad->SetLeftMargin(0.11);
    gPad->SetRightMargin(0.05); 
    gPad->SetBottomMargin(0.1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gr->GetXaxis()->SetLimits(1.0, 2.0);
    gr->GetYaxis()->SetRangeUser(1E-61, 6E-60);
    
    gr->GetXaxis()->SetLabelFont(43);
    gr->GetXaxis()->SetLabelSize(22);
    gr->GetXaxis()->SetTitleFont(43);
    gr->GetXaxis()->SetTitleSize(22);
    gr->GetXaxis()->SetTitleOffset(0.9);
    gr->GetYaxis()->SetLabelFont(43);
    gr->GetYaxis()->SetLabelSize(22);
    gr->GetYaxis()->SetTitleFont(43);
    gr->GetYaxis()->SetTitleSize(22);
    gr->GetYaxis()->SetTitleOffset(0.9);
    gr->GetXaxis()->SetTitle("M [M_{#odot}]");
    gr->GetYaxis()->SetTitle("Z_{npe} [erg]");

    c1->SaveAs(pdf_path.c_str());
}