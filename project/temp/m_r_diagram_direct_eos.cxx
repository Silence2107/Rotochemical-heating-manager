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
#include <string>
#include <fstream>

#if RHM_HAS_ROOT
#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TFile.h>
#endif

int main(int argc, char **argv)
{
    argparse::ArgumentParser parser("m_r_diagram", "Maps central densities selection to corresponding mass and radius based on EoS", "Argparse powered by SiLeader");

#if RHM_REQUIRES_INPUTFILE
    parser.addArgument({"--inputfile"}, "json input file path (required)");
#endif
#if RHM_HAS_ROOT
    parser.addArgument({"--pdf_path"}, "pdf output file path (optional, default: M-R-diagram.pdf)");
    parser.addArgument({"--rootfile_path"}, "root output file path (optional, default: None)");
#endif
    parser.addArgument({"--left_fraction"}, "least fraction of central density to consider (optional, default: 0.001)");
    parser.addArgument({"--right_fraction"}, "greatest fraction of central density to consider (optional, default: 0.999)");
    parser.addArgument({"--selection_size"}, "number of points to discretize density interval (optional, default: 1000)");
    auto args = parser.parseArgs(argc, argv);

    using namespace instantiator;
#if RHM_REQUIRES_INPUTFILE
    instantiator::instantiate_system(args.get<std::string>("inputfile"));
#endif

#if RHM_HAS_ROOT
    std::string pdf_path = args.safeGet<std::string>("pdf_path", "M-R-diagram.pdf");
    TFile *rootfile = nullptr;
    if (args.has("rootfile_path"))
        rootfile = new TFile(args.get<std::string>("rootfile_path").c_str(), "RECREATE");
#endif

    double left_fraction = args.safeGet<double>("left_fraction", 0.001);
    double right_fraction = args.safeGet<double>("right_fraction", 0.999);
    size_t selection_size = args.safeGet<size_t>("selection_size", 1000);

    // RUN --------------------------------------------------------------------------

    // EoS definition

    auto eos_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, double, double>(
        [&](std::vector<std::vector<double>> &cache, double rho)
        {
            if (rho < edensity_low || rho > edensity_upp)
                RHM_THROW(std::runtime_error, "Data request out of range.");
            if (cache.empty() || cache[0].size() != discr_size_EoS)
            {                                                                                        // then fill/refill cache
                cache = std::vector<std::vector<double>>(2, std::vector<double>(discr_size_EoS, 0)); // initialize 2xdiscr_size_EoS matrix
                std::vector<double> x(discr_size_EoS, 0);
                for (size_t i = 1; i < discr_size_EoS - 1; ++i)
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

    // returns {r, m} pair at given center density. (Hopefully) cleans up all global cache that may spoil further calls
    auto get_m_r_at_density = [&](double edensity)
    {
        // TOV solver
        auto tov_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                                        const std::function<double(double)> &, double, double, double, double, size_t>(tov_solver::tov_solution_direct_eos);
        auto tov = [&tov_cached, &eos_cached, edensity](double r)
        {
            // TOV solution cached
            return tov_cached(eos_cached, r, edensity, radius_step, surface_pressure, tov_adapt_limit);
        };

        double r_ns = tov(0.0)[4];
        double m_ns = tov(r_ns)[0];
        // think twice here if you need to clean up any global cache
        // We memorized P(rho), but cleaning it is doing extra unnecessary work

        return std::vector<double>({r_ns, m_ns});
    };

    size_t indent = 20;
    std::vector<double> x, y, z;
    // assemble data for different center densities
    std::cout << std::left << std::setw(indent) << "rho fraction" << std::setw(indent) << "rho [df. units]" << std::setw(indent) << "M [Ms]" << std::setw(indent) << "R [km]" << '\n';
    for (size_t count = 0; count < selection_size; ++count)
    {
        using namespace constants::conversion;
        double frac = left_fraction + count * (right_fraction - left_fraction) / (selection_size - 1);
        double edensity = frac * (edensity_upp - edensity_low) + edensity_low;
        auto point = get_m_r_at_density(edensity);
        x.push_back(point[0] / km_gev);
        y.push_back(point[1] * gev_over_msol);
        z.push_back(edensity / energy_density_conversion);
        std::cout << std::left << std::setw(indent) << frac << std::setw(indent) << edensity / energy_density_conversion << std::setw(indent) << y.back() << std::setw(indent) << x.back() << "\n";
    }

#if RHM_HAS_ROOT
    // draw
    TCanvas *c1 = new TCanvas("c1", "c1");
    auto gr = new TGraph(x.size(), x.data(), y.data());

    gr->GetXaxis()->SetTitle("R [km]");
    gr->GetYaxis()->SetTitle("M [Ms]");
    if (rootfile)
    {
        auto gr_rho = new TGraph(z.size(), z.data(), x.data());
        gr_rho->GetXaxis()->SetTitle("#rho [datafile units]");
        gr_rho->GetYaxis()->SetTitle("R [km]");
        auto gr_m = new TGraph(z.size(), z.data(), y.data());
        gr_m->GetXaxis()->SetTitle("#rho [datafile units]");
        gr_m->GetYaxis()->SetTitle("M [Ms]");
        rootfile->cd();
        rootfile->WriteObject(gr, "m_r_diagram");
        rootfile->WriteObject(gr_rho, "rho_r_diagram");
        rootfile->WriteObject(gr_m, "m_rho_diagram");
        rootfile->Close();
    }
    gr->SetLineColor(kBlue);
    gr->Draw("AL");
    // title offset
    gr->GetYaxis()->SetTitleOffset(1.5);
    gr->GetXaxis()->SetLimits(0, 20.0);
    gr->GetYaxis()->SetRangeUser(0, 3.0);
    // gPad->SetLogx();
    // gPad->SetLogy();

    // auto legend = new TLegend(0.1, 0.1, 0.38, 0.38);
    // legend->AddEntry(gr, "RH Manager", "l");

    // legend->Draw();

    c1->SaveAs(pdf_path.c_str());
#endif
}