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
    argparse::ArgumentParser parser("m_r_diagram_density_print", "Maps central pressures selection to corresponding mass and radius based on EoS", "Argparse powered by SiLeader");

    parser.addArgument({"--inputfile"}, "json input file path (required)");
#if RHM_HAS_ROOT
    parser.addArgument({"--pdf_path"}, "pdf output file path (optional, default: M-R-diagram.pdf)");
    parser.addArgument({"--rootfile_path"}, "root output file path (optional, default: None)");
#endif
    parser.addArgument({"--left_fraction"}, "least fraction of central pressure to consider (optional, default: 0.0001)");
    parser.addArgument({"--right_fraction"}, "greatest fraction of central pressure to consider (optional, default: 0.999)");
    parser.addArgument({"--selection_size"}, "number of points to discretize pressure interval (optional, default: 5000)");
    parser.addArgument({"--restrict_stable"}, "whether to exit early upon reaching dM/dP < 0 with M > 1.8Ms (optional, value-free, default: false)", argparse::ArgumentType::StoreTrue);
    auto args = parser.parseArgs(argc, argv);

    using namespace instantiator;
    instantiator::instantiate_system(args.get<std::string>("inputfile"), {});

#if RHM_HAS_ROOT
    std::string pdf_path = args.safeGet<std::string>("pdf_path", "M-R-diagram.pdf");
    TFile *rootfile = nullptr;
    if (args.has("rootfile_path"))
        rootfile = new TFile(args.get<std::string>("rootfile_path").c_str(), "RECREATE");
#endif

    double left_fraction = args.safeGet<double>("left_fraction", 0.0001);
    double right_fraction = args.safeGet<double>("right_fraction", 0.999);
    size_t selection_size = args.safeGet<size_t>("selection_size", 5000);
    bool restrict_stable_branch = args.has("restrict_stable");

    // RUN --------------------------------------------------------------------------

    // EoS definition
    
    auto eos_inv_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, double, double>(
        [&](std::vector<std::vector<double>> &cache, double p)
        {
            if (p < pressure_low || p > pressure_upp)
                RHM_ERROR("Data request out of range.");
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

    // returns {r, m} pair at given center pressure. (Hopefully) cleans up all global cache that may spoil further calls
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

    size_t indent = 20;
    std::vector<double> x, y, z;
    // assemble data for different center pressures
    std::cout << std::left << std::setw(indent) << "pressure fraction" << std::setw(indent) << "nbar [fm-3]" << std::setw(indent) << "M [Ms]" << std::setw(indent) << "R [km]" << '\n';
    for (size_t count = 0; count < selection_size; ++count)
    {
        using namespace constants::conversion;
        double frac = left_fraction * pow((right_fraction / left_fraction), count / (selection_size - 1.0));
        double pressure = frac * (pressure_upp - pressure_low) + pressure_low;
        // locate baryonic density with bisection
        double nbar_left = nbar_low, nbar_right = nbar_upp, nbar_mid = 0;
        while (nbar_right - nbar_left > nbar_low)
        {
            nbar_mid = (nbar_left + nbar_right) / 2;
            if (pressure_of_nbar(nbar_mid) < pressure)
                nbar_left = nbar_mid;
            else
                nbar_right = nbar_mid;
        }
        auto point = get_m_r_at_pressure(pressure);
        if(restrict_stable_branch && count != 0)
            // Early exit if dM/dP < 0 with M/Ms > 1.8
            if (y.back() > point[1] * gev_over_msol && point[1] * gev_over_msol > 1.8)
                break;
        x.push_back(point[0] / km_gev);
        y.push_back(point[1] * gev_over_msol);
        z.push_back(nbar_mid * constants::conversion::fm3_gev3);
        std::cout << std::left << std::setw(indent) << frac << std::setw(indent) << z.back() << std::setw(indent) << y.back() << std::setw(indent) << x.back() << "\n";
    }

#if RHM_HAS_ROOT
    // draw
    TCanvas *c1 = new TCanvas("c1", "c1");
    auto gr = new TGraph(x.size(), x.data(), y.data());

    gr->GetXaxis()->SetTitle("R [km]");
    gr->GetYaxis()->SetTitle("M [Ms]");
    if (rootfile)
    {
        auto gr_m = new TGraph(z.size(), z.data(), y.data());
        gr_m->GetXaxis()->SetTitle("n [fm^{-3}]");
        gr_m->GetYaxis()->SetTitle("M [Ms]");
        rootfile->cd();
        rootfile->WriteObject(gr, "m_r_diagram");
        rootfile->WriteObject(gr_m, "m_n_diagram");
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
