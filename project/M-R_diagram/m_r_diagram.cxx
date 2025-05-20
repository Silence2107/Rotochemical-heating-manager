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
    std::string program_name = "m_r_diagram";
    argparse::ArgumentParser parser(program_name, "Maps central pressures selection to corresponding mass and radius based on EoS", "Argparse powered by SiLeader");

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

    auxiliaries::io::Logger logger(program_name);

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

    logger.log([]()
               { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
               [&]()
               { return "Instantiation complete"; });

    // RUN --------------------------------------------------------------------------

    // EoS definition

    auto eos_inv_cached = edensity_of_pressure;

    // returns {r, m, n_c} triple at given center pressure. (Hopefully) cleans up all global cache that may spoil further calls
    auto get_m_r_at_pressure = [&](double pressure)
    {
        // TOV solver
        auto tov_df = tov_solver::tov_solution(eos_inv_cached, pressure, radius_step, surface_pressure, pressure_low, tov_adapt_limit);
        
        double r_ns = tov_df[0].back();
        double m_ns = tov_df[1].back();
        double central_nbar = nbar_of_pressure(pressure);
        // think twice here if you need to clean up any global cache
        // rho(P), nb(P) is the same between runs, so cleaning them is unnecessary work

        return std::vector<double>({r_ns, m_ns, central_nbar});
    };

    size_t indent = 20;
    std::vector<double> radii, masses, densities;
    // assemble data for different center pressures
    logger.log([]()
               { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
               [&]()
               {
                   std::stringstream ss;
                   ss << std::scientific << std::setprecision(3) << "Pressure fraction is exp mapped [" << left_fraction << ", " << right_fraction << ", " << selection_size << "]";
                   return ss.str();
               });
    std::cout << std::left << std::setw(indent) << "pressure_fraction" << std::setw(indent) << "density_c[fm-3]" << std::setw(indent) << "M[Ms]" << std::setw(indent) << "R[km]" << '\n';
    for (size_t count = 0; count < selection_size; ++count)
    {
        using namespace constants::conversion;
        double frac = left_fraction * pow((right_fraction / left_fraction), count / (selection_size - 1.0));
        double pressure = frac * (pressure_upp - pressure_low) + pressure_low;
        auto point = get_m_r_at_pressure(pressure);
        if (restrict_stable_branch && count != 0)
            // Early exit if dM/dP < 0 with M/Ms > 1.8
            if (masses.back() > point[1] * gev_over_msol && point[1] * gev_over_msol > 1.8)
            {
                logger.log([]()
                           { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
                           [&]()
                           { return std::string("Reached potential unstable configuration at ") + std::to_string(point[1] * gev_over_msol) + " Ms"; });
                break;
            }
        radii.push_back(point[0] / km_gev);
        masses.push_back(point[1] * gev_over_msol);
        densities.push_back(point[2] * constants::conversion::fm3_gev3);
        logger.log([&]()
                   { return count % 100 == 0; }, auxiliaries::io::Logger::LogLevel::kInfo,
                   [&]()
                   { return std::to_string(count) + " counts past"; }, "M-R loop");
        logger.log([&]()
                   { return radii.back() > 50.0; }, auxiliaries::io::Logger::LogLevel::kDebug,
                   [&]()
                   { return "NS radius exceeds 50 km. Consider entering stiffer area (perhaps raise --left_fraction)"; }, "M-R loop");
        logger.log([&]()
                   { return masses.back() > 3.0; }, auxiliaries::io::Logger::LogLevel::kDebug,
                   [&]()
                   { return "NS mass exceeds 3 Ms. Consider halting the calculation (perhaps pass --restrict_stable or reduce --right_fraction)"; }, "M-R loop");
        logger.log([&]()
                   { return true; }, auxiliaries::io::Logger::LogLevel::kTrace,
                   [&]()
                   { 
                        std::stringstream ss;
                        ss << std::to_string(count) + " counts past. ";
                        ss << "R[km] = " << radii.back() << ", M[Ms] = " << masses.back();
                        return ss.str(); }, "M-R loop");
        std::cout << std::left << std::setw(indent) << frac << std::setw(indent) << densities.back() << std::setw(indent) << masses.back() << std::setw(indent) << radii.back() << "\n";
    }

#if RHM_HAS_ROOT
    // draw
    TCanvas *c1 = new TCanvas("c1", "c1");
    auto gr = new TGraph(radii.size(), radii.data(), masses.data());

    gr->GetXaxis()->SetTitle("R [km]");
    gr->GetYaxis()->SetTitle("M [Ms]");
    if (rootfile)
    {
        auto gr_n = new TGraph(densities.size(), densities.data(), radii.data());
        gr_n->GetXaxis()->SetTitle("nb [fm-3]");
        gr_n->GetYaxis()->SetTitle("R [km]");
        auto gr_m = new TGraph(densities.size(), densities.data(), masses.data());
        gr_m->GetXaxis()->SetTitle("nb [fm-3]");
        gr_m->GetYaxis()->SetTitle("M [Ms]");
        rootfile->cd();
        rootfile->WriteObject(gr, "m_r_diagram");
        rootfile->WriteObject(gr_n, "r_n_diagram");
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
