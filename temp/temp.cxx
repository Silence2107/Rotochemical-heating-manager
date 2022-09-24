
#include "../include/auxiliaries.h"
#include "../include/eos_reader.h"
#include "../include/constants.h"

#include <TGraph.h>
#include <TCanvas.h>
#include <TAxis.h>
#include <TMultiGraph.h>
#include <TLegend.h>

#include <vector>
#include <fstream>
#include <iostream>
#include <cmath>
#include <chrono>


int main()
{
    std::ifstream fstr("../data/IST_NS.TXT");
    auto ist_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>, const std::vector<double> &, std::ifstream &>(eos_reader::predefined::ist_for_ns_cached);

    auto istdat1 = [&fstr, &ist_cached](const std::vector<double> &input)
    {
        return eos_reader::eos_data(input, ist_cached, fstr);
    };

    // a function that calculates the pressure for a given energy density
    auto eos1 = auxiliaries::CachedFunc<std::vector<std::vector<double>>, double, double, int>(
        [&istdat1](std::vector<std::vector<double>> &cache, double rho, int nbar_discretization)
        {
            // APR4 EoS with caching support
            // Recaching happens if cache is empty or in case nbar_discretization changes

            using namespace constants;
            using namespace constants::ist_ns;
            if (rho < 0 || rho > conversion::mev_over_fm3_gev4 * edensity_upp)
                throw std::runtime_error("Data request out of range; Encountered in main::eos_cached");
            if (rho <= conversion::mev_over_fm3_gev4 * edensity_low)
                return 0.0;
            if (cache.empty() || cache[0].size() != nbar_discretization)
            {                                                                                             // then fill/refill cache
                cache = std::vector<std::vector<double>>(2, std::vector<double>(nbar_discretization, 0)); // initialize 2xnbar_discretization matrix
                std::vector<double> x(nbar_discretization, 0);
                for (int i = 1; i < nbar_discretization - 1; ++i)
                { // cache EoS for further efficiency
                    x[i] = i * (nbar_upp - nbar_low) / nbar_discretization + nbar_low;
                    cache[0][i] = conversion::mev_over_fm3_gev4 * istdat1(std::vector<double>({x[i]}))[1];
                    cache[1][i] = conversion::mev_over_fm3_gev4 * istdat1(std::vector<double>({x[i]}))[0];
                }
                x[0] = nbar_low;
                x[x.size() - 1] = nbar_upp;
                cache[0][0] = conversion::mev_over_fm3_gev4 * edensity_low;
                cache[0][cache[0].size() - 1] = conversion::mev_over_fm3_gev4 * edensity_upp;
                cache[1][0] = conversion::mev_over_fm3_gev4 * pressure_low;
                cache[1][cache[1].size() - 1] = conversion::mev_over_fm3_gev4 * pressure_upp;
            }
            return auxiliaries::interpolate(cache[0], cache[1], auxiliaries::InterpolationMode::kLinear, rho);
        });
    
    auto istdat2 = [&fstr](const std::vector<double> &input)
    {
        return eos_reader::eos_data(input, eos_reader::predefined::ist_for_ns, fstr);
    };

    // a function that calculates the pressure for a given energy density
    auto eos2 = auxiliaries::CachedFunc<std::vector<std::vector<double>>, double, double, int>(
        [&istdat2](std::vector<std::vector<double>> &cache, double rho, int nbar_discretization)
        {
            // APR4 EoS with caching support
            // Recaching happens if cache is empty or in case nbar_discretization changes

            using namespace constants;
            using namespace constants::ist_ns;
            if (rho < 0 || rho > conversion::mev_over_fm3_gev4 * edensity_upp)
                throw std::runtime_error("Data request out of range; Encountered in main::eos_cached");
            if (rho <= conversion::mev_over_fm3_gev4 * edensity_low)
                return 0.0;
            if (cache.empty() || cache[0].size() != nbar_discretization)
            {                                                                                             // then fill/refill cache
                cache = std::vector<std::vector<double>>(2, std::vector<double>(nbar_discretization, 0)); // initialize 2xnbar_discretization matrix
                std::vector<double> x(nbar_discretization, 0);
                for (int i = 1; i < nbar_discretization - 1; ++i)
                { // cache EoS for further efficiency
                    x[i] = i * (nbar_upp - nbar_low) / nbar_discretization + nbar_low;
                    cache[0][i] = conversion::mev_over_fm3_gev4 * istdat2(std::vector<double>({x[i]}))[1];
                    cache[1][i] = conversion::mev_over_fm3_gev4 * istdat2(std::vector<double>({x[i]}))[0];
                }
                x[0] = nbar_low;
                x[x.size() - 1] = nbar_upp;
                cache[0][0] = conversion::mev_over_fm3_gev4 * edensity_low;
                cache[0][cache[0].size() - 1] = conversion::mev_over_fm3_gev4 * edensity_upp;
                cache[1][0] = conversion::mev_over_fm3_gev4 * pressure_low;
                cache[1][cache[1].size() - 1] = conversion::mev_over_fm3_gev4 * pressure_upp;
            }
            return auxiliaries::interpolate(cache[0], cache[1], auxiliaries::InterpolationMode::kLinear, rho);
        });

    //compare timings of both EoSs with std::chrono
    auto start = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 40000; ++i)
    {
        using namespace constants::ist_ns;
        eos1(constants::conversion::mev_over_fm3_gev4 * (edensity_low + i * (edensity_upp - edensity_low) / 40000), 2000);
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = end - start;
    std::cout << "Time for eos1: " << elapsed.count() << '\n';
    auto start2 = std::chrono::high_resolution_clock::now();
    for (int i = 0; i < 40000; ++i)
    {
        using namespace constants::ist_ns;
        eos2(constants::conversion::mev_over_fm3_gev4 * (edensity_low + i * (edensity_upp - edensity_low) / 40000), 2000);
    }
    auto end2 = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed2 = end2 - start2;
    std::cout << "Time for eos2: " << elapsed2.count() << '\n';


    // plot both EoSs side by side on TGraph
    TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
    c1->SetLogy();
    c1->SetLogx();
    TGraph *g1 = new TGraph();
    TGraph *g2 = new TGraph();
    g1->SetLineColor(kRed);
    g2->SetLineColor(kBlue);
    g1->SetLineWidth(2);
    g2->SetLineWidth(2);
    g1->SetMarkerColor(kRed);
    g2->SetMarkerColor(kBlue);
    g1->SetMarkerStyle(20);
    g2->SetMarkerStyle(20);
    g1->SetMarkerSize(0.5);
    g2->SetMarkerSize(0.5);
    TMultiGraph *mg = new TMultiGraph();
    mg->Add(g1, "2CP");
    mg->Add(g2, "2CP");
    mg->GetXaxis()->SetTitle("Energy density [GeV^4]");
    mg->GetYaxis()->SetTitle("Pressure [GeV^4]");
    // Set ranges for x and y axes with lowest and highest values
    mg->GetXaxis()->SetLimits(constants::conversion::mev_over_fm3_gev4*constants::ist_ns::edensity_low, constants::conversion::mev_over_fm3_gev4*constants::ist_ns::edensity_upp);
    mg->GetYaxis()->SetRangeUser(constants::conversion::mev_over_fm3_gev4*constants::ist_ns::pressure_low, constants::conversion::mev_over_fm3_gev4*constants::ist_ns::pressure_upp);
    for (int i = 1; i < 4000; ++i)
    {
        using namespace constants::ist_ns;
        // log scale for x axis
        double x = constants::conversion::mev_over_fm3_gev4 * pow(10, log10(edensity_low) + (log10(edensity_upp) - log10(edensity_low)) * i / 4000);
        g1->SetPoint(i, x, eos1(x, 4000));
        g2->SetPoint(i, x, eos2(x, 4000));
    }
    mg->Draw("APL");
    TLegend *leg = new TLegend(0.1, 0.7, 0.48, 0.9);
    leg->AddEntry(g1, "EOS1", "l");
    leg->AddEntry(g2, "EOS2", "l");
    leg->Draw();
    c1->SaveAs("eos_comparison.pdf");

    return 0;
}