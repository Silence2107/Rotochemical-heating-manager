#include <iostream>
#include <fstream>
#include <exception>
#include <vector>
#include <string>
#include <sstream>
#include <iomanip>

#include "../include/eos_reader.h"  // allows to read EoS datafiles
#include "../include/tov_solver.h"  // contains TOV solver
#include "../include/constants.h"   // contains constants
#include "../include/auxiliaries.h" // contains auxiliary functionality

#include <TCanvas.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TMultiGraph.h>
#include <TMath.h>
#include <TStyle.h>
#include <TLegend.h>

int main()
{
    std::ifstream fstr("../data/IST_NS.TXT");
    auto istdat = [&fstr](const std::vector<double> &input)
    {
        return eos_reader::eos_data(input, eos_reader::EOS::IST, fstr);
    };

    auto eos_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>, double, double, int>(
        [&istdat](std::vector<std::vector<double>> &cache, double rho, int nbar_discretization)
        {
            // IST EoS with caching support
            // EoS : energy density rho ([0:~0.12]GeV^4) -> pressure ([0:~1.65]GeV^4)
            // Recaching happens if cache is empty or in case nbar_discretization changes

            using namespace constants;
            if (rho < 0 || rho > conversion::mev_over_fm3_gev4 * ist_ns::edensity_upp)
                throw std::runtime_error("Data request out of range; Encountered in main::eos_cached");
            if (rho <= conversion::mev_over_fm3_gev4 * ist_ns::edensity_low)
                return 0.0;
            if (cache.empty() || cache[0].size() != nbar_discretization)
            {                                                                                             // then fill/refill cache
                cache = std::vector<std::vector<double>>(2, std::vector<double>(nbar_discretization, 0)); // initialize 2xnbar_discretization matrix
                std::vector<double> x(nbar_discretization, 0);
                for (int i = 1; i < nbar_discretization - 1; ++i)
                { // cache EoS for further efficiency
                    x[i] = i * (ist_ns::nbar_upp - ist_ns::nbar_low) / nbar_discretization + ist_ns::nbar_low;
                    cache[0][i] = conversion::mev_over_fm3_gev4 * istdat(std::vector<double>({x[i]}))[1];
                    cache[1][i] = conversion::mev_over_fm3_gev4 * istdat(std::vector<double>({x[i]}))[0];
                }
                x[0] = ist_ns::nbar_low;
                x[x.size() - 1] = ist_ns::nbar_upp;
                cache[0][0] = conversion::mev_over_fm3_gev4 * ist_ns::edensity_low;
                cache[0][cache[0].size() - 1] = conversion::mev_over_fm3_gev4 * ist_ns::edensity_upp;
                cache[1][0] = conversion::mev_over_fm3_gev4 * ist_ns::pressure_low;
                cache[1][cache[1].size() - 1] = conversion::mev_over_fm3_gev4 * ist_ns::pressure_upp;
            }
            return auxiliaries::interpolate(cache[0], cache[1], auxiliaries::InterpolationMode::kLinear, rho);
        });

    auto eos = [&eos_cached](double rho, int nbar_discretization = 2000) -> double
    {
        return eos_cached(rho, nbar_discretization);
    }; // now we have callable cachable eos

    // Possible applications follow (uncomment for usage)

    /*
    TCanvas *c = new TCanvas("c", "c", 600, 600);
    auto tov_solution = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                                const std::function<double(double)> &, double, double, double, double>(tov_solver::tov_solution);

    int n = 2000; // discretization
    std::vector<double> r(n), m(n);
    double init_density = 0.99 * constants::ist_ns::edensity_upp * constants::conversion::mev_over_fm3_gev4,
           step_density = 0.001 * constants::ist_ns::edensity_upp * constants::conversion::mev_over_fm3_gev4;
    double R_ns = tov_solution(eos, 0, init_density, 10.0 / n * constants::conversion::km_gev, step_density)[4];
    tov_solution.erase();
    for (int i = 0; i < n; ++i)
    {
        using namespace constants;

        r[i] = 1.0 * i / n * R_ns;
        m[i] = tov_solution(eos, r[i], init_density, 1.0 / n * R_ns, step_density)[0];
        r[i] /= conversion::km_gev;
        m[i] *= conversion::gev_over_msol;
    }
    TGraph *tg = new TGraph(n, &r[0], &m[0]);
    tg->Draw("ACP");
    // tg->GetXaxis()->SetLimits(1E-5, 1E-1);
    // tg->GetYaxis()->SetRangeUser(1E-8, 1E-1);
    // gPad->SetLogx();
    // gPad->SetLogy();
    c->SaveAs("pictures/mass_of_r_ist.pdf");*/

    /*
    TCanvas *c = new TCanvas("c", "c", 600, 600);
    auto tov_solution = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                                    const std::function<double(double)> &, double, double, double, double>(tov_solver::tov_solution);
    int n = 2000; // discretization
    std::vector<double> rho(n), p(n);
    double init_density = 0.999 * constants::ist_ns::edensity_upp * constants::conversion::mev_over_fm3_gev4,
           step_density = 1 * constants::ist_ns::edensity_low * constants::conversion::mev_over_fm3_gev4;
    double R_ns = tov_solution(eos, 0, init_density, 10.0 / n * constants::conversion::km_gev, step_density)[4];
    tov_solution.erase();
    for (int i = 0; i < n; ++i)
    {
        using namespace constants;
        rho[i] = tov_solution(eos, 1.0 * i / n * R_ns, init_density, 1.0 / n * R_ns, step_density)[1];
        p[i] = tov_solution(eos, 1.0 * i / n * R_ns, init_density, 1.0 / n * R_ns, step_density)[3];
    }
    TGraph *tg = new TGraph(n, &rho[0], &p[0]);
    tg->Draw("ACP");
    tg->GetXaxis()->SetLimits(constants::conversion::mev_over_fm3_gev4 * constants::ist_ns::edensity_low, constants::conversion::mev_over_fm3_gev4 * constants::ist_ns::edensity_upp);
    tg->GetYaxis()->SetRangeUser(constants::conversion::mev_over_fm3_gev4 * constants::ist_ns::pressure_low, constants::conversion::mev_over_fm3_gev4 * constants::ist_ns::pressure_upp);
    gPad->SetLogx();
    gPad->SetLogy();
    c->SaveAs("pictures/p_of_rho_ist.pdf");*/

    /*
    TCanvas *c = new TCanvas("c", "c", 600, 600);
    std::ifstream istmr("data/M-R_IST_NS.TXT");
    std::vector<double> rdat, mdat, r, m;
    std::string line;                            // variable created for reading from file
    std::getline(istmr >> std::ws, line);        // skip 1 line
    while (std::getline(istmr >> std::ws, line)) // std::ws for omission of external whitespaces
    {
        for (int i = 0; i < line.length() - 1; ++i)
        {
            if (line[i] <= 32 && line[i + 1] <= 32) // all symbols below or equal to 32 are technical or whitespace
            {
                line.erase(i, 1);
                --i;
            }
        } // strip whitespaces
        std::stringstream strstr(line);
        std::string word;
        std::getline(strstr, word, ' '); // contains density
        std::getline(strstr, word, ' '); // contains r
        rdat.push_back(std::stod(word));
        std::getline(strstr, word, ' '); // contains m
        mdat.push_back(std::stod(word));
    }
    TGraph *gr1 = new TGraph(rdat.size(), &rdat[0], &mdat[0]);
    auto tov_solution = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                                    const std::function<double(double)> &, double, double, double, double>(tov_solver::tov_solution);
    // for (double init_density = 0.01 * constants::apr4::edensity_upp; init_density < constants::apr4::edensity_upp; init_density += 0.01 * constants::apr4::edensity_upp)
    for (double init_density = 0.008 * constants::conversion::mev_over_fm3_gev4 * constants::ist_ns::edensity_upp;
         init_density < 1 * constants::conversion::mev_over_fm3_gev4 * constants::ist_ns::edensity_upp;
         init_density += 0.001 * constants::conversion::mev_over_fm3_gev4 * constants::ist_ns::edensity_upp)
    {
        r.push_back(tov_solution(eos, 0, init_density, 1.0 / 50 * constants::conversion::km_gev, 0.0000001 * constants::conversion::mev_over_fm3_gev4 * constants::ist_ns::edensity_upp)[4]);
        m.push_back(tov_solution(eos, r.back(), init_density, 1.0 / 50 * constants::conversion::km_gev, 0.0000001 * constants::conversion::mev_over_fm3_gev4 * constants::ist_ns::edensity_upp)[0]);
        r.back() /= constants::conversion::km_gev;        // in km
        m.back() *= constants::conversion::gev_over_msol; // in M_sol
    }
    TGraph *gr2 = new TGraph(r.size(), &r[0], &m[0]);
    TMultiGraph *jogr = new TMultiGraph();
    jogr->Add(gr1, "2CP");
    jogr->Add(gr2, "2CP");

    jogr->SetTitle("M-R relation for IST NS");
    jogr->GetXaxis()->SetTitle("r [km]");
    jogr->GetYaxis()->SetTitle("M [Msol]");
    jogr->Draw("AP");
    c->SaveAs("pictures/mass-radius-diagram_ist.pdf");*/

    /*
    TCanvas *c = new TCanvas("c", "c", 800, 800);
    int size = 2000;
    std::vector<std::vector<double>> x(3, std::vector<double>(size));		// coordinates
    std::vector<std::vector<double>> ephi(3, std::vector<double>(size));	// e^phi time metric functions
    std::vector<std::vector<double>> elambda(3, std::vector<double>(size)); // e^lambda space metric functions
    std::vector<double> init_densities({3.8553E-3, 4.945E-3, 7.035E-3});	// initial densities, GeV^4 ; correspond to {1.2,1.6,2.0} M_sun

    std::vector<double> colors({4, 3, 2}); // used for different init_densities
    std::vector<double> styles({1, 9});	   // used for different metric functions
    auto tov_solution = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                                    const std::function<double(double)> &, double, double, double, double>(tov_solver::tov_solution);
    for (int i = 0; i < init_densities.size(); ++i)
    {
        double R_ns = tov_solution(eos, 0, init_densities[i], 1.0 / size * 5E19, 1E-6)[4];
        // std::cout << tov_solution(eos, R_ns, init_densities[i], 1.0 / size * 5E19, 1E-6)[0] * 1.79 / 2 * 1E-57 << '\n';
        tov_solution.erase();
        double G = 6.709E-39; // gravitational constant, natural units
        for (int j = 0; j < size; ++j)
        {
            x[i][j] = static_cast<double>(j + 1) / size * R_ns;
            ephi[i][j] = TMath::Exp(tov_solution(eos, x[i][j], init_densities[i], 1.0 / size * R_ns, 1E-6)[2]);
            elambda[i][j] = TMath::Power(1 - 2 * G * tov_solution(eos, x[i][j], init_densities[i], 1.0 / size * R_ns, 1E-6)[0] / x[i][j], -1.0 / 2);
            // std::cout << elambda[i][j] << '\n';
            x[i][j] *= 2E-19 / 11.46; // in (11.46 km) units, which is R_ns at M_ns = 2.0 M_sun
        }
    } // data filled

    TGraph ***metric_functions_graphs = new TGraph **[2];
    for (int i = 0; i < 2; ++i)
        metric_functions_graphs[i] = new TGraph *[3];
    for (int i = 0; i < 2; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            if (i == 0)
                metric_functions_graphs[i][j] = new TGraph(x[j].size(), &x[j][0], &ephi[j][0]);
            if (i == 1)
                metric_functions_graphs[i][j] = new TGraph(x[j].size(), &x[j][0], &elambda[j][0]);
            // at this point, each graph in metric_functions_graphs is filled

            // cosmetics follows
            metric_functions_graphs[i][j]->SetLineWidth(4);
            metric_functions_graphs[i][j]->SetLineStyle(styles[i]);
            // metric_functions_graphs[i][j]->SetLineColor(colors[j]);
            metric_functions_graphs[i][j]->SetLineColorAlpha(colors[j], 0.85);
        }
    }

    TMultiGraph *joint_metric_functions = new TMultiGraph();
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 3; ++j)
            joint_metric_functions->Add(metric_functions_graphs[i][j], "2C");

    gPad->SetTicks();
    gPad->SetTopMargin(0.05);
    gPad->SetLeftMargin(0.11);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    TLegend *legend = new TLegend(0.15, 0.63, 0.6, 0.93);
    legend->SetBorderSize(0);
    legend->SetTextFont(113);
    legend->SetTextSize(27);
    legend->SetFillStyle(0);
    legend->SetMargin(0.35);
    legend->SetHeader("APR4 EoS, metric functions");
    legend->SetTextFont(43);
    legend->AddEntry(metric_functions_graphs[0][0], "exp(#Phi)", "L");
    legend->AddEntry(metric_functions_graphs[1][0], "exp(#Lambda)", "L");
    legend->SetFillStyle(0);
    TLegend *sublegend = new TLegend(0.15, 0.45, 0.9, 0.53);
    sublegend->SetBorderSize(0);
    sublegend->SetTextFont(113);
    sublegend->SetTextSize(27);
    sublegend->SetFillStyle(0);
    sublegend->SetMargin(0.35);
    sublegend->SetNColumns(3);
    sublegend->SetTextFont(43);
    sublegend->AddEntry(metric_functions_graphs[0][0], "1.2M_{#odot}", "L");
    sublegend->AddEntry(metric_functions_graphs[0][1], "1.6M_{#odot}", "L");
    sublegend->AddEntry(metric_functions_graphs[0][2], "2.0M_{#odot}", "L");
    sublegend->SetFillStyle(3001);
    joint_metric_functions->Draw("A");
    joint_metric_functions->GetYaxis()->SetLabelFont(43);
    joint_metric_functions->GetYaxis()->SetLabelSize(26);
    joint_metric_functions->GetYaxis()->SetTitleFont(43);
    joint_metric_functions->GetYaxis()->SetTitleSize(33);
    joint_metric_functions->GetYaxis()->SetTitleOffset(0.9);
    joint_metric_functions->GetXaxis()->SetLabelFont(43);
    joint_metric_functions->GetXaxis()->SetLabelSize(26);
    joint_metric_functions->GetXaxis()->SetTitleFont(43);
    joint_metric_functions->GetXaxis()->SetTitleSize(23);
    joint_metric_functions->GetXaxis()->SetTitleOffset(1.5);
    joint_metric_functions->GetXaxis()->SetTitle("#frac{r}{R_{2.0}}");
    joint_metric_functions->GetYaxis()->SetTitle("g_{00}, g_{rr}");
    joint_metric_functions->GetXaxis()->SetLimits(0, 0.9);
    sublegend->Draw();
    legend->Draw();
    c->SaveAs("pictures/various_metric_functions.pdf");*/

    // We have to find B_ij coefficients, whose definition follows B_ij = \int_0^R 4pir^2 e^{\Lambda(r)-\Phi(r)} dn_i/d\mu_j dr
    // It is not challenging to extract R, \Lambda(r), \Phi(r) from TOV solver
    // dn_i/d\mu_j requires some extra work though; I have to calculate these derivatives along constant chemical potentials and put their exact values afterwards
    //, where n_B should be a function of r (use TOV for this purpose)
    // Here n_i = Y_i n_B, what I can easily extract from EoS
    // and mu_j = q_j mu_q + B_j mu_b with mu_q = -\cbrt{3pi^2 n_B Ye} and mu_b = (rho+P)/n_B - mu_q (Yp-Ye)
    // convention : [0] = e, [1] = mu, [2] = n, [3] = p

    
    std::vector<std::vector<double>> Bij(4, std::vector<double>(4)); // 4x4 matrix that were looking for

    auto tov = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                       const std::function<double(double)> &, double, double, double, double>(tov_solver::tov_solution); // TOV solver

    double initial_density = 0.07 * constants::ist_ns::edensity_upp * constants::conversion::mev_over_fm3_gev4,
           density_step = 0.0001 * constants::conversion::mev_over_fm3_gev4 * constants::ist_ns::edensity_upp; // central energy density in GeV^4
    std::vector<double> q({-1, -1, 0, 1}), B({0, 0, 1, 1});                                                    // relative electric/baryonic charges

    TCanvas *c = new TCanvas("c", "c", 600, 600);

    // calculate nbar of r
    auto nbar = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>, double, double>(
        [&istdat, initial_density, density_step, &eos, &tov](std::vector<std::vector<double>> &cache, double r, double r_step)
        {
            // cache contains {r, n_B(r)} arrays; recaching is not supported at the moment, call ::erase instead
            // returns [0] -> nbar at given point and [1] -> radius at which crust begins
            // warning : initial_density should not be less (or even of order of) that density at core_limit; exception would be thrown otherwise
            using namespace constants;

            if (cache.empty())
            {
                double R_ns = tov(eos, 0, initial_density, r_step, density_step)[4];
                cache = std::vector<std::vector<double>>(2, std::vector<double>());
                for (double r_current = 0; r_current < R_ns; r_current += r_step)
                    cache[0].push_back(r_current);
                for (int i = 0; i < cache[0].size(); ++i)
                {
                    double r_current = cache[0][i];
                    // now we somehow have to find corresponding n_B
                    // let's stick to densities
                    double density_at_r = tov(eos, r_current, initial_density, r_step, density_step)[1];
                    // density_prec = TMath::Abs(tov(eos, r_current + r_step, initial_density, r_step, density_step)[1] - density_at_r); // we need these for bisection; maybe better use density_step
                    if (density_at_r <= conversion::mev_over_fm3_gev4 * ist_ns::edensity_core_limit)
                    { // until we reach core limit
                        // cache[1][i] = 0.0;
                        // continue;
                        break; // finish writing
                    }
                    double nbar_low = 1 * ist_ns::nbar_low, nbar_upp = 1 * ist_ns::nbar_upp; // we need these for bisection search; in fm-3 units for now
                    double nbar_mid = (nbar_low + nbar_upp) / 2.0;
                    // while (conversion::g_over_cm3_gev4 * TMath::Abs(istdat(std::vector<double>({nbar_low}))[0] - istdat(std::vector<double>({nbar_upp}))[0]) >= density_prec)
                    while (TMath::Abs(nbar_upp - nbar_low) > ist_ns::nbar_low)
                    {
                        // while we are too far from appropriate precision for nbar estimate
                        // recalculate via bisection method
                        nbar_mid = (nbar_low + nbar_upp) / 2.0;
                        if (conversion::mev_over_fm3_gev4 * istdat(std::vector<double>({nbar_mid}))[0] > density_at_r)
                            nbar_upp = nbar_mid;
                        else
                            nbar_low = nbar_mid;
                    }
                    cache[1].push_back(nbar_mid);
                }
                cache[0].resize(cache[1].size()); // truncate radii array so that to fit to nbar
            }
            return std::vector<double>({auxiliaries::interpolate(cache[0], cache[1], auxiliaries::InterpolationMode::kLinear, r), cache[0].back()});
        });

    // now we have to calculate dn_i/d\mu_j at constant mu_k=/=j. Original data does not contain this information,
    // so we need to calculate it via external source.
    std::ifstream n_i_of_mu_j("../data/DENSITIES.TXT");
    if (!n_i_of_mu_j.is_open()) // this file is large, so we better check if it is open
        throw std::runtime_error("Could not open DENSITIES.TXT");

    // a function that reads densities at (mu_n-m_N = 5i, mu_p-m_N = 5j, mu_e = 5k) [MeV] and returns them in fm-3 units
    auto n_i_of_mu_j_extractor = [&n_i_of_mu_j](size_t i, size_t j, size_t k)
    {
        // skip i + 302j + 301*302k - 1 lines
        for (size_t l = 0; l < i + j * 302 + k * 301 * 302; ++l)
            n_i_of_mu_j.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
        // read next line
        std::string line;
        std::getline(n_i_of_mu_j, line);
        // clean the line up
        line = auxiliaries::retrieve_cleared_line(line);
        // read densities from it
        std::stringstream ss(line);
        std::vector<double> output(3);
        std::string word;
        getline(ss, word, ' '); // skip chemical potentials data which we know
        getline(ss, word, ' '); // -//-
        getline(ss, word, ' '); // -//-
        getline(ss, word, ' '); // neutron density
        output[0] = std::stod(word);
        getline(ss, word, ' '); // proton density
        output[1] = std::stod(word);
        getline(ss, word, ' '); // electron density
        output[2] = std::stod(word);

        // set file reader to initial state
        n_i_of_mu_j.clear();
        n_i_of_mu_j.seekg(0, std::ios::beg);

        // return densities in fm-3 units
        return output;
    };


    using namespace constants::ist_ns;
    using namespace constants;
    for (double nbar_loc = nbar_upp / 3; nbar_loc > nbar_core_limit; nbar_loc += (nbar_upp / 10 - nbar_core_limit) / 1000.0)
    {
        double muq = istdat(std::vector<double>({nbar_loc}))[4], mub = istdat(std::vector<double>({nbar_loc}))[2];
        double n_n = istdat(std::vector<double>({nbar_loc}))[6], n_p = istdat(std::vector<double>({nbar_loc}))[7],
               n_e = istdat(std::vector<double>({nbar_loc}))[5];
        size_t i_mu = floor((mub - scientific::M_N * conversion::gev_over_mev) / 5), j_mu = floor((mub + muq - scientific::M_N * conversion::gev_over_mev) / 5), k_mu = floor(muq / 5);
        auto n_vect = n_i_of_mu_j_extractor(i_mu, j_mu, k_mu);
        //compare n_vect with densities defined above
        double n_n_ext = n_vect[0], n_p_ext = n_vect[1], n_e_ext = n_vect[2];
        std::cout << "nbar = " << nbar_loc << "; n_n = " << n_n << " against n_n_ext = " << n_n_ext << "; n_p = " << n_p << " against n_p_ext = " << n_p_ext << "; n_e = " << n_e << " against n_e_ext = " << n_e_ext << std::endl;
    }

    /*
    double R_ns = tov(eos, 0, initial_density, 1.0 / 1000 * constants::conversion::km_gev, density_step)[4];
    tov.erase();
    double Rcore = nbar(0, R_ns / 1000)[1];
    std::vector<double> r, nb;
    int n = 1000;
    for (int i = 0; i < n; ++i)
    {
        r.push_back(1.0 * i / n * Rcore);
        nb.push_back(nbar(r.back(), R_ns / 1000)[0]);
    }
    TGraph *gr = new TGraph(n, &r[0], &nb[0]);
    gr->Draw("ACP");
    c->SaveAs("pictures/nbar_of_r_ist.pdf");*/

    /*
    auto dn_over_dmu = [&istdat, &n_i_of_mu_j_extractor, q, B, initial_density, &nbar](double r, int i, int j, double r_step)
    {
        // calculates dni/dmj
        using namespace constants;

        double nb = nbar(r, r_step)[0];

        double n_low, n_upp; // GeV^3
        double muq;          // MeV
        double mub;          // MeV

        if (i == 1 || j == 1)
            return 0.0; // no muons in the EoS

        muq = istdat({nb})[4];
        mub = istdat({nb})[2];

        // now we can calculate dn/dmu over given mu_j direction
        // for this purpose, we should simply find derivative of nbar with respect to mu_j

        // first, we need to find lower bound indices for these potentials
        size_t i_mu = floor((mub - scientific::M_N * conversion::gev_over_mev) / 5), j_mu = floor((mub + muq - scientific::M_N * conversion::gev_over_mev) / 5), k_mu = floor(muq / 5);

        switch (i)
        {
        case 0:
            n_low = n_i_of_mu_j_extractor(i_mu, j_mu, k_mu)[2];
            switch (j)
            {
            case 0:
                n_upp = n_i_of_mu_j_extractor(i_mu, j_mu, k_mu + 1)[2];
                break;
            case 2:
                n_upp = n_i_of_mu_j_extractor(i_mu + 1, j_mu, k_mu)[2];
                break;
            case 3:
                n_upp = n_i_of_mu_j_extractor(i_mu, j_mu + 1, k_mu)[2];
                break;
            default:
                break;
            }
            break;
        case 2:
            n_low = n_i_of_mu_j_extractor(i_mu, j_mu, k_mu)[0];
            switch (j)
            {
            case 0:
                n_upp = n_i_of_mu_j_extractor(i_mu, j_mu, k_mu + 1)[0];
                break;
            case 2:
                n_upp = n_i_of_mu_j_extractor(i_mu + 1, j_mu, k_mu)[0];
                break;
            case 3:
                n_upp = n_i_of_mu_j_extractor(i_mu, j_mu + 1, k_mu)[0];
                break;
            default:
                break;
            }
            break;
        case 3:
            n_low = n_i_of_mu_j_extractor(i_mu, j_mu, k_mu)[1];
            switch (j)
            {
            case 0:
                n_upp = n_i_of_mu_j_extractor(i_mu, j_mu, k_mu + 1)[1];
                break;
            case 2:
                n_upp = n_i_of_mu_j_extractor(i_mu + 1, j_mu, k_mu)[1];
                break;
            case 3:
                n_upp = n_i_of_mu_j_extractor(i_mu, j_mu + 1, k_mu)[1];
                break;
            default:
                break;
            }
            break;
        default:
            break;
        }

        return (n_upp - n_low) / 5;
    };

    auto exp_lambda_minus_phi = [&tov, &eos, initial_density, density_step](double r, double r_step)
    {
        double m = tov(eos, r, initial_density, r_step, density_step)[0];
        double exp_lambda = TMath::Power(1 - 2 * constants::scientific::G * m / r, -1.0 / 2),
               phi = tov(eos, r, initial_density, r_step, density_step)[2];
        return exp_lambda * TMath::Exp(-phi);
    };

    double R_ns = tov(eos, 0, initial_density, 1.0 / 1000 * constants::conversion::km_gev, density_step)[4];
    tov.erase();
    double r_step = R_ns / 1000.0;
    double Rcore = nbar(0, R_ns / 1000)[1];
    for (int i = 0; i < 4; ++i)
    {
        for (int j = 0; j < 4; ++j)
        {
            using namespace constants;
            double r = r_step / 2;
            while (r < Rcore - 1.5 * r_step)
            {
                Bij[i][j] += r_step / 2 * 4 * scientific::Pi * (r * r * exp_lambda_minus_phi(r, r_step) * dn_over_dmu(r, i, j, r_step) + (r + r_step) * (r + r_step) * exp_lambda_minus_phi(r + r_step, r_step) * dn_over_dmu(r + r_step, i, j, r_step));
                r += r_step;
            }
            std::cout << Bij[i][j] << " ";
        }
        std::cout << "\n";
    }*/
}