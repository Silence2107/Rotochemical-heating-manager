
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
    argparse::ArgumentParser parser("cooling_curve", "Evaluates temperature-time dependency based on EoS", "Argparse powered by SiLeader");

    parser.addArgument({"--inputfile"}, "json input file path (optional)");
    #if RHM_HAS_ROOT
    parser.addArgument({"--pdf_path"}, "pdf output file path (optional, default: Cooling.pdf)");
    parser.addArgument({"--rootfile_path"}, "root output file path (optional, default: None)");
    #endif
    auto args = parser.parseArgs(argc, argv);

    using namespace instantiator;
    if (args.has("inputfile"))
        instantiator::instantiate_system(args.get<std::string>("inputfile"));

    #if RHM_HAS_ROOT
    std::string pdf_path = args.safeGet<std::string>("pdf_path", "Cooling.pdf");
    TFile *rootfile = nullptr;
    if (args.has("rootfile_path"))
        rootfile = new TFile(args.get<std::string>("rootfile_path").c_str(), "RECREATE");
    #endif

    // RUN --------------------------------------------------------------------------

    // EoS definition

    auto eos_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, double, double>(
        [&](std::vector<std::vector<double>> &cache, double rho)
        {
            if (rho < 0 || rho > edensity_upp)
                THROW(std::runtime_error, "Data request out of range.");
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

    // TOV solver

    auto tov_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                                    const std::function<double(double)> &, double, double, double, double>(tov_solver::tov_solution);
    auto tov = [&tov_cached, &eos_cached](double r)
    {
        // TOV solution cached
        return tov_cached(eos_cached, r, center_density, radius_step, density_step);
    };

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

    auto exp_phi = [&tov](double r)
    {
        return std::exp(tov(r)[2]);
    };

    auto exp_lambda = [&tov](double r)
    {
        return pow(1 - 2 * constants::scientific::G * tov(r)[0] / r, -0.5);
    };

    // cooling settings

    auto te_tb = [&r_ns, &m_ns](double T_binf)
    {
        return auxiliaries::phys::te_tb_relation(T_binf, r_ns, m_ns, crust_eta);
    };

    // internal emission
    auto hadron_durca_emissivity = cooling::predefined::neutrinic::hadron_durca_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_core_limit, exp_phi, superfluid_n_1s0,
        superfluid_p_1s0, superfluid_n_3p2, superfluid_p_temp, superfluid_n_temp);

    auto hadron_murca_emissivity = cooling::predefined::neutrinic::hadron_murca_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_core_limit, exp_phi, superfluid_n_1s0,
        superfluid_p_1s0, superfluid_n_3p2, superfluid_p_temp, superfluid_n_temp);

    auto hadron_bremsstrahlung_emissivity = cooling::predefined::neutrinic::hadron_bremsstrahlung_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, ion_volume_fr, nbar_core_limit, exp_phi, superfluid_n_1s0,
        superfluid_p_1s0, superfluid_n_3p2, superfluid_p_temp, superfluid_n_temp);

    auto hadron_PBF_emissivity = cooling::predefined::neutrinic::hadron_pbf_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_core_limit, exp_phi, superfluid_n_1s0,
        superfluid_p_1s0, superfluid_n_3p2, superfluid_p_temp, superfluid_n_temp);

    bool has_quarks = false;
    for (auto it = m_stars_of_nbar.begin(); it != m_stars_of_nbar.end(); ++it)
    {
        if (it->first.classify() == auxiliaries::phys::Species::ParticleClassification::kQuark)
        {
            has_quarks = true;
            break;
        }
    }

    auto quark_durca_emissivity = (has_quarks ? cooling::predefined::neutrinic::quark_durca_emissivity(
                                                    k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap)
                                              : [](double, double, double)
                                       { return 0.0; });

    auto quark_murca_emissivity = (has_quarks ? cooling::predefined::neutrinic::quark_murca_emissivity(
                                                    k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap)
                                              : [](double, double, double)
                                       { return 0.0; });

    auto quark_bremsstrahlung_emissivity = (has_quarks ? cooling::predefined::neutrinic::quark_bremsstrahlung_emissivity(
                                                             k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap)
                                                       : [](double, double, double)
                                                { return 0.0; });

    auto electron_bremsstrahlung_emissivity = cooling::predefined::neutrinic::electron_bremsstrahlung_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi);

    auto Q_nu = [&](double r, double t, double T)
    {
        using namespace constants::scientific;
        double result = 0;

        for (auto it = m_stars_of_nbar.begin(); it != m_stars_of_nbar.end(); ++it)
        {
            auto key = it->first;
            if (key.classify() == auxiliaries::phys::Species::ParticleClassification::kLepton)
            {
                result += hadron_murca_emissivity(r, key, t, T);
                result += hadron_durca_emissivity(r, key, t, T);
            }
            if (key.classify() == auxiliaries::phys::Species::ParticleClassification::kBaryon)
            {
                result += hadron_PBF_emissivity(r, key, t, T);
            }
        }
        result += hadron_bremsstrahlung_emissivity(r, t, T);
        result += quark_durca_emissivity(r, t, T) + quark_murca_emissivity(r, t, T) + quark_bremsstrahlung_emissivity(r, t, T);
        result += electron_bremsstrahlung_emissivity(r, t, T);
        return result * exp_phi(r) * exp_phi(r);
    };

    // microscopics
    auto fermi_specific_heat_dens = auxiliaries::phys::fermi_specific_heat_density(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_core_limit, exp_phi, superfluid_n_1s0,
        superfluid_p_1s0, superfluid_n_3p2, superfluid_p_temp, superfluid_n_temp, superconduct_q_gap);

    auto thermal_conductivity = auxiliaries::phys::thermal_conductivity_FI(energy_density_of_nbar,
                                                                           nbar, exp_phi);

    // equilibrium cooling settings

    // photon luminosity
    auto photon_luminosity = cooling::predefined::photonic::surface_luminosity(r_ns, m_ns, crust_eta);

    // neutrino luminosity
    auto neutrino_luminosity = auxiliaries::math::integrate_volume<double, double>(
        std::function<double(double, double, double)>(Q_nu), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p, radius_step);

    auto heat_capacity = auxiliaries::math::integrate_volume<double, double>(
        std::function<double(double, double, double)>(fermi_specific_heat_dens), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p, radius_step);

    auto cooling_rhs = [&heat_capacity, &photon_luminosity, &neutrino_luminosity](double t, double T)
    {
        // std::cout << t << " " << photon_luminosity(t, T) << " " << neutrino_luminosity(t, T) << " " << heat_capacity(t, T) << '\n';
        return -(photon_luminosity(t, T) + neutrino_luminosity(t, T)) / heat_capacity(t, T);
    };

    // solve cooling equations
    double exp_phi_at_R = pow(1 - 2 * constants::scientific::G * m_ns / r_ns, 0.5);

    // tabulate initial profile and radii
    std::vector<double> radii, profile;
    for(double r = cooling_radius_step / 2.0; r < r_ns; r += cooling_radius_step)
    {
        radii.push_back(r);
        profile.push_back(initial_t_profile_inf(r, exp_phi_at_R));
    }

    // plot the solution (assumes exp_rate_estim > 1)
    std::vector<double> x;
    std::vector<double> y;
    x.reserve(cooling_n_points_estimate);
    y.reserve(cooling_n_points_estimate);
    std::cout << "M/Msol " << m_ns * constants::conversion::gev_over_msol << std::endl;
    std::cout << "t [years] "
              << "\tTe^inf [K] "
              << "\tL_ph [erg/s] "
              << "\tL_nu [erg/s] " << std::endl;
    x.push_back(t_init);
    y.push_back(profile.end()[-2]);
    double t_step = base_t_step;

    while (x.back() < t_end)
    {
        double next_T; // predicted T

        // non-equilibrium stage
        if (!switch_to_equilibrium(x.back(), profile))
        {
            auto t_l_profiles = cooling::solver::nonequilibrium_cooling(
                x.back(), t_step, Q_nu, fermi_specific_heat_dens, thermal_conductivity,
                exp_lambda, exp_phi, radii, profile, te_tb, cooling_newton_step_eps, cooling_newton_max_iter);
            next_T = t_l_profiles[0].end()[-2];
            double max_diff = std::abs((y.back() - next_T) / y.back());
            if (max_diff > 0.05)
            {
                t_step /= 2.0;
                continue;
            }
            profile = t_l_profiles[0];
        }

        // equilibrium stage
        else
        {
            next_T = cooling::solver::equilibrium_cooling(x.back(), t_step, cooling_rhs, y.back(), cooling_newton_step_eps, cooling_newton_max_iter);
            double max_diff = std::abs((y.back() - next_T) / y.back());
            if (max_diff > 0.05)
            {
                t_step /= 2.0;
                continue;
            }
        }
        y.push_back(next_T);
        x.push_back(x.back() + t_step);
        t_step *= exp_rate_estim;

        // print in understandable units
        std::cout << 1.0E6 * x.back() / (constants::conversion::myr_over_s * constants::conversion::gev_s) << "\t" << auxiliaries::phys::te_tb_relation(y.back(), r_ns, m_ns, crust_eta) * exp_phi_at_R * constants::conversion::gev_over_k << "\t" << photon_luminosity(x.back(), y.back()) * constants::conversion::gev_s / constants::conversion::erg_over_gev << "\t" << neutrino_luminosity(x.back(), y.back()) * constants::conversion::gev_s / constants::conversion::erg_over_gev << "\t" << '\n';
    }
    // rescale to same units as in the printout
    for (size_t i = 0; i < x.size(); ++i)
    {
        x[i] *= 1.0E6 / (constants::conversion::myr_over_s * constants::conversion::gev_s);
        y[i] = auxiliaries::phys::te_tb_relation(y[i], r_ns, m_ns, crust_eta) * exp_phi_at_R * constants::conversion::gev_over_k;
    }

    #if RHM_HAS_ROOT
    // draw
    TCanvas *c1 = new TCanvas("c1", "c1");
    gPad->SetLogy();
    gPad->SetLogx();
    gPad->SetTicks();
    gPad->SetTopMargin(0.05);
    gPad->SetLeftMargin(0.11);
    gPad->SetRightMargin(0.05);
    gPad->SetBottomMargin(0.1);
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);

    auto gr = new TGraph(x.size(), x.data(), y.data());
    if (rootfile)
    {
        gr->Write();
        rootfile->Close();
    }
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(2);
    gr->SetLineStyle(1);
    gr->Draw("AL");
    gr->GetYaxis()->SetTitleOffset(1.5);
    gr->GetXaxis()->SetTitle("t [yr]");
    gr->GetYaxis()->SetTitle("T^{#infty}_{s} [K]");
    gr->GetYaxis()->SetLabelFont(43);
    gr->GetYaxis()->SetLabelSize(22);
    gr->GetYaxis()->SetTitleFont(43);
    gr->GetYaxis()->SetTitleSize(26);
    gr->GetYaxis()->SetTitleOffset(0.5);
    gr->GetXaxis()->SetLabelFont(43);
    gr->GetXaxis()->SetLabelSize(22);
    gr->GetXaxis()->SetTitleFont(43);
    gr->GetXaxis()->SetTitleSize(26);
    gr->GetXaxis()->SetTitleOffset(0.9);
    gr->GetYaxis()->SetRangeUser(7e2, 7e6);
    gr->GetXaxis()->SetLimits(1e-12, 1e7);

    auto legend = new TLegend(0.15, 0.1, 0.43, 0.38);
    legend->AddEntry(gr, "RH Manager", "l");
    legend->SetBorderSize(0);
    legend->SetTextFont(43);
    legend->SetTextSize(27);
    legend->SetFillStyle(0);
    legend->SetMargin(0.35);

    legend->Draw();

    c1->SaveAs(pdf_path.c_str());
    #endif
}