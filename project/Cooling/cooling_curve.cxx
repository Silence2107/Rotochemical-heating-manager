
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
    argparse::ArgumentParser parser("cooling_curve", "Evaluates surface temperature time dependency based on EoS", "Argparse powered by SiLeader");

    parser.addArgument({"--inputfile"}, "json input file path (required)");
#if RHM_HAS_ROOT
    parser.addArgument({"--pdf_path"}, "pdf output file path (optional, default: Cooling.pdf)");
    parser.addArgument({"--rootfile_path"}, "root output file path (optional, default: None)");
#endif
    auto args = parser.parseArgs(argc, argv);

    using namespace instantiator;
    instantiator::instantiate_system(args.get<std::string>("inputfile"), {"TOV", "COOL"});

    auxiliaries::io::Logger logger(__func__);

#if RHM_HAS_ROOT
    std::string pdf_path = args.safeGet<std::string>("pdf_path", "Cooling.pdf");
    TFile *rootfile = nullptr;
    if (args.has("rootfile_path"))
        rootfile = new TFile(args.get<std::string>("rootfile_path").c_str(), "RECREATE");
#endif

    logger.log([]()
               { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
               [&]()
               { return "Instantiation complete"; });

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
        [&](std::vector<std::vector<double>> &cache, double r)
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
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp);

    auto hadron_murca_emissivity = cooling::predefined::neutrinic::hadron_murca_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp);

    auto hadron_bremsstrahlung_emissivity = cooling::predefined::neutrinic::hadron_bremsstrahlung_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, ion_volume_fr, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp);

    auto hadron_PBF_emissivity = cooling::predefined::neutrinic::hadron_pbf_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp);

    auto quark_ud_durca_emissivity = cooling::predefined::neutrinic::quark_ud_durca_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto quark_us_durca_emissivity = cooling::predefined::neutrinic::quark_us_durca_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto quark_ud_murca_emissivity = cooling::predefined::neutrinic::quark_ud_murca_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto quark_us_murca_emissivity = cooling::predefined::neutrinic::quark_us_murca_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto quark_bremsstrahlung_emissivity = cooling::predefined::neutrinic::quark_bremsstrahlung_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

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
        result += quark_ud_durca_emissivity(r, t, T) + quark_us_durca_emissivity(r, t, T) +
                  quark_ud_murca_emissivity(r, t, T) + quark_us_murca_emissivity(r, t, T) +
                  quark_bremsstrahlung_emissivity(r, t, T);
        result += electron_bremsstrahlung_emissivity(r, t, T);
        return result * exp_phi(r) * exp_phi(r);
    };

    // microscopics
    auto fermi_specific_heat_dens = auxiliaries::phys::fermi_specific_heat_density(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp, superconduct_q_gap);

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
    for (double r = cooling_radius_step / 2.0; r < r_ns; r += cooling_radius_step)
    {
        radii.push_back(r);
        profile.push_back(initial_t_profile_inf(r, r_ns, exp_phi, nbar));
    }
    double t_step = base_t_step,
           t_curr = t_init,
           temp_curr = profile.end()[-2];

    // solution arrays
    std::vector<double> time, surface_temp;
    std::vector<std::vector<double>> others(2);
    time.reserve(cooling_n_points_estimate);
    surface_temp.reserve(cooling_n_points_estimate);
    others[0].reserve(cooling_n_points_estimate);
    others[1].reserve(cooling_n_points_estimate);

    logger.log([]()
               { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
               [&]()
               {
                   double conv = 1.0E6 / (constants::conversion::myr_over_s * constants::conversion::gev_s);
                   std::stringstream ss;
                   ss << std::scientific << std::setprecision(3) << "Time[yr] array is exp mapped [" << base_t_step * conv << ", " << t_end * conv << ", " << cooling_n_points_estimate << "] with possible adaptions";
                   return ss.str();
               });

    size_t indent = 20;
    std::cout << "M = " << m_ns * constants::conversion::gev_over_msol << " [Ms]\n";
    std::cout << std::left << std::setw(indent) << "t[years] "
              << std::setw(indent) << "Te^inf[K] "
              << std::setw(indent) << "L^inf_ph[erg/s] "
              << std::setw(indent) << "L^inf_nu[erg/s] " << '\n';

    while (t_curr < t_end)
    {
        double next_T; // predicted T
        bool reached_adaption_limit = false; // control for adaptive solver
        bool reached_negative_temperature = false; // exclude NaN generation to "negative" temperature
        double neutrino_lum = 0.0; // placeholder for neutrino luminosity estimate

        // non-equilibrium stage
        if (!switch_to_equilibrium(t_curr, profile))
        {
            auto t_l_profiles = cooling::solver::nonequilibrium_cooling(
                t_curr, t_step, Q_nu, fermi_specific_heat_dens, thermal_conductivity,
                exp_lambda, exp_phi, radii, profile, te_tb, cooling_newton_step_eps, cooling_newton_max_iter);
            next_T = t_l_profiles[0].end()[-2];
            reached_adaption_limit = t_l_profiles[2][0];
            reached_negative_temperature = t_l_profiles[2][1];
            double max_diff = 0;
            for (size_t i = 0; i < radii.size() - 1; ++i)
            {
                // excluding surface point
                max_diff = std::max(max_diff, fabs(t_l_profiles[0][i] - profile[i]) / profile[i]);
            }
            if (max_diff > cooling_max_diff_per_t_step || reached_adaption_limit || reached_negative_temperature)
            {
                t_step /= 2;
                logger.log([&]()
                           { return max_diff > cooling_max_diff_per_t_step; }, auxiliaries::io::Logger::LogLevel::kDebug,
                           [&]()
                           { 
                                std::stringstream ss;
                                ss << std::scientific << std::setprecision(3) << "Target profile difference exceeded (" << max_diff << " > " << cooling_max_diff_per_t_step << "), adapting (raise CoolingSolver.StepTolerance/NewtonTolerance if this halts progress)";
                                return ss.str();
                            }, "non-eq. cooling");
                logger.log([&]()
                           { return reached_adaption_limit; }, auxiliaries::io::Logger::LogLevel::kDebug,
                           [&]()
                           { return "Adaption limit reached, adapting (raise CoolingSolver.NewtonMaxIter/NewtonTolerance if this reoccurs)"; }, "non-eq. cooling");
                logger.log([&]()
                           { return reached_negative_temperature; }, auxiliaries::io::Logger::LogLevel::kDebug,
                           [&]()
                           { return "Negative temperature reached, adapting (lower CoolingSolver.NewtonTolerance if this reoccurs)"; }, "non-eq. cooling");
                continue;
            }
            profile = t_l_profiles[0];
            // To calculate neutrino luminosity I refrain to interpolating T(r) to avoid performance hindering. I instead carefully dublicate integration procedure
            for (size_t count = 0; count < radii.size() - 1; ++count)
            {
                neutrino_lum += 4 * constants::scientific::Pi * radii[count] * radii[count] * exp_lambda(radii[count]) * (radii[count + 1] - radii[count]) * Q_nu(radii[count], t_curr + t_step, profile[count]);
            }
        }

        // equilibrium stage
        else
        {
            auto equilibrium_data = cooling::solver::equilibrium_cooling(t_curr, t_step, cooling_rhs, temp_curr, cooling_newton_step_eps, cooling_newton_max_iter);
            next_T = equilibrium_data[0];
            reached_adaption_limit = equilibrium_data[1];
            reached_negative_temperature = equilibrium_data[2];
            double max_diff = std::abs((temp_curr - next_T) / temp_curr);
            if (max_diff > cooling_max_diff_per_t_step || reached_adaption_limit || reached_negative_temperature)
            {
                t_step /= 2.0;
                logger.log([&]()
                           { return max_diff > cooling_max_diff_per_t_step; }, auxiliaries::io::Logger::LogLevel::kDebug,
                           [&]()
                           { 
                                std::stringstream ss;
                                ss << std::scientific << std::setprecision(3) << "Target temperature difference exceeded (" << max_diff << " > " << cooling_max_diff_per_t_step << "), adapting (raise CoolingSolver.StepTolerance/NewtonTolerance if this halts progress)";
                                return ss.str();
                            }, "eq. cooling");
                logger.log([&]()
                           { return reached_adaption_limit; }, auxiliaries::io::Logger::LogLevel::kDebug,
                           [&]()
                           { return "Adaption limit reached, adapting (raise CoolingSolver.NewtonMaxIter/NewtonTolerance if this reoccurs)"; }, "eq. cooling");
                logger.log([&]()
                           { return reached_negative_temperature; }, auxiliaries::io::Logger::LogLevel::kDebug,
                           [&]()
                           { return "Negative temperature reached, adapting (lower CoolingSolver.NewtonTolerance if this reoccurs)"; }, "eq. cooling");
                continue;
            }
            neutrino_lum = neutrino_luminosity(t_curr + t_step, next_T);
        }
        temp_curr = next_T;
        t_curr += t_step;
        t_step *= exp_rate_estim;

        // save in understandable units
        time.push_back(1.0E6 * t_curr / (constants::conversion::myr_over_s * constants::conversion::gev_s));
        surface_temp.push_back(auxiliaries::phys::te_tb_relation(temp_curr, r_ns, m_ns, crust_eta) * exp_phi_at_R * constants::conversion::gev_over_k);
        others[0].push_back(photon_luminosity(t_curr, temp_curr) * constants::conversion::gev_s / constants::conversion::erg_over_gev);
        others[1].push_back(neutrino_lum * constants::conversion::gev_s / constants::conversion::erg_over_gev);

        // print
        std::cout << std::left << std::setw(indent) << time.back() << std::setw(indent) << surface_temp.back() << std::setw(indent) << others[0].back() << std::setw(indent) << others[1].back() << '\n';
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

    auto gr = new TGraph(time.size(), time.data(), surface_temp.data());
    gr->GetXaxis()->SetTitle("t [yr]");
    gr->GetYaxis()->SetTitle("T^{#infty}_{s} [K]");
    if (rootfile)
    {
        auto gr_l_gamma = new TGraph(time.size(), time.data(), others[0].data());
        gr_l_gamma->GetXaxis()->SetTitle("t [yr]");
        gr_l_gamma->GetYaxis()->SetTitle("L^{#infty}_{#gamma} [erg/s]");
        auto gr_l_nu = new TGraph(time.size(), time.data(), others[1].data());
        gr_l_nu->GetXaxis()->SetTitle("t [yr]");
        gr_l_nu->GetYaxis()->SetTitle("L^{#infty}_{#nu} [erg/s]");
        rootfile->cd();
        rootfile->WriteObject(gr, "cooling_curve");
        rootfile->WriteObject(gr_l_gamma, "l_gamma");
        rootfile->WriteObject(gr_l_nu, "l_nu");
        rootfile->Close();
    }
    gr->SetLineColor(kBlue);
    gr->SetLineWidth(2);
    gr->SetLineStyle(1);
    gr->Draw("AL");
    gr->GetYaxis()->SetTitleOffset(1.5);
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
    gr->GetYaxis()->SetLimits(surface_temp.front(), surface_temp.back());
    gr->GetXaxis()->SetLimits(time.front(), time.back());

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