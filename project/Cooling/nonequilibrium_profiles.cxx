
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
#include <sstream>
#include <fstream>

#if RHM_HAS_ROOT
#include <TCanvas.h>
#include <TGraph.h>
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TFile.h>
#endif

int main(int argc, char **argv)
{
    std::string program_name = "nonequilibrium_profiles";
    argparse::ArgumentParser parser(program_name, "Evaluates radial temperature profiles with time based on EoS", "Argparse powered by SiLeader");

    parser.addArgument({"--inputfile"}, "json input file path (required)");
#if RHM_HAS_ROOT
    parser.addArgument({"--pdf_path"}, "pdf output file path (optional, default: CoolingProfiles.pdf)");
    parser.addArgument({"--rootfile_path"}, "root output file path (optional, default: None)");
#endif
    parser.addArgument({"--write_exp"}, "multiplier between consecutive write times (optional, default: 10.0)");
    parser.addArgument({"--no_intermediate_print"}, "whether to print minimal information at each time step (optional, value-free, default: print)", argparse::ArgumentType::StoreTrue);

    auto args = parser.parseArgs(argc, argv);

    using namespace instantiator;
    instantiator::instantiate_system(args.get<std::string>("inputfile"), {"TOV", "COOL"});

    auxiliaries::io::Logger logger(program_name);

#if RHM_HAS_ROOT
    std::string pdf_path = args.safeGet<std::string>("pdf_path", "CoolingProfiles.pdf");
    TFile *rootfile = nullptr;
    if (args.has("rootfile_path"))
        rootfile = new TFile(args.get<std::string>("rootfile_path").c_str(), "RECREATE");
#endif

    double write_time_expansion = args.safeGet<double>("write_exp", 10.0);
    bool print_all_time = !args.has("no_intermediate_print");

    logger.log([]()
               { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
               [&]()
               { return "Instantiation complete"; });

    // RUN --------------------------------------------------------------------------

    // EoS definition

    auto eos_inv_cached = edensity_of_pressure;

    // TOV solver

    auto tov_cached = auxiliaries::math::CachedFunc<std::vector<auxiliaries::math::Interpolator>, std::vector<double>,
                                                    const std::function<double(double)> &, double, double, double, double,
                                                    double, size_t, auxiliaries::math::Interpolator::InterpolationMode>(tov_solver::tov_solution);
    auto tov = [&tov_cached, &eos_inv_cached](double r)
    {
        // TOV solution cached
        return tov_cached(eos_inv_cached, r, center_pressure, radius_step, surface_pressure, pressure_low, tov_adapt_limit, radial_interp_mode);
    };

    auto nbar = [&](double r)
    {
        return nbar_of_pressure(tov(r)[3]);
    };

    // Cooling

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

    // photon luminosity related quantities
    auto te_tb = [&r_ns, &m_ns](double T_binf)
    {
        return auxiliaries::phys::te_tb_relation(T_binf, r_ns, m_ns, crust_eta);
    };

    // neutrino luminosity related quantities
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

    // macroscopics
    auto fermi_specific_heat_dens = auxiliaries::phys::fermi_specific_heat_density(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp, superconduct_q_gap);

    auto thermal_conductivity = auxiliaries::phys::thermal_conductivity_FI(energy_density_of_nbar,
                                                                           nbar, exp_phi);

    // solve cooling equations
    double exp_phi_at_R = pow(1 - 2 * constants::scientific::G * m_ns / r_ns, 0.5);

    // tabulate initial profile and radii
    std::vector<double> radii, profile;
    for (double r = cooling_radius_step / 2.0; r < r_ns; r += cooling_radius_step)
    {
        radii.push_back(r);
        profile.push_back(initial_t_profile_inf(r, r_ns, exp_phi, nbar));
    }

    // solution arrays (in understandable units)
    std::vector<std::vector<double>> saved_profiles;
    std::vector<double> saved_times;
    std::vector<double> saved_radii;
    for (size_t i = 0; i < radii.size(); ++i)
    {
        saved_radii.push_back(radii[i] / constants::conversion::km_gev);
    }
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

    double t_curr = 0, time_step = base_t_step;
    double write_time = base_t_step;
    size_t iter_count = 0;

    if (print_all_time)
    {
        std::cout << std::left << std::setw(indent) << "t[years]" << std::setw(indent) << "Te^inf[K]" << std::setw(indent) << "T^inf(r=0)[K]" << std::setw(indent) << "T^inf(r=R)[K]" << '\n';
    }

    while (t_curr + time_step < t_end)
    {
        using namespace constants::conversion;

        bool reached_adaption_limit = false;       // control for adaptive solver
        bool reached_negative_temperature = false; // exclude NaN generation to "negative" temperature
        // update
        while (true)
        {
            auto t_l_profiles = cooling::solver::nonequilibrium_cooling(
                t_curr, time_step, Q_nu, fermi_specific_heat_dens, thermal_conductivity,
                exp_lambda, exp_phi, radii, profile, te_tb, cooling_newton_step_eps, cooling_newton_max_iter);
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
                time_step /= 2;
                logger.log([&]()
                           { return max_diff > cooling_max_diff_per_t_step; }, auxiliaries::io::Logger::LogLevel::kDebug,
                           [&]()
                           {
                               std::stringstream ss;
                               ss << std::scientific << std::setprecision(3) << "At t = " << 1.0E6 * t_curr / (constants::conversion::myr_over_s * constants::conversion::gev_s) << " [yr] target profile difference exceeded(" << max_diff << " > " << cooling_max_diff_per_t_step << "), adapting (raise CoolingSolver.StepTolerance/NewtonTolerance if this halts progress)";
                               return ss.str();
                           },
                           "non-eq. cooling");
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
            break;
        }

        t_curr += time_step;
        time_step *= exp_rate_estim;
        logger.log([&]()
                   { return true; }, auxiliaries::io::Logger::LogLevel::kTrace,
                   [&]()
                   { 
                        std::stringstream ss;
                        ss << std::scientific << std::setprecision(3) << std::to_string(iter_count + 1) + " counts past" << " with t = " << t_curr * 1E6 / (myr_over_s * gev_s) << " [yr]";
                        return ss.str(); },
                   "T(r, t) save loop");
        ++iter_count;

        if (print_all_time)
        {
            std::cout << std::left << std::setw(indent) << t_curr * 1E6 / (myr_over_s * gev_s) << std::setw(indent) << te_tb(profile.end()[-2]) * exp_phi_at_R * gev_over_k << std::setw(indent) << profile[0] * gev_over_k << std::setw(indent) << profile.end()[-2] * gev_over_k << '\n';
        }

        if (t_curr >= write_time)
        {
            using namespace constants::conversion;
            saved_profiles.push_back(profile);
            for (size_t i = 0; i < radii.size(); ++i)
            {
                saved_profiles.back()[i] *= gev_over_k;
            }
            write_time *= write_time_expansion;
            saved_times.push_back(t_curr * 1E6 / (myr_over_s * gev_s));

            logger.log([&]()
                       { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
                       [&]()
                       {
                           std::stringstream ss;
                           ss << std::scientific << std::setprecision(3) << "Profile saved at t = " << saved_times.back() << " [yr]";
                           return ss.str();
                       },
                       "T(r, t) save loop");
        }
    }

    // print accumulated profiles
    std::cout << std::left << std::setw(indent) << "r[km]";
    for (size_t i = 0; i < saved_times.size(); ++i)
    {
        std::stringstream ss;
        ss << std::scientific << std::setprecision(3) << saved_times[i] << "[yr]";
        std::cout << std::setw(indent) << ss.str();
    }
    std::cout << '\n';
    for (size_t i = 0; i < saved_radii.size(); ++i)
    {
        std::cout << std::left << std::setw(indent) << saved_radii[i];
        for (size_t j = 0; j < saved_times.size(); ++j)
        {
            std::cout << std::setw(indent) << saved_profiles[j][i];
        }
        std::cout << '\n';
    }

#if RHM_HAS_ROOT
    // draw
    TCanvas *c1 = new TCanvas("c1", "c1");
    TMultiGraph *mg = new TMultiGraph("mg", "mg");

    mg->GetXaxis()->SetTitle("r [km]");
    mg->GetYaxis()->SetTitle("T^{#infty} [K]");
    // title offset
    mg->GetYaxis()->SetTitleOffset(1.5);
    mg->GetXaxis()->SetLimits(0, r_ns / constants::conversion::km_gev);
    mg->GetYaxis()->SetRangeUser(5.5E6, 5.5E9);
    for (size_t i = 0; i < saved_profiles.size(); ++i)
    {
        auto gr = new TGraph(saved_radii.size(), saved_radii.data(), saved_profiles[i].data());
        gr->SetLineColor(i + 1);
        gr->SetLineWidth(2.5);
        mg->Add(gr, "L");
    }
    if (rootfile)
    {
        rootfile->cd();
        rootfile->WriteObject(mg, "profiles");
        rootfile->Close();
    }
    mg->Draw("A");
    gPad->SetLogy();

    c1->SaveAs(pdf_path.c_str());
#endif
}