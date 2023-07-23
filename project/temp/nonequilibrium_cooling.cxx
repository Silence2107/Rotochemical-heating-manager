
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
#include <TMultiGraph.h>
#include <TAxis.h>
#include <TLegend.h>
#include <TFile.h>

int main(int argc, char **argv)
{
    if (argc == 1)
    {
        std::cout << "Usage: " << argv[0] << " <pdf_path=Cooling.pdf> <rootfile_path=None>" << std::endl;
    }
    std::string pdf_path = (argc > 1) ? argv[1] : "Cooling.pdf";
    bool rootfile_creation = (argc > 2);
    using namespace inputfile;

    // RUN --------------------------------------------------------------------------

    // EoS definition

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
    auto hadron_durca_emissivity = cooling::predefined::neutrinic::hadron_durca_emissivity(k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_core_limit, exp_phi, superfluid_n_1s0, superfluid_p_1s0, superfluid_n_3p2, superfluid_p_temp, superfluid_n_temp);

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

    // macroscopics
    auto fermi_specific_heat_dens = auxiliaries::phys::fermi_specific_heat_density(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_core_limit, exp_phi, superfluid_n_1s0,
        superfluid_p_1s0, superfluid_n_3p2, superfluid_p_temp, superfluid_n_temp, superconduct_q_gap);

    auto thermal_conductivity = auxiliaries::phys::thermal_conductivity_FI(energy_density_of_nbar,
                                                                           nbar, exp_phi);

    auto initial_profile = [&nbar, &exp_phi](double r)
    {
        using namespace constants::conversion;
        /*if (nbar(r) > nbar_core_limit)
            return 1E10 * exp_phi(r) / gev_over_k;
        else if (nbar(r) > nbar_crust_limit)
            return 8E9 * exp_phi(r) / gev_over_k;
        else*/
        return 5E9 / gev_over_k;
    };

    // tabulate initial profile and radii
    double cooling_radius_step = 10 * radius_step;
    std::vector<double> radii, profile;
    for(double r = cooling_radius_step / 2.0; r < r_ns; r += cooling_radius_step)
    {
        radii.push_back(r);
        profile.push_back(initial_profile(r));
    }

    std::vector<std::vector<double>> xs, ys;
    std::cout << "M/Msol " << m_ns * constants::conversion::gev_over_msol << std::endl;

    double t_curr = 0, time_step = base_t_step;
    double write_time = base_t_step;
    while (t_curr + time_step < t_end)
    {
        std::cout << "t = " << t_curr << ", prof[0] = " << profile[0] << ", prof[-1] = " << profile.back() << '\n';

        if (t_curr >= write_time)
        {
            using namespace constants::conversion;
            xs.push_back(radii);
            ys.push_back(profile);
            for (size_t i = 0; i < radii.size(); ++i)
            {
                xs.back()[i] /= km_gev;
                ys.back()[i] *= gev_over_k;
            }
            write_time *= 10.0;
        }

        bool error = false;

        // update
        while(true)
        {
            try
            {
                auto new_profile = cooling::solver::nonequilibrium_cooling_2(
                    t_curr, t_curr + time_step, Q_nu, fermi_specific_heat_dens, thermal_conductivity,
                    exp_lambda, exp_phi, radii, profile, te_tb);
                double max_diff = 0;
                for (size_t i = 0; i < radii.size(); ++i)
                {
                    max_diff = std::max(max_diff, fabs(new_profile[i] - profile[i])/profile[i]);
                    profile[i] = new_profile[i];
                }
                std::cout << "max_diff = " << max_diff << '\n';
                if (max_diff > 0.05)
                    {
                        //time_step *= 0.5;
                        //exp_rate_estim = sqrt(exp_rate_estim);
                        //std::cout << "Adapting time step \n";
                    }
                break;
            }
            catch (std::runtime_error &e)
            {
                std::cout << e.what() << '\n';
                using namespace constants::conversion;
                xs.push_back(radii);
                ys.push_back(profile);
                for (size_t i = 0; i < radii.size(); ++i)
                {
                    xs.back()[i] /= km_gev;
                    ys.back()[i] *= gev_over_k;
                }
                error = true;
                break;
            }
        }
        
        if(error)
            break;

        t_curr += time_step;
        time_step *= exp_rate_estim;
    }

    // draw
    TCanvas *c1 = new TCanvas("c1", "c1");
    TMultiGraph *mg = new TMultiGraph("mg", "mg");
    for(size_t i = 0; i < xs.size(); ++i)
    {
        auto gr = new TGraph(xs[i].size(), xs[i].data(), ys[i].data());
        gr->SetLineColor(i + 1);
        gr->SetLineWidth(2.5);
        mg->Add(gr, "L");
    }
    if (rootfile_creation)
    {
        std::string rootfile_path = argv[2];
        TFile *f = new TFile(rootfile_path.c_str(), "RECREATE");
        mg->Write();
        f->Close();
    }
    mg->Draw("A");
    // title offset
    mg->GetYaxis()->SetTitleOffset(1.5);
    mg->GetYaxis()->SetRangeUser(6E6, 6E9);
    gPad->SetLogy();

    mg->GetXaxis()->SetTitle("r [km]");
    mg->GetYaxis()->SetTitle("T^{#infty} [K]");

    c1->SaveAs(pdf_path.c_str());
}