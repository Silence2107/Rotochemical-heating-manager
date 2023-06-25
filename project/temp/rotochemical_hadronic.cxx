
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

    auto eos_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>, double, double>(
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

    auto tov_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                              const std::function<double(double)> &, double, double, double, double>(tov_solver::tov_solution);
    auto tov = [&tov_cached, &eos_cached](double r)
    {
        // TOV solution cached
        return tov_cached(eos_cached, r, center_density, radius_step, density_step);
    };

    double r_crust;
    bool crust_found = false;

    auto nbar = auxiliaries::CachedFunc<std::vector<std::vector<double>>, double, double>(
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
    
    // Prerequisites for rotochemical heating

    // We need Z_npe, I_e, \Omega evolution and cooling eq. etc to perform test run

    // Z_npe is given (for now) by estimate for B_ee, in GeV
    double z_npe = 2E-60 * constants::conversion::erg_over_gev;

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

    double omega_k_sqr = pow(2.0 / 3, 3.0) * constants::scientific::G * m_ns / (r_ns * r_ns * r_ns);
    auto integrand = [&](double r)
    {
        auto Y_e = Y_i_functions_of_nbar[constants::scientific::electron];
        return -nbar(r) * 1.0 / omega_k_sqr * (Y_e(nbar(r + radius_step)) - Y_e(nbar(r))) / (tov(r + radius_step)[3] / tov(r)[3] - 1);
    };
    // I_e is given by the following integral
    double I_e = auxiliaries::integrate_volume<>(std::function<double(double)>(integrand), 0.0, r_crust, exp_lambda, auxiliaries::IntegrationMode::kGaussLegendre_12p)();

    // Dipole ansatz for \Omega evolution
    double p_0 = 1.0E-3 * constants::conversion::gev_s,
           p_0_dot = 1.0E-20;
    auto omega_omega_dot = [p_0, p_0_dot](double t)
    {
        using namespace constants::scientific;
        double braking_index = 3; // 3 for pure magnetic dipole radiation; must not be 1 in my ansatz
        return -(4 * Pi * Pi * p_0_dot) / pow(p_0, 3) *
               pow(1 + (braking_index - 1) * p_0_dot / p_0 * t, (1 + braking_index) / (1 - braking_index));
    };

    // RH control functions
    auto enhance_durca = [](double xi)
    {
        using constants::scientific::Pi;
        return 1 + 1071.0 / 457 * pow(xi / Pi, 2.0) + 315.0 / 457 * pow(xi / Pi, 4.0) + 21.0 / 457 * pow(xi / Pi, 6.0);
    };
    auto enhance_murca = [](double xi)
    {
        using constants::scientific::Pi;
        return 1 + 22020.0 / 11513 * pow(xi / Pi, 2.0) + 5670.0 / 11513 * pow(xi / Pi, 4.0) + 420.0 / 11513 * pow(xi / Pi, 6.0) + 9.0 / 11513 * pow(xi / Pi, 8.0);
    };
    auto nonequilibrium_diff_func_durca = [](double xi)
    {
        using constants::scientific::Pi;
        return 1 / Pi * (714.0 / 457 * pow(xi / Pi, 1.0) + 420.0 / 457 * pow(xi / Pi, 3.0) + 42.0 / 457 * pow(xi / Pi, 5.0));
    };
    auto nonequilibrium_diff_func_murca = [](double xi)
    {
        using constants::scientific::Pi;
        return 1 / Pi * (14680.0 / 11513 * pow(xi / Pi, 1.0) + 7560.0 / 11513 * pow(xi / Pi, 3.0) + 840.0 / 11513 * pow(xi / Pi, 5.0) + 24.0 / 11513 * pow(xi / Pi, 7.0));
    };

    // Cooling (is comprised of evolution of T and eta_e)

    // photon luminosity
    auto photon_luminosity = cooling::predefined::photonic::surface_luminosity(r_ns, m_ns, crust_eta);

    // neutrino luminosity
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
        if (it->first.classify() == auxiliaries::Species::ParticleClassification::kQuark)
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

    // pure rotochemical heating rate with (incorrect) accounting for superfluid suppression
    auto rotochemical_heating_e = [&](double r, double t, double T, double eta_e)
    {
        using namespace constants::scientific;
        double ratediff_durca = (hadron_durca_emissivity(r, electron, t, T) + quark_durca_emissivity(r, t, T)) / T * nonequilibrium_diff_func_durca(eta_e / T),
               ratediff_murca = (hadron_murca_emissivity(r, electron, t, T) + quark_murca_emissivity(r, t, T)) / T * nonequilibrium_diff_func_murca(eta_e / T);
        return eta_e * (ratediff_durca + ratediff_murca);
    };

    auto Q_nu = [&](double r, double t, double T, double eta_e)
    {
        using namespace constants::scientific;
        double result = 0;

        // hadronic part with non-equilibrium corrections
        for (auto it = m_stars_of_nbar.begin(); it != m_stars_of_nbar.end(); ++it)
        {
            auto key = it->first;
            if (key.classify() == auxiliaries::Species::ParticleClassification::kLepton)
            {
                if (key == electron)
                {
                    result += hadron_murca_emissivity(r, key, t, T) * enhance_murca(eta_e / T);
                    result += hadron_durca_emissivity(r, key, t, T) * enhance_durca(eta_e / T);
                }
                else
                {
                    result += hadron_murca_emissivity(r, key, t, T);
                    result += hadron_durca_emissivity(r, key, t, T);
                }
            }
            if (key.classify() == auxiliaries::Species::ParticleClassification::kBaryon)
            {
                result += hadron_PBF_emissivity(r, key, t, T);
            }
        }
        result += hadron_bremsstrahlung_emissivity(r, t, T);
        // quark part
        result += quark_durca_emissivity(r, t, T) + quark_murca_emissivity(r, t, T) + quark_bremsstrahlung_emissivity(r, t, T);
        result += electron_bremsstrahlung_emissivity(r, t, T);
        // pure rotochemical heating (must be subtracted)
        result -= rotochemical_heating_e(r, t, T, eta_e);

        return result * exp_phi(r) * exp_phi(r);
    };

    auto neutrino_luminosity = auxiliaries::integrate_volume<double, double, double>(
        std::function<double(double, double, double, double)>(Q_nu), 0, r_ns, exp_lambda, auxiliaries::IntegrationMode::kGaussLegendre_12p, radius_step);

    // specific heat
    auto fermi_specific_heat_dens = cooling::predefined::auxiliary::fermi_specific_heat_density(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_core_limit, exp_phi, superfluid_n_1s0,
        superfluid_p_1s0, superfluid_n_3p2, superfluid_p_temp, superfluid_n_temp, superconduct_q_gap);

    auto heat_capacity = auxiliaries::integrate_volume<double, double>(
        std::function<double(double, double, double)>(fermi_specific_heat_dens), 0, r_ns, exp_lambda, auxiliaries::IntegrationMode::kGaussLegendre_12p, radius_step);

    auto cooling_rhs = [&heat_capacity, &photon_luminosity, &neutrino_luminosity](double t, double T, double eta_e)
    {
        // std::cout << t << " " << photon_luminosity(t, T) << " " << neutrino_luminosity(t, T) << " " << heat_capacity(t, T) << '\n';
        return -(photon_luminosity(t + t_init, T) + neutrino_luminosity(t + t_init, T, eta_e)) / heat_capacity(t + t_init, T);
    };

    // we also need rhs for eta_e evolution
    auto eta_e_evo_rhs = [&](double t, double T, double eta_e)
    {
        using namespace constants::scientific;
        auto ratediff_local = [&](double r)
        {
            return ((hadron_durca_emissivity(r, electron, t, T) + quark_durca_emissivity(r, t, T)) / T * nonequilibrium_diff_func_durca(eta_e / T) +
                    (hadron_murca_emissivity(r, electron, t, T) + quark_murca_emissivity(r, t, T)) / T * nonequilibrium_diff_func_murca(eta_e / T)) *
                   pow(exp_phi(r), 2.0);
        };
        return -z_npe * (auxiliaries::integrate_volume<>(std::function<double(double)>(ratediff_local), 0, r_crust, exp_lambda, auxiliaries::IntegrationMode::kGaussLegendre_12p)() - 2 * omega_omega_dot(t) * I_e);
    };

    // solve cooling equation

    // we of course need a new solver to capture simultaneous evolution of T and eta_e

    using interpolator_t = std::function<double(const std::vector<double> &, const std::vector<double> &, double)>;

    auto cooling_solver = auxiliaries::CachedFunc<std::vector<std::vector<double>>, std::vector<double>, double, const std::function<double(double, double, double)> &, const std::function<double(double, double, double)> &, double, double, double, double,
                                                  const interpolator_t &, const interpolator_t &>(
        [&](std::vector<std::vector<double>> &cache, double t, const std::function<double(double, double, double)> &rhs_T, const std::function<double(double, double, double)> &rhs_eta_e, double T_init, double eta_e_init, double base_time_step, double exp_rate, const interpolator_t &interpolator_T, const interpolator_t &interpolator_eta_e)
        {
            // times must be positive
            if (t < 0 || base_time_step <= 0)
                throw std::invalid_argument("Evolution time must be positive");
            // inverse euler solver; I want this method to be stable for any time step, including huge ones
            auto back_euler_step = [&rhs_T, &rhs_eta_e](double t, double T, double eta_e, double time_step)
            {
                // solve T_{n+1} - T_n - dt * a^n * F(t_{n+1}, T_{n+1}, eta_{n+1}) = 0 ;
                // eta_{n+1} - eta_n - dt * a^n * G(t_{n+1}, T_{n+1}, eta_{n+1}) = 0 with Newton's steps
                double eps = 1e-5;
                size_t max_iter = 100, iter = 0;
                double T_new = T,
                       eta_new = eta_e;
                double F, G;
                do
                {
                    ++iter;
                    double f = rhs_T(t + time_step, T_new, eta_new),
                           g = rhs_eta_e(t + time_step, T_new, eta_new);
                    auto temp_step = time_step / (t + time_step) * T_new,
                         eta_step = time_step / (t + time_step) * T_new / 10000;
                    F = (T_new - T - time_step * f);
                    G = (eta_new - eta_e - time_step * g);
                    double dF_over_dT = 1 - time_step * (rhs_T(t + time_step, T_new + temp_step, eta_new) - f) / temp_step,
                           dG_over_dT = -time_step * (rhs_eta_e(t + time_step, T_new + temp_step, eta_new) - g) / temp_step,
                           dF_over_deta = -time_step * (rhs_T(t + time_step, T_new, eta_new + eta_step) - f) / eta_step,
                           dG_over_deta = 1 - time_step * (rhs_eta_e(t + time_step, T_new, eta_new + eta_step) - g) / eta_step;
                    T_new -= (dG_over_deta * F - dF_over_deta * G) / (dF_over_dT * dG_over_deta - dF_over_deta * dG_over_dT);
                    eta_new -= (dF_over_dT * G - dG_over_dT * F) / (dF_over_dT * dG_over_deta - dF_over_deta * dG_over_dT);
                    if (T_new < 0)
                        throw std::runtime_error("Reached negative temperature with current method; Encountered in stationary_cooling_cached");
                    if (iter > max_iter)
                        break;
                } while (std::abs(F) > eps * T);
                return std::vector<double>({T_new, eta_new});
            };

            // if we conduct a first calculation (or time exceeded existing cache), we need to fill the cache
            if (cache.empty() || cache[0].back() < t)
            {
                // create two empty vectors
                cache = std::vector<std::vector<double>>(3, std::vector<double>());
                // fill the cache[0] with time exp. growing values, and the cache[1] with the corresponding temperatures
                cache[0].push_back(0.0);
                cache[1].push_back(T_init);
                cache[2].push_back(eta_e_init);
                double t_curr = 0.0, time_step = base_time_step;
                do
                {
                    auto back_euler_step_result = back_euler_step(t_curr, cache[1].back(), cache[2].back(), time_step);
                    cache[1].push_back(back_euler_step_result[0]);
                    cache[2].push_back(back_euler_step_result[1]);
                    t_curr += time_step;
                    cache[0].push_back(t_curr);
                    time_step *= exp_rate;
                    // std::cout << cache[0].back() << " " << cache[1].back() << '\n';
                } while (t > t_curr);
            }
            // now we're sure the time is in the cache, we just interpolate
            return std::vector<double>({interpolator_T(cache[0], cache[1], t), interpolator_eta_e(cache[0], cache[2], t)});
        });

    double exp_phi_at_R = pow(1 - 2 * constants::scientific::G * m_ns / r_ns, 0.5);

    double T_init = T_init_local * exp_phi_at_R;

    auto eta_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                           double, const std::vector<double> &, const std::vector<double> &,
                                                           auxiliaries::InterpolationMode, double, bool, bool>(auxiliaries::interpolate_cached);
    auto eta_interpolator = [&eta_interpolator_cached](const std::vector<double> &x, const std::vector<double> &y, double val)
    {
        return eta_interpolator_cached(x, y, auxiliaries::InterpolationMode::kCubic, val, false, true);
    };

    // invoke the solver once to cache the solution
    cooling_solver(t_end - t_init, cooling_rhs, eta_e_evo_rhs, T_init, 0.0, base_t_step, exp_rate_estim, cooling_interpolator, eta_interpolator);

    // plot the solution (assumes exp_rate_estim > 1)
    std::vector<double> x(cooling_n_points_estimate, 0);
    std::vector<double> y(cooling_n_points_estimate, 0);
    std::cout << "M/Msol " << m_ns * constants::conversion::gev_over_msol << std::endl;
    std::cout << "t [years] "
              << "\tTe^inf [K] "
              << "\tL_ph [erg/s] "
              << "\tL_nu [erg/s] " << std::endl;
    for (int i = 0;; ++i)
    {
        x[i] = base_t_step * (pow(exp_rate_estim, i + 1) - 1) / (exp_rate_estim - 1);
        if (x[i] > t_end)
        {
            x.resize(i);
            y.resize(i);
            break;
        }
        auto res = cooling_solver(x[i], cooling_rhs, eta_e_evo_rhs, T_init, 0.0, base_t_step, exp_rate_estim, cooling_interpolator, eta_interpolator);
        y[i] = res[0];
        // print in understandable units
        std::cout << 1.0E6 * (x[i] + t_init) / (constants::conversion::myr_over_s * constants::conversion::gev_s) << "\t" << cooling::predefined::auxiliary::te_tb_relation(y[i], r_ns, m_ns, crust_eta) * exp_phi_at_R * constants::conversion::gev_over_k << "\t" << photon_luminosity(x[i], y[i]) * constants::conversion::gev_s / constants::conversion::erg_over_gev << "\t" << neutrino_luminosity(x[i], y[i], res[1]) * constants::conversion::gev_s / constants::conversion::erg_over_gev << "\t" << '\n';
        // rescale
        x[i] += t_init;
        x[i] *= 1.0E6 / (constants::conversion::myr_over_s * constants::conversion::gev_s);
        y[i] = cooling::predefined::auxiliary::te_tb_relation(y[i], r_ns, m_ns, crust_eta) * exp_phi_at_R * constants::conversion::gev_over_k;
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
    gr->Draw("AL");
    // title offset
    gr->GetYaxis()->SetTitleOffset(1.5);
    gPad->SetLogx();
    gPad->SetLogy();

    gr->GetXaxis()->SetTitle("t [yr]");
    gr->GetYaxis()->SetTitle("T [K]");

    auto legend = new TLegend(0.1, 0.1, 0.38, 0.38);
    legend->AddEntry(gr, "RH Manager", "l");

    legend->Draw();

    c1->SaveAs(pdf_path.c_str());
}