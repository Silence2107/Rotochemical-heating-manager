
#include "../include/auxiliaries.h"
#include "../include/cooling.h"
#include "../include/eos_reader.h"
#include "../include/constants.h"
#include "../include/tov_solver.h"
#include "../include/inputfile.hpp"

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

int main()
{
    using namespace inputfile;

    // Cooling solver setup
    auto cooling_interpolator_cached = auxiliaries::CachedFunc<std::function<double(double)>,
                                                               double, const std::vector<double> &, const std::vector<double> &,
                                                               auxiliaries::InterpolationMode, double, bool>(auxiliaries::interpolate_cached);
    auto cooling_interpolator = [&cooling_interpolator_cached](const std::vector<double> &x, const std::vector<double> &y, double val)
    {
        return cooling_interpolator_cached(x, y, auxiliaries::InterpolationMode::kCubic, val, false);
    };

    // RUN --------------------------------------------------------------------------

    // EoS definition

    auto eos_cached = auxiliaries::CachedFunc<std::vector<std::vector<double>>, double, double>(
        [&](std::vector<std::vector<double>> &cache, double rho)
        {
            if (rho < 0 || rho > edensity_upp * energy_density_conversion)
                throw std::runtime_error("Data request out of range; Encountered in main::eos_cached");
            if (rho <= edensity_low * energy_density_conversion)
                return 0.0;
            if (cache.empty() || cache[0].size() != discr_size_EoS)
            {                                                                                        // then fill/refill cache
                cache = std::vector<std::vector<double>>(2, std::vector<double>(discr_size_EoS, 0)); // initialize 2xdiscr_size_EoS matrix
                std::vector<double> x(discr_size_EoS, 0);
                for (int i = 1; i < discr_size_EoS - 1; ++i)
                { // cache EoS for further efficiency
                    x[i] = i * (nbar_upp - nbar_low) / discr_size_EoS + nbar_low;
                    cache[0][i] = energy_density_conversion * energy_density_of_nbar(x[i]);
                    cache[1][i] = pressure_conversion * pressure_of_nbar(x[i]);
                }
                x[0] = nbar_low;
                x[x.size() - 1] = nbar_upp;
                cache[0][0] = energy_density_conversion * edensity_low;
                cache[0][cache[0].size() - 1] = energy_density_conversion * edensity_upp;
                cache[1][0] = pressure_conversion * pressure_low;
                cache[1][cache[1].size() - 1] = pressure_conversion * pressure_upp;
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

    auto nbar = auxiliaries::CachedFunc<std::vector<std::vector<double>>, double, double>(
        [&](std::vector<std::vector<double>> &cache, double r)
        {
            // cache contains {r, n_B(r)} arrays; recaching is not supported at the moment, call ::erase instead
            // return nbar(r) for given r (in datafile units)

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

                    double nbar_left = nbar_low, nbar_right = nbar_upp; // we need these for bisection search; in fm-3 units for now
                    double nbar_mid = (nbar_left + nbar_right) / 2.0;
                    while (fabs(nbar_right - nbar_left) > nbar_low)
                    {
                        // while we are too far from appropriate precision for nbar estimate
                        // recalculate via bisection method
                        nbar_mid = (nbar_left + nbar_right) / 2.0;
                        if (energy_density_conversion * energy_density_of_nbar(nbar_mid) > density_at_r)
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

    // ready to run

    double r_ns = tov(0.0)[4];
    double m_ns = tov(r_ns)[0];
    double eta = 2.26E-18;

    auto exp_phi = [&tov](double r)
    {
        return std::exp(tov(r)[2]);
    };

    auto exp_lambda = [&tov](double r)
    {
        return pow(1 - 2 * constants::scientific::G * tov(r)[0] / r, -0.5);
    };

    // photon luminosity
    auto photon_luminosity = cooling::predefined::photonic::surface_luminosity(r_ns, m_ns, eta);

    // neutrino luminosity
    auto hadron_durca_emissivity = [&](double r, const std::string& lepton_flavour, double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::conversion;

        double nbar_val = nbar(r);
        double pf_l = k_fermi_of_nbar.at(lepton_flavour)(nbar_val),
                pf_n = k_fermi_of_nbar.at("neutron")(nbar_val),
                pf_p = k_fermi_of_nbar.at("proton")(nbar_val);
        double mst_n = m_stars_of_nbar.at("neutron")(nbar_val),
                mst_p = m_stars_of_nbar.at("proton")(nbar_val),
                mst_l = m_stars_of_nbar.at(lepton_flavour)(nbar_val);
        if (pf_l + pf_p - pf_n <= 0)
            return 0.0;
        double dens = (4.24E27 / 1.68E54) * (mst_n / M_N) * (mst_p / M_N) * mst_l *
                        pow(T/exp_phi(r), 6) * pow(gev_over_k, 6) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);
        return dens;
    };

    auto hadron_murca_emissivity = [&](double r, const std::string& lepton_flavour, double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::conversion;

        double nbar_val = nbar(r);
        double pf_l = k_fermi_of_nbar.at(lepton_flavour)(nbar_val),
                pf_n = k_fermi_of_nbar.at("neutron")(nbar_val),
                pf_p = k_fermi_of_nbar.at("proton")(nbar_val);
        double mst_n = m_stars_of_nbar.at("neutron")(nbar_val),
                mst_p = m_stars_of_nbar.at("proton")(nbar_val),
                mst_l = m_stars_of_nbar.at(lepton_flavour)(nbar_val);
        double alpha = 1.76 - 0.63 * pow(N_sat / (nbar_val * nbar_conversion), 2.0 / 3), beta = 0.68,
               v_fl = pf_l / mst_l;
        double dens = (8.05E21 / 1.68E72) * v_fl * pow(mst_n / M_N, 3) * (mst_p / M_N) * pf_p *
                      pow(T/exp_phi(r), 8) * alpha * beta *
                      pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);
        if (pf_l + 3 * pf_p - pf_n > 0)
            dens += (8.05E21 / (8 * 1.68E72)) * (pow(pf_l + 3 * pf_p - pf_n, 2) / mst_l) * pow(mst_p / M_N, 3) * (mst_n / M_N) *
                    pow(T/exp_phi(r), 8) * alpha * beta *
                    pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);
        return dens;
    };

    auto hadron_bremsstrahlung_emissivity = [&](double r, double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::conversion;

        double nbar_val = nbar(r);
        double pf_n = k_fermi_of_nbar.at("neutron")(nbar_val),
                pf_p = k_fermi_of_nbar.at("proton")(nbar_val);
        double mst_n = m_stars_of_nbar.at("neutron")(nbar_val),
                mst_p = m_stars_of_nbar.at("proton")(nbar_val);
        double alpha_nn = 0.59, alpha_np = 1.06, alpha_pp = 0.11,
               beta_nn = 0.56, beta_np = 0.66, beta_pp = 0.7;
        int n_flavours = 3;
        double dens_nn = (7.5E19 / 1.68E72) * pow(mst_n / M_N, 4) * pf_n * n_flavours *
                         pow(T/exp_phi(r), 8) * alpha_nn * beta_nn *
                         pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3),
               dens_pp = (7.5E19 / 1.68E72) * pow(mst_p / M_N, 4) * pf_p * n_flavours *
                         pow(T/exp_phi(r), 8) * alpha_pp * beta_pp *
                         pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3),
               dens_np = (1.5E20 / 1.68E72) * pow(mst_p / M_N, 2) * pow(mst_n / M_N, 2) * pf_p * n_flavours *
                         pow(T/exp_phi(r), 8) * alpha_np * beta_np *
                         pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);

        return dens_nn + dens_np + dens_pp;
    };

    auto Q_nu = [&](double r, double t, double T)
    {
        double result = 0;
        for (auto it = m_stars_of_nbar.begin(); it != m_stars_of_nbar.end(); ++it)
        {
            auto key = it->first;
            if (key == "electron" || key == "muon" || key == "tau")
            {
                result += hadron_murca_emissivity(r, key, t, T);
                result += hadron_durca_emissivity(r, key, t, T);
            }
        }
        result += hadron_bremsstrahlung_emissivity(r, t, T);
        return result * exp_phi(r) * exp_phi(r);
    };

    auto neutrino_luminosity = auxiliaries::integrate_volume<double, double>(
        std::function<double(double, double, double)>(Q_nu), 0, r_ns, radius_step, exp_lambda);

    // specific heat
    auto fermi_specific_heat_dens = [&](double r, double t, double T)
    {
        using namespace constants::scientific;

        double cv_dens = 0;
        for (auto it = m_stars_of_nbar.begin(); it != m_stars_of_nbar.end(); ++it)
        {
            auto key = it->first;
            double nbar_val = nbar(r);
            double m_star = m_stars_of_nbar.at(key)(nbar_val);
            double k_fermi = k_fermi_of_nbar.at(key)(nbar_val);
            double exp_min_phi = 1.0 / exp_phi(r);
            cv_dens += m_star * k_fermi / 3.0 * exp_min_phi * T;
        }
        return cv_dens;
    };

    auto heat_capacity = auxiliaries::integrate_volume<double, double>(
        std::function<double(double, double, double)>(fermi_specific_heat_dens), 0, r_ns, radius_step, exp_lambda);

    auto cooling_rhs = [&heat_capacity, &photon_luminosity, &neutrino_luminosity](double t, double T)
    {
        //std::cout << t << " " << photon_luminosity(t, T) << " " << neutrino_luminosity(t, T) << " " << heat_capacity(t, T) << '\n';
        return -(photon_luminosity(t, T) + neutrino_luminosity(t, T)) / heat_capacity(t, T);
    };

    // solve cooling equation

    auto cooling_solver = auxiliaries::CachedFunc<std::vector<std::vector<double>>, double, double, const std::function<double(double, double)> &, double, double, double,
                                                  const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &>(cooling::solver::stationary_cooling_cached);

    double exp_phi_at_R = pow(1 - 2 * constants::scientific::G * m_ns / r_ns, 0.5);

    // evolution options
    double t_end = 6.1E-0 * constants::conversion::myr_over_s * constants::conversion::gev_s,
           T_init = 5E9 * exp_phi_at_R / constants::conversion::gev_over_k,
           base_t_step = 1.0E-18 * constants::conversion::myr_over_s * constants::conversion::gev_s;

    double exp_rate_estim = pow(t_end / base_t_step, 1.0 / 1000) * pow((pow(t_end / base_t_step, 1.0 / 1000) - 1), 1.0 / 1000);

    // invoke the solver once to cache the solution
    cooling_solver(t_end, cooling_rhs, T_init, base_t_step, exp_rate_estim, cooling_interpolator);

    // plot the solution
    std::vector<double> x(1000, 0);
    std::vector<double> y(1000, 0);
    std::cout << "M/Msol " << m_ns * constants::conversion::gev_over_msol << std::endl;
    std::cout << "t [years] " << "\tTe^inf [K] " << "\tL_ph [erg/s] " << "\tL_nu [erg/s] " << std::endl;
    for (int i = 0;; ++i)
    {
        x[i] = base_t_step * (pow(exp_rate_estim, i+1) - 1) / (exp_rate_estim - 1);
        if (x[i] > t_end)
        {
            x.resize(i);
            y.resize(i);
            break;
        }
        y[i] = cooling_solver(x[i], cooling_rhs, T_init, base_t_step, exp_rate_estim, cooling_interpolator);
        // print luminosities in humanic units
        std::cout << 1.0E6 * x[i] / (constants::conversion::myr_over_s * constants::conversion::gev_s) << "\t" << cooling::predefined::auxiliary::te_tb_relation(y[i], r_ns, m_ns, eta) * exp_phi_at_R * constants::conversion::gev_over_k << "\t" << 
            photon_luminosity(x[i], y[i]) * constants::conversion::gev_s / constants::conversion::erg_over_gev << "\t" << neutrino_luminosity(x[i], y[i]) * constants::conversion::gev_s / constants::conversion::erg_over_gev << "\t" <<
            heat_capacity(x[i], y[i]) << '\n';
        //rescale 
        x[i] *= 1.0E6 / (constants::conversion::myr_over_s * constants::conversion::gev_s);
        y[i] = cooling::predefined::auxiliary::te_tb_relation(y[i], r_ns, m_ns, eta) * exp_phi_at_R * constants::conversion::gev_over_k;
    }

    // compare with nscool
    std::ifstream apr_nscool("../../data/1.4.dat");
    std::vector<double> x_nscool, y_nscool;
    // iterate over file, but skip 25 lines
    for (int i = 0; i < 25; ++i)
    {
        std::string line;
        std::getline(apr_nscool, line);
    }
    while (apr_nscool.good())
    {
        std::string line;
        std::getline(apr_nscool, line);
        line = auxiliaries::retrieve_cleared_line(line);
        std::stringstream ss(line);
        double step, t, T;
        ss >> step >> t >> T;
        x_nscool.push_back(t);
        y_nscool.push_back(T);
    }

    TCanvas *c1 = new TCanvas("c1", "c1");
    auto gr = new TGraph(x.size(), x.data(), y.data());
    gr->SetLineColor(kBlue);
    gr->Draw("AL");
    // title offset
    gr->GetYaxis()->SetTitleOffset(1.5);
    gPad->SetLogx();
    gPad->SetLogy();

    gr->GetXaxis()->SetTitle("t [yr]");
    gr->GetYaxis()->SetTitle("T [K]");

    auto gr_ns_cool = new TGraph(x_nscool.size(), x_nscool.data(), y_nscool.data());
    gr_ns_cool->SetLineColor(kRed);
    gr_ns_cool->Draw("L");

    auto legend = new TLegend(0.1, 0.1, 0.38, 0.38);
    legend->AddEntry(gr, "RH Manager", "l");
    legend->AddEntry(gr_ns_cool, "NSCool", "l");

    legend->Draw();

    c1->SaveAs("cooling.pdf");

    // plot nbar of r in the whole star
    /*
    std::vector<double> x(1000, 0);
    std::vector<double> y(1000, 0);
    for (int i = 0; i < 1000; ++i)
    {
        x[i] = i * r_ns / 1000.0;
        y[i] = nbar(x[i]);
        x[i] /= constants::conversion::km_gev;
    }

    TCanvas *c1 = new TCanvas("c1", "c1");
    auto gr = new TGraph(1000, x.data(), y.data());
    gr->Draw("AL");
    // title offset
    gr->GetYaxis()->SetTitleOffset(1.5);

    gr->GetXaxis()->SetTitle("r [km]");
    gr->GetYaxis()->SetTitle("n_B [fm^{-3}]");
    c1->SaveAs("nbar.pdf");
    */

    // plot photon luminosity for temperatures between 1 and 100 MeV

    /*
    std::vector<double> x(100, 0);
    std::vector<double> y(100, 0);
    for (int i = 0; i < 100; ++i)
    {
        x[i] = 1 + i * 0.99;
        y[i] = photon_luminosity(0, x[i] / 1000.0) * constants::conversion::gev_s / constants::conversion::erg_over_gev;
    }

    TCanvas *c1 = new TCanvas("c1", "c1");
    auto gr = new TGraph(100, x.data(), y.data());
    gr->Draw("AL");
    // title offset
    gr->GetYaxis()->SetTitleOffset(1.5);

    gr->GetXaxis()->SetTitle("T [MeV]");
    gr->GetYaxis()->SetTitle("L [erg/s]");
    c1->SaveAs("photon_luminosity.pdf");
    */
}