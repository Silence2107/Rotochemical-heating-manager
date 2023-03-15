
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
    auto hadron_durca_luminosity_electron = auxiliaries::CachedFunc<std::vector<double>, std::function<double(double, double)>,
                                                                    const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, double, double>(cooling::predefined::neutrinic::hadron_durca_luminocity_cached);
    auto hadron_durca_luminosity_muon = auxiliaries::CachedFunc<std::vector<double>, std::function<double(double, double)>,
                                                                const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, double, double>(cooling::predefined::neutrinic::hadron_durca_luminocity_cached);

    auto hadron_murca_luminosity_electron = auxiliaries::CachedFunc<std::vector<double>, std::function<double(double, double)>,
                                                                   const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, double, double, double>(cooling::predefined::neutrinic::hadron_murca_luminocity_cached);
    auto hadron_murca_luminosity_muon = auxiliaries::CachedFunc<std::vector<double>, std::function<double(double, double)>,
                                                                const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, double, double, double>(cooling::predefined::neutrinic::hadron_murca_luminocity_cached);

    auto hadron_bremsstrahlung_luminosity = auxiliaries::CachedFunc<std::vector<double>, std::function<double(double, double)>,
                                                                 const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, double, double>(cooling::predefined::neutrinic::hadron_bremsstrahlung_luminocity_cached);
    auto neutrino_luminosity = [&](double t, double T)
    {
        return hadron_durca_luminosity_electron(
                   m_stars_of_nbar.at("neutron"), m_stars_of_nbar.at("proton"), m_stars_of_nbar.at("electron"), k_fermi_of_nbar.at("neutron"), k_fermi_of_nbar.at("proton"), k_fermi_of_nbar.at("electron"), nbar, exp_lambda, exp_phi, r_ns, radius_step)(t, T) +
               hadron_durca_luminosity_muon(
                   m_stars_of_nbar.at("neutron"), m_stars_of_nbar.at("proton"), m_stars_of_nbar.at("muon"), k_fermi_of_nbar.at("neutron"), k_fermi_of_nbar.at("proton"), k_fermi_of_nbar.at("muon"), nbar, exp_lambda, exp_phi, r_ns, radius_step)(t, T) +
               hadron_murca_luminosity_electron(
                   m_stars_of_nbar.at("neutron"), m_stars_of_nbar.at("proton"), m_stars_of_nbar.at("electron"), k_fermi_of_nbar.at("neutron"), k_fermi_of_nbar.at("proton"), k_fermi_of_nbar.at("electron"), nbar, exp_lambda, exp_phi, r_ns, radius_step, nbar_conversion)(t, T) +
               hadron_murca_luminosity_muon(
                   m_stars_of_nbar.at("neutron"), m_stars_of_nbar.at("proton"), m_stars_of_nbar.at("muon"), k_fermi_of_nbar.at("neutron"), k_fermi_of_nbar.at("proton"), k_fermi_of_nbar.at("muon"), nbar, exp_lambda, exp_phi, r_ns, radius_step, nbar_conversion)(t, T) +
               hadron_bremsstrahlung_luminosity(
                   m_stars_of_nbar.at("neutron"), m_stars_of_nbar.at("proton"), k_fermi_of_nbar.at("neutron"), k_fermi_of_nbar.at("proton"), nbar, exp_lambda, exp_phi, r_ns, radius_step)(t, T);
    };

    // specific heat
    auto fermi_specific_heat = auxiliaries::CachedFunc<std::vector<double>, std::function<double(double, double)>,
                                                       const std::map<std::string, std::function<double(double)>> &, const std::map<std::string, std::function<double(double)>> &, const std::function<double(double)> &, const std::function<double(double)> &, const std::function<double(double)> &, double, double>(cooling::predefined::auxiliary::fermi_specific_heat_cached);

    auto heat_capacity = fermi_specific_heat(
        m_stars_of_nbar, k_fermi_of_nbar, nbar, exp_lambda, exp_phi, r_ns, radius_step);

    auto cooling_rhs = [&heat_capacity, &photon_luminosity, &neutrino_luminosity](double t, double T)
    {
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
    for (int i = 0; i <= 1000; ++i)
    {
        x[i] = base_t_step * (pow(exp_rate_estim, i+1) - 1) / (exp_rate_estim - 1);
        y[i] = cooling_solver(x[i], cooling_rhs, T_init, base_t_step, exp_rate_estim, cooling_interpolator);
        // print luminosities in humanic units
        std::cout << 1.0E6 * x[i] / (constants::conversion::myr_over_s * constants::conversion::gev_s) << " " << cooling::predefined::auxiliary::te_tb_relation(y[i], r_ns, m_ns, eta) * exp_phi_at_R * constants::conversion::gev_over_k << " " << 
            photon_luminosity(x[i], y[i]) * constants::conversion::gev_s / constants::conversion::erg_over_gev << " " << neutrino_luminosity(x[i], y[i]) * constants::conversion::gev_s / constants::conversion::erg_over_gev << std::endl;
        //rescale 
        x[i] *= 1.0E6 / (constants::conversion::myr_over_s * constants::conversion::gev_s);
        y[i] = cooling::predefined::auxiliary::te_tb_relation(y[i], r_ns, m_ns, eta) * exp_phi_at_R * constants::conversion::gev_over_k;
    }

    // compare with nscool
    std::ifstream apr_2_1_nscool("../../data/Teff_Try.dat");
    std::vector<double> x_nscool, y_nscool;
    // iterate over file, but skip 25 lines
    for (int i = 0; i < 25; ++i)
    {
        std::string line;
        std::getline(apr_2_1_nscool, line);
    }
    while (apr_2_1_nscool.good())
    {
        std::string line;
        std::getline(apr_2_1_nscool, line);
        line = auxiliaries::retrieve_cleared_line(line);
        std::stringstream ss(line);
        double step, t, T;
        ss >> step >> t >> T;
        x_nscool.push_back(t);
        y_nscool.push_back(T);
    }

    TCanvas *c1 = new TCanvas("c1", "c1");
    auto gr = new TGraph(1000, x.data(), y.data());
    gr->SetLineColor(kBlue);
    gr->Draw("AL");
    // title offset
    gr->GetYaxis()->SetTitleOffset(1.5);
    gPad->SetLogx();

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