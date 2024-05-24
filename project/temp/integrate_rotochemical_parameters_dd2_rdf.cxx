
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
#include <algorithm>

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
    argparse::ArgumentParser parser("integrate_rotochemical_parameters_dd2_rdf", "Print I, Z, W quantities in memorized order based on DD2-RDF EoS", "Argparse powered by SiLeader");

#if RHM_REQUIRES_INPUTFILE
    parser.addArgument({"--inputfile"}, "json input file path (required)");
#endif
#if RHM_HAS_ROOT
    parser.addArgument({"--pdf_path"}, "pdf output file path (optional, default: Cooling.pdf)");
    parser.addArgument({"--rootfile_path"}, "root output file path (optional, default: None)");
#endif
    parser.addArgument({"--center_pressure"}, "override for center pressure fraction (optional, default: fetched from inputfile)");

    auto args = parser.parseArgs(argc, argv);

    using namespace instantiator;
#if RHM_REQUIRES_INPUTFILE
    instantiator::instantiate_system(args.get<std::string>("inputfile"));
#endif

#if RHM_HAS_ROOT
    std::string pdf_path = args.safeGet<std::string>("pdf_path", "Cooling.pdf");
    TFile *rootfile = nullptr;
    if (args.has("rootfile_path"))
        rootfile = new TFile(args.get<std::string>("rootfile_path").c_str(), "RECREATE");
#endif

    double center_pressure;
    if (!args.has("center_pressure"))
        center_pressure = instantiator::center_pressure;
    else
        center_pressure = args.get<double>("center_pressure") * (instantiator::pressure_upp - instantiator::pressure_low) + instantiator::pressure_low;

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
                for (size_t i = 1; i < discr_size_EoS - 1; ++i)
                { // cache EoS for further efficiency
                    x[i] = i * (nbar_upp - nbar_low) / discr_size_EoS + nbar_low;
                    cache[0][i] = pressure_of_nbar(x[i]);
                    cache[1][i] = energy_density_of_nbar(x[i]);
                }
                x[0] = nbar_low;
                x[x.size() - 1] = nbar_upp;
                cache[0][0] = pressure_low;
                cache[0][cache[1].size() - 1] = pressure_upp;
                cache[1][0] = edensity_low;
                cache[1][cache[0].size() - 1] = edensity_upp;
                eos_interpolator_cached.erase(); // clean up cached interpolator
            }
            return eos_interpolator(cache[0], cache[1], p);
        });

    // TOV solver

    auto tov_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                                    const std::function<double(double)> &, double, double, double, double, size_t>(tov_solver::tov_solution);
    auto tov = [&tov_cached, &eos_inv_cached, center_pressure](double r)
    {
        // TOV solution cached
        return tov_cached(eos_inv_cached, r, center_pressure, radius_step, surface_pressure, tov_adapt_limit);
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

    // integrate the bij integrand over volume

    // instantiate bij; bee is only integrated over hadronic phase
    double b_ee = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                             { return dne_to_dmue(nbar(r)) / exp_phi(r) * (number_densities_of_nbar[constants::species::neutron] != 0); }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kRectangular)(),
           b_em = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                             { return dne_to_dmum(nbar(r)) / exp_phi(r); }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)(),
           b_mm = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                             { return dnm_to_dmum(nbar(r)) / exp_phi(r); }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)(),
           b_uu = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                             { return dnu_to_dmuu(nbar(r)) / exp_phi(r); }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)(),
           b_us = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                             { return dnu_to_dmus(nbar(r)) / exp_phi(r); }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)(),
           b_ss = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                             { return dns_to_dmus(nbar(r)) / exp_phi(r); }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)();
    // std::cout << "b_ee = " << b_ee << " b_em = " << b_em << " b_mm = " << b_mm << " b_uu = " << b_uu << " b_us = " << b_us << " b_ss = " << b_ss << std::endl;
    // We have explicit expressions for inverted bij, so let's calculate them
    double z_npe, z_npm, z_np, z_due, z_us, z_sue;
    /*
    z_npe = -b_mm / [b_em * b_em - b_ee * b_mm],
    z_np = b_em / [b_em * b_em - b_ee * b_mm],
    z_npm = -b_ee / [b_em * b_em - b_ee * b_mm],
    z_due = -b_ss / [b_us * b_us - b_uu * b_ss],
    z_us = b_us / [b_us * b_us - b_uu * b_ss],
    z_sue = (b_uu + b_us) / [b_us * b_us - b_uu * b_ss];
    */

    // partial pivoting would be nicer, but let's stick to comparing with zeros for now
    if (b_ee * b_mm - b_em * b_em == 0)
    {
        // electron-muon subsystem is singular
        if (b_em == 0)
        {
            // then there's definitely no cross term
            z_np = 0;
            if (b_ee == 0)
            {
                // then there's no electron imbalance
                z_npe = 0;
            }
            else
            {
                z_npe = 1 / b_ee;
            }

            if (b_mm == 0)
            {
                // then there's no muon imbalance
                z_npm = 0;
            }
            else
            {
                z_npm = 1 / b_mm;
            }
        }
        else
        {
            // exotic case when particle deviations are proportional
            RHM_THROW(std::runtime_error, "Bij matrix is singular beyond expectations. Probable cause: provided system features unexpected conservation law.");
        }
    }
    else
    {
        z_npe = -b_mm / (b_em * b_em - b_ee * b_mm);
        z_np = b_em / (b_em * b_em - b_ee * b_mm);
        z_npm = -b_ee / (b_em * b_em - b_ee * b_mm);
    }

    if (b_us * b_us - b_uu * b_ss == 0)
    {
        // uquark-squark subsystem is singular
        if (b_ss == 0)
        {
            // then b_us is zero as well
            // there's no squark imbalance
            z_us = 0;
            z_sue = 0;
            if (b_uu == 0)
            {
                // then there's no uquark imbalance
                z_due = 0;
            }
            else
            {
                z_due = 1 / b_uu;
            }
        }
        else if (b_uu == 0)
        {
            // then there's no uquark imbalance
            z_due = 0;
            z_us = 0;
            if (b_ss == 0)
            {
                // then there's no squark imbalance
                z_sue = 0;
            }
            else
            {
                z_sue = -1 / b_ss;
            }
        }
        else
        {
            // exotic case when particle deviations are proportional
            RHM_THROW(std::runtime_error, "Bij matrix is singular beyond expectations. Probable cause: provided system features unexpected conservation law.");
        }
    }
    else
    {
        z_due = -b_ss / (b_us * b_us - b_uu * b_ss);
        z_us = b_us / (b_us * b_us - b_uu * b_ss);
        z_sue = (b_uu + b_us) / (b_us * b_us - b_uu * b_ss);
    }

    // Now we can construct the Z matrix that transforms particle deviations to chemical imbalances.
    auxiliaries::math::MatrixD zij(4, 4, 0);

    zij.at(0, 0) = -z_npe;
    zij.at(0, 1) = -z_np;
    zij.at(1, 0) = -z_np;
    zij.at(1, 1) = -z_npm;
    zij.at(2, 2) = -z_due;
    zij.at(2, 3) = -z_us;
    zij.at(3, 2) = (z_us - z_due);
    zij.at(3, 3) = -z_sue;

    // Let's now clear empty rows, if any
    auto rh_particles = std::vector<auxiliaries::phys::Species>{
        constants::species::electron, constants::species::muon,
        constants::species::uquark, constants::species::squark};
    auto supported_rh_particles = rh_particles;

    for (size_t row = 0; row < zij.rows(); ++row)
    {
        bool empty_row = true;
        for (size_t col = 0; col < zij.columns(); ++col)
        {
            if (zij.at(row, col) != 0)
            {
                empty_row = false;
                break;
            }
        }
        if (empty_row)
        {
            auto zij_new = auxiliaries::math::MatrixD(zij.rows() - 1, zij.columns(), 0);
            for (size_t i = 0; i < zij_new.rows(); ++i)
            {
                for (size_t j = 0; j < zij_new.columns(); ++j)
                {
                    zij_new.at(i, j) = zij.at(i + (i >= row), j);
                }
            }
            zij = zij_new;
            rh_particles.erase(rh_particles.begin() + row);
            --row;
        }
    }

    // I_omega estimates

    std::map<auxiliaries::phys::Species, double> i_omegas;

    for (auto species = number_densities_of_nbar.begin(); species != number_densities_of_nbar.end(); ++species)
    {
        auto n_i = species->second;
        double omega_k_sqr = pow(2.0 / 3, 3.0) * constants::scientific::G * m_ns / (r_ns * r_ns * r_ns);
        auto integrand = [&](double r)
        {
            return -nbar(r) * 1.0 / omega_k_sqr * (n_i(nbar(r + radius_step)) / nbar(r + radius_step) - n_i(nbar(r)) / nbar(r)) / (tov(r + radius_step)[3] / tov(r)[3] - 1);
        };
        i_omegas[species->first] = auxiliaries::math::integrate_volume<>(std::function<double(double)>(integrand), 0.0, r_ns - radius_step, exp_lambda, auxiliaries::math::IntegrationMode::kRectangular)() / (constants::conversion::gev_s * constants::conversion::gev_s);
    }
    for (auto rh_species = supported_rh_particles.begin(); rh_species != supported_rh_particles.end(); ++rh_species)
    {
        if (i_omegas.find(*rh_species) == i_omegas.end())
        {
            i_omegas[*rh_species] = 0;
        }
    }

    size_t indent = 20;
    // std::cout << std::left << std::setw(indent) << "I_e" << std::setw(indent) << "I_n" << std::setw(indent) << "I_p" << std::setw(indent) << "I_u" << std::setw(indent) << "I_d" << std::setw(indent) << "Z_npe" << std::setw(indent) << "Z_due" << std::setw(indent) << "W_npe" << std::setw(indent) << "W_due" << std::endl;
    std::cout << std::left << std::setw(indent) << nbar(0) * constants::conversion::fm3_gev3 << std::setw(indent) << m_ns * constants::conversion::gev_over_msol << std::setw(indent) << i_omegas[constants::species::electron] << std::setw(indent) << i_omegas[constants::species::neutron] << std::setw(indent) << i_omegas[constants::species::proton] << std::setw(indent) << i_omegas[constants::species::uquark] << std::setw(indent) << i_omegas[constants::species::dquark] << std::setw(indent) << z_npe / constants::conversion::erg_over_gev << std::setw(indent) << z_due / constants::conversion::erg_over_gev << std::endl;
}