
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
    argparse::ArgumentParser parser("cooling_npemuds_rotochemical", "Solves the cooling equation coupled with chemical imbalances based on EoS", "Argparse powered by SiLeader");

#if RHM_REQUIRES_INPUTFILE
    parser.addArgument({"--inputfile"}, "json input file path (required)");
#endif
#if RHM_HAS_ROOT
    parser.addArgument({"--pdf_path"}, "pdf output file path (optional, default: Cooling-with-RH.pdf)");
    parser.addArgument({"--rootfile_path"}, "root output file path (optional, default: None)");
#endif
    parser.addArgument({"--save_chemical_imbalances"}, "print & write chemical imbalances (optional, value-free, default: false)", argparse::ArgumentType::StoreTrue);

    auto args = parser.parseArgs(argc, argv);

    using namespace instantiator;
#if RHM_REQUIRES_INPUTFILE
    instantiator::instantiate_system(args.get<std::string>("inputfile"), {"TOV", "COOL", "RH"});
#endif

#if RHM_HAS_ROOT
    std::string pdf_path = args.safeGet<std::string>("pdf_path", "Cooling-with-RH.pdf");
    TFile *rootfile = nullptr;
    if (args.has("rootfile_path"))
        rootfile = new TFile(args.get<std::string>("rootfile_path").c_str(), "RECREATE");
#endif

    bool save_chemical_imbalances = args.has("save_chemical_imbalances");

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

    // integrate the bij integrand over volume

    // instantiate bij; bee is only integrated over hadronic phase
    double b_ee = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                             {  
                    using namespace constants::species;
                    try
                    {
                        return dni_to_dmuj[{electron, electron}](nbar(r)) / exp_phi(r) * (number_densities_of_nbar[neutron](nbar(r)) != 0);
                    }
                    catch(const std::exception& e)
                    {
                        return 0.0;
                    } }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kRectangular)(),
           b_em = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                             {  
                    using namespace constants::species;
                    try
                    {
                        return dni_to_dmuj[{electron, muon}](nbar(r)) / exp_phi(r);
                    }
                    catch(const std::exception& e)
                    {
                        return 0.0;
                    } }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kRectangular)(),
           b_mm = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                                {  
                        using namespace constants::species;
                        try
                        {
                            return dni_to_dmuj[{muon, muon}](nbar(r)) / exp_phi(r);
                        }
                        catch(const std::exception& e)
                        {
                            return 0.0;
                        } }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kRectangular)(),
           b_uu = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                             {
                    using namespace constants::species;
                    try
                    {
                        return dni_to_dmuj[{uquark, uquark}](nbar(r)) / exp_phi(r);
                    }
                    catch(const std::exception& e)
                    {
                        return 0.0;
                    } }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kRectangular)(),
           b_us = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                             {
                    using namespace constants::species;
                    try
                    {
                        return dni_to_dmuj[{uquark, squark}](nbar(r)) / exp_phi(r);
                    }
                    catch(const std::exception& e)
                    {
                        return 0.0;
                    } }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kRectangular)(),
           b_ss = auxiliaries::math::integrate_volume<>(
               std::function<double(double)>([&nbar, &exp_phi](double r)
                                             { 
                    using namespace constants::species;
                    try
                    {
                        return dni_to_dmuj[{squark, squark}](nbar(r)) / exp_phi(r);
                    }
                    catch(const std::exception& e)
                    {
                        return 0.0;
                    } }),
               0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kRectangular)();
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
    zij.at(3, 2) = z_us - z_due;
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

    // cooling settings

    // internal emission
    auto hadron_durca_enhanced_emissivity = cooling::predefined::rotochemical::hadron_durca_enhanced_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp);

    auto hadron_murca_enhanced_emissivity = cooling::predefined::rotochemical::hadron_murca_enhanced_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp);

    auto hadron_bremsstrahlung_emissivity = cooling::predefined::neutrinic::hadron_bremsstrahlung_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, ion_volume_fr, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp);

    auto hadron_PBF_emissivity = cooling::predefined::neutrinic::hadron_pbf_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp);

    auto quark_ud_durca_enhanced_emissivity = cooling::predefined::rotochemical::quark_ud_durca_enhanced_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto quark_us_durca_enhanced_emissivity = cooling::predefined::rotochemical::quark_us_durca_enhanced_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto quark_ud_murca_enhanced_emissivity = cooling::predefined::rotochemical::quark_ud_murca_enhanced_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto quark_us_murca_enhanced_emissivity = cooling::predefined::rotochemical::quark_us_murca_enhanced_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto quark_bremsstrahlung_emissivity = cooling::predefined::neutrinic::quark_bremsstrahlung_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto electron_bremsstrahlung_emissivity = cooling::predefined::neutrinic::electron_bremsstrahlung_emissivity(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi);

    // rate differences

    auto hadron_durca_rate_difference = cooling::predefined::rotochemical::hadron_durca_rate_difference(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp);

    auto hadron_murca_rate_difference = cooling::predefined::rotochemical::hadron_murca_rate_difference(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp);

    auto quark_ud_durca_rate_difference = cooling::predefined::rotochemical::quark_ud_durca_rate_difference(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto quark_us_durca_rate_difference = cooling::predefined::rotochemical::quark_us_durca_rate_difference(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto quark_ud_murca_rate_difference = cooling::predefined::rotochemical::quark_ud_murca_rate_difference(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

    auto quark_us_murca_rate_difference = cooling::predefined::rotochemical::quark_us_murca_rate_difference(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, exp_phi, superconduct_q_gap);

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
        i_omegas[species->first] = auxiliaries::math::integrate_volume<>(std::function<double(double)>(integrand), 0.0, r_ns - radius_step, exp_lambda, auxiliaries::math::IntegrationMode::kRectangular)();
    }
    for (auto rh_species = constants::species::known_particles.begin(); rh_species != constants::species::known_particles.end(); ++rh_species)
    {
        if (i_omegas.find(*rh_species) == i_omegas.end())
        {
            i_omegas[*rh_species] = 0;
        }
    }

    auto Q_nu = [&](double r, double t, double T, const std::map<auxiliaries::phys::Species, double> &etas)
    {
        using namespace constants::scientific;
        using namespace constants::species;
        double result = 0;

        // hadron processes
        for (auto it = m_stars_of_nbar.begin(); it != m_stars_of_nbar.end(); ++it)
        {
            auto key = it->first;
            if (key.classify() == auxiliaries::phys::Species::ParticleClassification::kLepton)
            {
                // ignore taus
                result += hadron_murca_enhanced_emissivity(r, key, t, T, etas.at(key));
                result += hadron_durca_enhanced_emissivity(r, key, t, T, etas.at(key));
            }
            if (key.classify() == auxiliaries::phys::Species::ParticleClassification::kBaryon)
            {
                result += hadron_PBF_emissivity(r, key, t, T);
            }
        }
        result += hadron_bremsstrahlung_emissivity(r, t, T);
        // quark processes
        result += quark_ud_durca_enhanced_emissivity(r, t, T, etas.at(uquark)) +
                  quark_us_durca_enhanced_emissivity(r, t, T, etas.at(squark)) +
                  quark_ud_murca_enhanced_emissivity(r, t, T, etas.at(uquark)) +
                  quark_us_murca_enhanced_emissivity(r, t, T, etas.at(squark)) +
                  quark_bremsstrahlung_emissivity(r, t, T);
        // extra
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
    auto neutrino_luminosity = auxiliaries::math::integrate_volume<double, double, const std::map<auxiliaries::phys::Species, double> &>(
        std::function<double(double, double, double, const std::map<auxiliaries::phys::Species, double> &)>(Q_nu), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p, radius_step);

    auto heat_capacity = auxiliaries::math::integrate_volume<double, double>(
        std::function<double(double, double, double)>(fermi_specific_heat_dens), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p, radius_step);

    auto rotochemical_vector_rhs = [&](double t, const std::vector<double> &funcs)
    {
        using namespace constants::species;
        std::vector<double> evaluated_rhs;
        // create vector of etas with zero values for trivial etas; funcs[1+] are nontrivial etas
        std::map<auxiliaries::phys::Species, double> etas;
        // funcs[0] is temperature
        double T = funcs[0];

        for (size_t species_count = 0; species_count < supported_rh_particles.size(); ++species_count)
        {
            auto pos = std::find(rh_particles.begin(), rh_particles.end(), supported_rh_particles[species_count]);
            if (pos != rh_particles.end())
            {
                etas[supported_rh_particles[species_count]] = funcs[pos - rh_particles.begin() + 1];
            }
            else
            {
                etas[supported_rh_particles[species_count]] = 0;
            }
        }

        // rotochemical heating in npemuds
        auto diff_npl_inf = [&](double r, const auxiliaries::phys::Species &lepton)
        {
            return exp_phi(r) * (number_densities_of_nbar[constants::species::neutron](nbar(r)) != 0) *
                   (hadron_durca_rate_difference(r, lepton, t, T, etas.at(lepton)) +
                    hadron_murca_rate_difference(r, lepton, t, T, etas.at(lepton)));
        };
        auto diff_due_inf = [&](double r, const auxiliaries::phys::Species &species)
        {
            return exp_phi(r) * (quark_ud_durca_rate_difference(r, t, T, etas.at(species)) +
                                 quark_ud_murca_rate_difference(r, t, T, etas.at(species)));
        };
        auto diff_sue_inf = [&](double r, const auxiliaries::phys::Species &species)
        {
            return exp_phi(r) * (quark_us_durca_rate_difference(r, t, T, etas.at(species)) +
                                 quark_us_murca_rate_difference(r, t, T, etas.at(species)));
        };
        auto npe_reaction_change = auxiliaries::math::integrate_volume<const auxiliaries::phys::Species &>(std::function<double(double, const auxiliaries::phys::Species &)>(diff_npl_inf), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)(electron),
             npm_reaction_change = auxiliaries::math::integrate_volume<const auxiliaries::phys::Species &>(std::function<double(double, const auxiliaries::phys::Species &)>(diff_npl_inf), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)(muon),
             due_reaction_change = auxiliaries::math::integrate_volume<const auxiliaries::phys::Species &>(std::function<double(double, const auxiliaries::phys::Species &)>(diff_due_inf), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)(uquark),
             sue_reaction_change = auxiliaries::math::integrate_volume<const auxiliaries::phys::Species &>(std::function<double(double, const auxiliaries::phys::Species &)>(diff_sue_inf), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)(squark);

        auto e_particle_change_H = npe_reaction_change,
             m_particle_change = npm_reaction_change,
             u_particle_change = due_reaction_change + sue_reaction_change,
             s_particle_change = -sue_reaction_change;

        auto rotochemical_luminosity = e_particle_change_H * etas.at(electron) +
                                       m_particle_change * etas.at(muon) +
                                       u_particle_change * etas.at(uquark) +
                                       s_particle_change * (etas.at(uquark) - etas.at(squark));

        // rhs for the temperature balance equation
        evaluated_rhs.push_back(-(photon_luminosity(t, T) + neutrino_luminosity(t, T, etas) - rotochemical_luminosity) / heat_capacity(t, T));

        // make particle change relevant for convoluting with Z matrix
        std::vector<double> del_particle_changes{e_particle_change_H - (i_omegas.at(electron) - i_omegas.at(uquark)) * omega_sqr_dot(t),
                                                 m_particle_change - i_omegas.at(muon) * omega_sqr_dot(t),
                                                 u_particle_change - i_omegas.at(uquark) * omega_sqr_dot(t),
                                                 s_particle_change - i_omegas.at(squark) * omega_sqr_dot(t)};

        // rhs for the etas
        for (size_t species_count = 0; species_count < rh_particles.size(); ++species_count)
        {
            auto species = rh_particles[species_count];
            double value = 0;

            for (size_t j = 0; j < zij.columns(); ++j)
            {
                value += zij.at(species_count, j) * del_particle_changes[j];
            }
            evaluated_rhs.push_back(value);
        }
        return evaluated_rhs;
    };

    // solve cooling equations
    double exp_phi_at_R = pow(1 - 2 * constants::scientific::G * m_ns / r_ns, 0.5);

    // initial temperature profile
    std::vector<double> radii, profile;
    for (double r = cooling_radius_step / 2.0; r < r_ns; r += cooling_radius_step)
    {
        radii.push_back(r);
        profile.push_back(initial_t_profile_inf(r, r_ns, exp_phi, nbar));
    }
    double t_step = base_t_step,
           t_curr = t_init;
    // instantiate initial values for (T, etas)
    std::vector<double> previous_values(1 + rh_particles.size(), 0);
    previous_values[0] = profile.end()[-2];
    // initial chemical imbalances are zero

    // solution arrays
    std::vector<double> time, surface_temp;
    std::vector<std::vector<double>> others(save_chemical_imbalances ? 2 + rh_particles.size() : 2);
    time.reserve(cooling_n_points_estimate);
    surface_temp.reserve(cooling_n_points_estimate);
    for (size_t i = 0; i < others.size(); ++i)
    {
        others[i].reserve(cooling_n_points_estimate);
    }

    size_t indent = 20;
    std::cout << "M = " << m_ns * constants::conversion::gev_over_msol << " [Ms]\n";
    std::cout << std::left << std::setw(indent) << "t [years] "
              << std::setw(indent) << "Te^inf [K] "
              << std::setw(indent) << "L^inf_ph [erg/s] "
              << std::setw(indent) << "L^inf_nu [erg/s] ";
    if (save_chemical_imbalances)
    {
        for (auto rh_species = rh_particles.begin(); rh_species != rh_particles.end(); ++rh_species)
        {
            std::stringstream ss;
            ss << "eta^inf_" << rh_species->name() << " [K]";
            std::cout << std::left << std::setw(indent) << ss.str();
        }
    }
    std::cout << '\n';

    while (t_curr < t_end)
    {
        auto values = cooling::solver::coupled_cooling(t_curr, t_step, rotochemical_vector_rhs, previous_values, cooling_newton_step_eps, cooling_newton_max_iter);
        double max_diff = 0;
        for (size_t i = 0; i < values.size(); ++i)
        {
            max_diff = (previous_values[i] != 0) ? std::max(max_diff, std::abs((values[i] - previous_values[i]) / std::max(previous_values[0], std::abs(previous_values[i])))) : max_diff;
        }
        if (max_diff > cooling_max_diff_per_t_step)
        {
            t_step /= 2.0;
            continue;
        }
        previous_values = values;
        t_curr += t_step;
        t_step *= exp_rate_estim;

        // print in understandable units
        // std::cout << 1.0E6 * t_curr / (constants::conversion::myr_over_s * constants::conversion::gev_s) << "\t" << auxiliaries::phys::te_tb_relation(y.back(), r_ns, m_ns, crust_eta) * exp_phi_at_R * constants::conversion::gev_over_k << "\t" << photon_luminosity(t_curr, y.back()) * constants::conversion::gev_s / constants::conversion::erg_over_gev << "\t" << (-rotochemical_vector_rhs(t_curr, previous_values)[0] * heat_capacity(t_curr, y.back()) - photon_luminosity(t_curr, y.back())) * constants::conversion::gev_s / constants::conversion::erg_over_gev << "\t" << '\n';
        // save in understandable units
        time.push_back(1.0E6 * t_curr / (constants::conversion::myr_over_s * constants::conversion::gev_s));
        surface_temp.push_back(auxiliaries::phys::te_tb_relation(values[0], r_ns, m_ns, crust_eta) * exp_phi_at_R * constants::conversion::gev_over_k);
        others[0].push_back(photon_luminosity(t_curr, values[0]) * constants::conversion::gev_s / constants::conversion::erg_over_gev);
        others[1].push_back((-rotochemical_vector_rhs(t_curr, previous_values)[0] * heat_capacity(t_curr, values[0]) - photon_luminosity(t_curr, values[0])) * constants::conversion::gev_s / constants::conversion::erg_over_gev);
        if (save_chemical_imbalances)
        {
            for (size_t i = 0; i < rh_particles.size(); ++i)
            {
                others[i + 2].push_back(values[i + 1] * constants::conversion::gev_over_k);
            }
        }
        // print
        std::cout << std::left << std::setw(indent) << time.back() << std::setw(indent) << surface_temp.back() << std::setw(indent) << others[0].back() << std::setw(indent) << others[1].back();
        if (save_chemical_imbalances)
        {
            for (size_t i = 0; i < rh_particles.size(); ++i)
            {
                std::cout << std::left << std::setw(indent) << others[i + 2].back();
            }
        }
        std::cout << '\n';
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
        if (save_chemical_imbalances)
        {
            for (size_t i = 0; i < rh_particles.size(); ++i)
            {
                auto gr_eta = new TGraph(time.size(), time.data(), others[i + 2].data());
                gr_eta->GetXaxis()->SetTitle("t [yr]");
                std::stringstream ss;
                ss << "#eta^{#infty}_{" << rh_particles[i].name() << "} [K]";
                gr_eta->GetYaxis()->SetTitle(ss.str().c_str());
                ss.str(std::string());
                ss << "eta_" << rh_particles[i].name();
                rootfile->WriteObject(gr_eta, ss.str().c_str());
            }
        }
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