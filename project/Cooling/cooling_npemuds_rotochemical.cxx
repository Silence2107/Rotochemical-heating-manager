
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
    std::string program_name = "cooling_npemuds_rotochemical";
    argparse::ArgumentParser parser(program_name, "Solves the cooling equation coupled with chemical imbalances based on EoS", "Argparse powered by SiLeader");

    parser.addArgument({"--inputfile"}, "json input file path (required)");
#if RHM_HAS_ROOT
    parser.addArgument({"--pdf_path"}, "pdf output file path (optional, default: Cooling-with-RH.pdf)");
    parser.addArgument({"--rootfile_path"}, "root output file path (optional, default: None)");
#endif
    parser.addArgument({"--save_chemical_imbalances"}, "print & write chemical imbalances (optional, value-free, default: false)", argparse::ArgumentType::StoreTrue);

    auto args = parser.parseArgs(argc, argv);

    using namespace instantiator;
    instantiator::instantiate_system(args.get<std::string>("inputfile"), {"TOV", "COOL", "RH"});

    auxiliaries::io::Logger logger(program_name);

#if RHM_HAS_ROOT
    std::string pdf_path = args.safeGet<std::string>("pdf_path", "Cooling-with-RH.pdf");
    TFile *rootfile = nullptr;
    if (args.has("rootfile_path"))
        rootfile = new TFile(args.get<std::string>("rootfile_path").c_str(), "RECREATE");
#endif

    bool save_chemical_imbalances = args.has("save_chemical_imbalances");

    logger.log([]()
               { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
               [&]()
               { return "Instantiation complete"; });
    // RUN --------------------------------------------------------------------------

    // EoS definition

    auto eos_inv_cached = edensity_of_pressure;

    // TOV solver

    auto tov_df = tov_solver::tov_solution(eos_inv_cached, center_pressure, radius_step, surface_pressure, pressure_low, tov_adapt_limit);

    std::vector<double> df_nbar(tov_df[0].size()), df_exp_phi(tov_df[0].size()), df_exp_lambda(tov_df[0].size()), df_pressure(tov_df[0].size()), df_edensity(tov_df[0].size());

    double r_ns = tov_df[0].back();
    double m_ns = tov_df[1].back();
    for (size_t i = 0; i < tov_df[0].size(); ++i)
    {
        df_pressure[i] = tov_df[2][i];
        df_nbar[i] = nbar_of_pressure(tov_df[2][i]);
        df_exp_phi[i] = std::exp(tov_df[3][i]);
        if (i == 0)
            df_exp_lambda[i] = 1.0;
        else
            df_exp_lambda[i] = std::pow(1 - 2 * constants::scientific::G * tov_df[1][i] / tov_df[0][i], -0.5);
        df_edensity[i] = edensity_of_pressure(tov_df[2][i]);
    }
    auto nbar = auxiliaries::math::Interpolator(tov_df[0], df_nbar, radial_interp_mode);
    auto exp_phi = auxiliaries::math::Interpolator(tov_df[0], df_exp_phi, radial_interp_mode);
    auto exp_lambda = auxiliaries::math::Interpolator(tov_df[0], df_exp_lambda, radial_interp_mode);
    auto pressure_of_r = auxiliaries::math::Interpolator(tov_df[0], df_pressure, radial_interp_mode);
    auto edensity_of_r = auxiliaries::math::Interpolator(tov_df[0], df_edensity, radial_interp_mode);

    // Identify particle list, supplied by the user
    std::vector<auxiliaries::phys::Species> rh_particles;

    for (auto it = dni_to_dmuj.begin(); it != dni_to_dmuj.end(); ++it)
    {
        auto key = it->first.first;
        if (std::find(rh_particles.begin(), rh_particles.end(), key) == rh_particles.end())
        {
            rh_particles.push_back(key);
        }
    }

    // integrate the bij integrand over volume
    auxiliaries::math::MatrixD bij(rh_particles.size(), rh_particles.size(), 0);

    for (size_t i = 0; i < rh_particles.size(); ++i)
    {
        for (size_t j = 0; j < rh_particles.size(); ++j)
        {
            auto integrand = [&](double r)
            {
                return dni_to_dmuj.at({rh_particles[i], rh_particles[j]})(nbar(r)) / exp_phi(r);
            };
            bij.at(i, j) = auxiliaries::math::integrate_volume<>(std::function<double(double)>(integrand), 0.0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kRectangular, radius_step / 10)();
        }
    }

    // Clear empty rows, if any
    for (size_t row = 0; row < rh_particles.size(); ++row)
    {
        bool empty_row = true;
        for (size_t col = 0; col < rh_particles.size(); ++col)
        {
            if (bij.at(row, col) != 0)
            {
                empty_row = false;
                break;
            }
        }
        if (empty_row)
        {
            auto bij_updated = auxiliaries::math::MatrixD(rh_particles.size() - 1, rh_particles.size() - 1, 0);
            for (size_t i = 0; i < bij_updated.rows(); ++i)
            {
                for (size_t j = 0; j < bij_updated.columns(); ++j)
                {
                    bij_updated.at(i, j) = bij.at(i + (i >= row), j + (j >= row));
                }
            }
            bij = bij_updated;
            rh_particles.erase(rh_particles.begin() + row);
            --row;
        }
    }

    // species to indices mapping
    std::map<auxiliaries::phys::Species, size_t> species_to_index;
    for (size_t i = 0; i < rh_particles.size(); ++i)
    {
        species_to_index[rh_particles[i]] = i;
    }

    // URCA reaction array
    std::vector<std::string> available_imbalances;
    // npe availability
    if (std::find(rh_particles.begin(), rh_particles.end(), constants::species::neutron) != rh_particles.end() &&
        std::find(rh_particles.begin(), rh_particles.end(), constants::species::proton) != rh_particles.end() &&
        std::find(rh_particles.begin(), rh_particles.end(), constants::species::electron) != rh_particles.end())
    {
        available_imbalances.push_back("npe");
    }
    // npmu availability
    if (std::find(rh_particles.begin(), rh_particles.end(), constants::species::neutron) != rh_particles.end() &&
        std::find(rh_particles.begin(), rh_particles.end(), constants::species::proton) != rh_particles.end() &&
        std::find(rh_particles.begin(), rh_particles.end(), constants::species::muon) != rh_particles.end())
    {
        available_imbalances.push_back("npm");
    }
    // due availability
    if (std::find(rh_particles.begin(), rh_particles.end(), constants::species::dquark) != rh_particles.end() &&
        std::find(rh_particles.begin(), rh_particles.end(), constants::species::uquark) != rh_particles.end() &&
        std::find(rh_particles.begin(), rh_particles.end(), constants::species::electron) != rh_particles.end())
    {
        available_imbalances.push_back("due");
    }
    // sue availability
    if (std::find(rh_particles.begin(), rh_particles.end(), constants::species::squark) != rh_particles.end() &&
        std::find(rh_particles.begin(), rh_particles.end(), constants::species::uquark) != rh_particles.end() &&
        std::find(rh_particles.begin(), rh_particles.end(), constants::species::electron) != rh_particles.end())
    {
        available_imbalances.push_back("sue");
    }

    // attempt to invert dni_to_dmuj
    if (bij.det() == 0.0)
        RHM_ERROR("Rotochemical B-matrix is singular, cannot invert.");
    auxiliaries::math::MatrixD zij = bij.inverse();

    auto zij_named = [&](const auxiliaries::phys::Species &i, const auxiliaries::phys::Species &j)
    {
        return zij.at(species_to_index.at(i), species_to_index.at(j));
    };

    // log zij
    logger.log([]()
               { return true; }, auxiliaries::io::Logger::LogLevel::kInfo,
                [&]()
                {
                    std::ostringstream oss;
                    oss << "\nRotochemical B-matrix:\n";
                    for (size_t i = 0; i < rh_particles.size(); ++i)
                    {
                        for (size_t j = 0; j < rh_particles.size(); ++j)
                            oss << std::setw(10) << std::setprecision(4) << bij.at(i, j) << " ";
                        oss << "\n";
                    }
                    oss << "Rotochemical Z-matrix:\n";
                    for (size_t i = 0; i < rh_particles.size(); ++i)
                    {
                        for (size_t j = 0; j < rh_particles.size(); ++j)
                            oss << std::setw(10) << std::setprecision(4) << zij.at(i, j) << " ";
                        oss << "\n";
                    }
                    return oss.str();
                });

    if (std::find(available_imbalances.begin(), available_imbalances.end(), "npe") != available_imbalances.end())
    {    
        auto z_npe = zij_named(constants::species::neutron, constants::species::neutron) + zij_named(constants::species::proton, constants::species::proton) + zij_named(constants::species::electron, constants::species::electron) - 2 * zij_named(constants::species::neutron, constants::species::proton);
        logger.log([=]()
                { return z_npe < 0.0; }, auxiliaries::io::Logger::LogLevel::kDebug,
                [&]()
                {
                    return "Rotochemical stability compromised as Z(npe) = Z_nn + Z_pp + Z_ee - 2Z_np < 0.0. Expect unexpected behavior in cooling curve.";
                });
    }

    // return 0;

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
        
        // integration of -P / Omega_k^2 n_B 4pir^2 exp(lambda(r)) / P'_r dY_i
        auto integrand_wrt_Y = [&](double r)
        {
            using namespace constants::scientific;
            // P'_r comes from TOV
            double p = pressure_of_r(r);
            double rho = edensity_of_r(r);
            double exp_lam = exp_lambda(r);
            double m_encl = r / (2 * G) * (1 - 1 / (exp_lam * exp_lam));
            double P_r_prime = -G * (rho + p) * (m_encl + 4 * Pi * r * r * r * p) * exp_lam * exp_lam / (r * r);
            return -p / omega_k_sqr * nbar(r) * 4 * Pi * r * r * exp_lam / P_r_prime;
        };
        // manual integration wrt Y_i
        for (double r = radius_step / 10; r < r_ns - radius_step; r += radius_step / 10)
        {
            double nbar_r = nbar(r),
                   nbar_r_plus = nbar(r + radius_step / 10);
            double Y_i_r = n_i(nbar_r) / nbar_r;
            double Y_i_r_plus = n_i(nbar_r_plus) / nbar_r_plus;
            double dY_i = Y_i_r_plus - Y_i_r;
            i_omegas[species->first] += integrand_wrt_Y(r) * dY_i;
        }
    }
    for (auto rh_species = constants::species::known_particles.begin(); rh_species != constants::species::known_particles.end(); ++rh_species)
    {
        if (i_omegas.find(*rh_species) == i_omegas.end())
        {
            i_omegas[*rh_species] = 0;
        }
    }
    
    auto Q_nu = [&](double r, double t, double T, const std::map<std::string, double> &etas)
    {
        using namespace constants::scientific;
        using namespace constants::species;
        double result = 0;

        // hadron processes
        for (auto it = m_stars_of_nbar.begin(); it != m_stars_of_nbar.end(); ++it)
        {
            auto key = it->first;
            if (key == electron)
            {
                result += hadron_murca_enhanced_emissivity(r, key, t, T, etas.at("npe"));
                result += hadron_durca_enhanced_emissivity(r, key, t, T, etas.at("npe"));
            }
            if (key == muon)
            {
                result += hadron_murca_enhanced_emissivity(r, key, t, T, etas.at("npm"));
                result += hadron_durca_enhanced_emissivity(r, key, t, T, etas.at("npm"));
            }
            if (key.classify() == auxiliaries::phys::Species::ParticleClassification::kBaryon)
            {
                result += hadron_PBF_emissivity(r, key, t, T);
            }
        }
        result += hadron_bremsstrahlung_emissivity(r, t, T);
        // quark processes
        result += quark_ud_durca_enhanced_emissivity(r, t, T, etas.at("due")) +
                  quark_us_durca_enhanced_emissivity(r, t, T, etas.at("sue")) +
                  quark_ud_murca_enhanced_emissivity(r, t, T, etas.at("due")) +
                  quark_us_murca_enhanced_emissivity(r, t, T, etas.at("sue")) +
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
    auto neutrino_luminosity = auxiliaries::math::integrate_volume<double, double, const std::map<std::string, double> &>(
        std::function<double(double, double, double, const std::map<std::string, double> &)>(Q_nu), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p, radius_step);

    auto heat_capacity = auxiliaries::math::integrate_volume<double, double>(
        std::function<double(double, double, double)>(fermi_specific_heat_dens), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p, radius_step);

    auto rotochemical_vector_rhs = [&](double t, const std::vector<double> &funcs)
    {
        using namespace constants::species;
        std::vector<double> evaluated_rhs;
        // funcs[1+] are nontrivial chem. imbalances, ordered by available_imbalances
        std::map<std::string, double> etas{
            {"npe", 0},
            {"npm", 0},
            {"due", 0},
            {"sue", 0}};
        for (size_t i = 0; i < available_imbalances.size(); ++i)
        {
            etas[available_imbalances[i]] = funcs[i + 1];
        }

        // funcs[0] is temperature
        double T = funcs[0];

        using namespace constants::species;

        // rotochemical heating in npemuds
        auto diff_npe_inf = [&](double r)
        {
            return exp_phi(r) * (hadron_durca_rate_difference(r, constants::species::electron, t, T, etas.at("npe")) +
                                 hadron_murca_rate_difference(r, constants::species::electron, t, T, etas.at("npe")));
        };
        auto diff_npm_inf = [&](double r)
        {
            return exp_phi(r) * (hadron_durca_rate_difference(r, constants::species::muon, t, T, etas.at("npm")) +
                                 hadron_murca_rate_difference(r, constants::species::muon, t, T, etas.at("npm")));
        };
        auto diff_due_inf = [&](double r)
        {
            return exp_phi(r) * (quark_ud_durca_rate_difference(r, t, T, etas.at("due")) +
                                 quark_ud_murca_rate_difference(r, t, T, etas.at("due")));
        };
        auto diff_sue_inf = [&](double r)
        {
            return exp_phi(r) * (quark_us_durca_rate_difference(r, t, T, etas.at("sue")) +
                                 quark_us_murca_rate_difference(r, t, T, etas.at("sue")));
        };
        auto npe_reaction_change = auxiliaries::math::integrate_volume<>(std::function<double(double)>(diff_npe_inf), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)(),
             npm_reaction_change = auxiliaries::math::integrate_volume<>(std::function<double(double)>(diff_npm_inf), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)(),
             due_reaction_change = auxiliaries::math::integrate_volume<>(std::function<double(double)>(diff_due_inf), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)(),
             sue_reaction_change = auxiliaries::math::integrate_volume<>(std::function<double(double)>(diff_sue_inf), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p)();

        // net particle change
        std::map<auxiliaries::phys::Species, double> particle_change;

        particle_change[electron] = npe_reaction_change + due_reaction_change + sue_reaction_change;
        particle_change[muon] = npm_reaction_change;
        particle_change[uquark] = due_reaction_change + sue_reaction_change;
        particle_change[squark] = -sue_reaction_change;
        particle_change[proton] = npe_reaction_change + npm_reaction_change;
        particle_change[neutron] = -npe_reaction_change - npm_reaction_change;
        particle_change[dquark] = -due_reaction_change;
        for (auto species = rh_particles.begin(); species != rh_particles.end(); ++species)
        {
            if (particle_change.find(*species) == particle_change.end())
            {
                particle_change[*species] = 0;
            }
        }

        double rotochemical_luminosity = particle_change[electron] * etas.at("npe") + particle_change[muon] * etas.at("npm") + particle_change[uquark] * (etas.at("due") - etas.at("npe")) + particle_change[squark] * (etas.at("due") - etas.at("sue"));

        // rhs for the temperature balance equation
        evaluated_rhs.push_back(-(photon_luminosity(t, T) + neutrino_luminosity(t, T, etas) - rotochemical_luminosity) / heat_capacity(t, T));

        // make particle change relevant for convoluting with Z matrix
        std::vector<double> del_particle_changes;
        for (const auto &species : rh_particles)
        {
            del_particle_changes.push_back(particle_change.at(species) - i_omegas.at(species) * omega_sqr_dot(t));
        }

        // rhs for the etas
        double npe_imb_eq = 0,
                npm_imb_eq = 0,
                due_imb_eq = 0,
                sue_imb_eq = 0;

        for (size_t j = 0; j < zij.columns(); ++j)
        {
            // if npe is available
            if (std::find(available_imbalances.begin(), available_imbalances.end(), "npe") != available_imbalances.end())
            {
                npe_imb_eq += (zij.at(species_to_index.at(neutron), j) - zij.at(species_to_index.at(proton), j) - zij.at(species_to_index.at(electron), j)) * del_particle_changes[j];
                if (j == zij.columns() - 1)
                    evaluated_rhs.push_back(npe_imb_eq);
            }
            // if npm is available
            if (std::find(available_imbalances.begin(), available_imbalances.end(), "npm") != available_imbalances.end())
            {
                npm_imb_eq += (zij.at(species_to_index.at(neutron), j) - zij.at(species_to_index.at(proton), j) - zij.at(species_to_index.at(muon), j)) * del_particle_changes[j];
                if (j == zij.columns() - 1)
                    evaluated_rhs.push_back(npm_imb_eq);
            }
            // if due is available
            if (std::find(available_imbalances.begin(), available_imbalances.end(), "due") != available_imbalances.end())
            {
                due_imb_eq += (zij.at(species_to_index.at(dquark), j) - zij.at(species_to_index.at(uquark), j) - zij.at(species_to_index.at(electron), j)) * del_particle_changes[j];
                if (j == zij.columns() - 1)
                    evaluated_rhs.push_back(due_imb_eq);
            }
            // if sue is available
            if (std::find(available_imbalances.begin(), available_imbalances.end(), "sue") != available_imbalances.end())
            {
                sue_imb_eq += (zij.at(species_to_index.at(squark), j) - zij.at(species_to_index.at(uquark), j) - zij.at(species_to_index.at(electron), j)) * del_particle_changes[j];
                if (j == zij.columns() - 1)
                    evaluated_rhs.push_back(sue_imb_eq);
            }
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
    std::vector<double> previous_values(1 + available_imbalances.size(), 0);
    previous_values[0] = profile.end()[-2];
    // initial chemical imbalances are zero

    // solution arrays
    std::vector<double> time, surface_temp;
    std::vector<std::vector<double>> others(save_chemical_imbalances ? 2 + available_imbalances.size() : 2);
    time.reserve(cooling_n_points_estimate);
    surface_temp.reserve(cooling_n_points_estimate);
    for (size_t i = 0; i < others.size(); ++i)
    {
        others[i].reserve(cooling_n_points_estimate);
    }

    size_t indent = 20;
    std::cout << "M = " << m_ns * constants::conversion::gev_over_msol << " [Ms]\n";
    std::cout << std::left << std::setw(indent) << "t[years] "
              << std::setw(indent) << "Te^inf[K] "
              << std::setw(indent) << "L^inf_ph[erg/s] "
              << std::setw(indent) << "L^inf_nu[erg/s] ";
    if (save_chemical_imbalances)
    {
        for (auto reac = available_imbalances.begin(); reac != available_imbalances.end(); ++reac)
        {
            std::stringstream ss;
            ss << "eta^inf_" << *reac << "[K]";
            std::cout << std::left << std::setw(indent) << ss.str();
        }
    }
    std::cout << '\n';

    while (t_curr < t_end)
    {
        auto coupled_data = cooling::solver::coupled_cooling(t_curr, t_step, rotochemical_vector_rhs, previous_values, cooling_newton_step_eps, cooling_newton_max_iter);
        auto values = coupled_data[0];
        bool reached_adaption_limit = coupled_data[1][0];       // control for adaptive solver
        bool reached_negative_temperature = coupled_data[1][1]; // exclude NaN generation to "negative" temperature

        double max_diff = 0;
        for (size_t i = 0; i < values.size(); ++i)
        {
            max_diff = (previous_values[i] != 0) ? std::max(max_diff, std::abs((values[i] - previous_values[i]) / std::max(previous_values[0], std::abs(previous_values[i])))) : max_diff;
        }
        if (max_diff > cooling_max_diff_per_t_step || reached_adaption_limit || reached_negative_temperature)
        {
            t_step /= 2.0;
            logger.log([&]()
                       { return max_diff > cooling_max_diff_per_t_step; }, auxiliaries::io::Logger::LogLevel::kDebug,
                       [&]()
                       {
                           std::stringstream ss;
                           ss << std::scientific << std::setprecision(3) << "At t = " << 1.0E6 * t_curr / (constants::conversion::myr_over_s * constants::conversion::gev_s) << " [yr] target temperature difference exceeded (" << max_diff << " > " << cooling_max_diff_per_t_step << "), adapting (raise CoolingSolver.StepTolerance/NewtonTolerance if this halts progress)";
                           return ss.str();
                       },
                       "RH cooling");
            logger.log([&]()
                       { return reached_adaption_limit; }, auxiliaries::io::Logger::LogLevel::kDebug,
                       [&]()
                       { return "Adaption limit reached, adapting (raise CoolingSolver.NewtonMaxIter/NewtonTolerance if this reoccurs)"; }, "RH cooling");
            logger.log([&]()
                       { return reached_negative_temperature; }, auxiliaries::io::Logger::LogLevel::kDebug,
                       [&]()
                       { return "Negative temperature reached, adapting (lower CoolingSolver.NewtonTolerance if this reoccurs)"; }, "RH cooling");
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
        others[1].push_back((-rotochemical_vector_rhs(t_curr, values)[0] * heat_capacity(t_curr, values[0]) - photon_luminosity(t_curr, values[0])) * constants::conversion::gev_s / constants::conversion::erg_over_gev);
        if (save_chemical_imbalances)
        {
            for (size_t i = 0; i < available_imbalances.size(); ++i)
            {
                others[i + 2].push_back(values[i + 1] * constants::conversion::gev_over_k);
            }
        }
        logger.log([&]()
                   { return time.size() % 100 == 0; }, auxiliaries::io::Logger::LogLevel::kInfo,
                   [&]()
                   {
                       std::stringstream ss;
                       ss << std::scientific << std::setprecision(3) << std::to_string(time.size()) + " counts past" << " with t = " << time.back() << " [yr]";
                       return ss.str();
                   },
                   "T(t) loop");
        logger.log([&]()
                   { return true; }, auxiliaries::io::Logger::LogLevel::kTrace,
                   [&]()
                   { 
                        std::stringstream ss;
                        ss << std::scientific << std::setprecision(3) << std::to_string(time.size()) + " counts past" << " with t = " << time.back() << " [yr]";
                        return ss.str(); },
                   "T(t) loop");
        // print
        std::cout << std::left << std::setw(indent) << time.back() << std::setw(indent) << surface_temp.back() << std::setw(indent) << others[0].back() << std::setw(indent) << others[1].back();
        if (save_chemical_imbalances)
        {
            for (size_t i = 0; i < available_imbalances.size(); ++i)
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
            for (size_t i = 0; i < available_imbalances.size(); ++i)
            {
                auto gr_eta = new TGraph(time.size(), time.data(), others[i + 2].data());
                gr_eta->GetXaxis()->SetTitle("t [yr]");
                std::stringstream ss;
                ss << "#eta^{#infty}_{" << available_imbalances[i] << "} [K]";
                gr_eta->GetYaxis()->SetTitle(ss.str().c_str());
                ss.str(std::string());
                ss << "eta_" << available_imbalances[i];
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