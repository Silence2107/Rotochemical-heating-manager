
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

int main(int argc, char **argv)
{
    std::string program_name = "dominant_process";
    argparse::ArgumentParser parser(program_name, "Determines the dominant cooling processes for every timeslice", "Argparse powered by SiLeader");

    parser.addArgument({"--inputfile"}, "json input file path (required)");
    parser.addArgument({"--fraction"}, "what fraction of total luminosity is assumed significant (double, defaults to 0.1)");
    auto args = parser.parseArgs(argc, argv);

    using namespace instantiator;
    instantiator::instantiate_system(args.get<std::string>("inputfile"), {"TOV", "COOL"});

    auxiliaries::io::Logger logger(program_name);

    double dominant_fraction = args.safeGet<double>("fraction", 0.1);

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
    // to find dominant processes, I will also have to split species-specific processes
    std::function<double(double, double, double)> hadron_murca_emissivity_electron = [&](double r, double t, double T)
    {
        return hadron_murca_emissivity(r, constants::species::electron, t, T);
    };

    std::function<double(double, double, double)> hadron_murca_emissivity_muon = [&](double r, double t, double T)
    {
        return hadron_murca_emissivity(r, constants::species::muon, t, T);
    };

    std::function<double(double, double, double)> hadron_durca_emissivity_electron = [&](double r, double t, double T)
    {
        return hadron_durca_emissivity(r, constants::species::electron, t, T);
    };

    std::function<double(double, double, double)> hadron_durca_emissivity_muon = [&](double r, double t, double T)
    {
        return hadron_durca_emissivity(r, constants::species::muon, t, T);
    };

    std::function<double(double, double, double)> hadron_PBF_emissivity_proton = [&](double r, double t, double T)
    {
        return hadron_PBF_emissivity(r, constants::species::proton, t, T);
    };

    std::function<double(double, double, double)> hadron_PBF_emissivity_neutron = [&](double r, double t, double T)
    {
        return hadron_PBF_emissivity(r, constants::species::proton, t, T);
    };

    std::map<std::string, std::function<double(double, double, double)>> processes =
        {{"hMUnpe", hadron_murca_emissivity_electron},
         {"hMUnpmu", hadron_murca_emissivity_muon},
         {"hDUnpe", hadron_durca_emissivity_electron},
         {"hDUnpmu", hadron_durca_emissivity_muon},
         {"hPBFp", hadron_PBF_emissivity_proton},
         {"hPBFn", hadron_PBF_emissivity_neutron},
         {"hhBREMS", hadron_bremsstrahlung_emissivity},
         {"qDUude", quark_ud_durca_emissivity},
         {"qDUuse", quark_us_durca_emissivity},
         {"qMUude", quark_ud_murca_emissivity},
         {"qMUuse", quark_us_murca_emissivity},
         {"qqBREMS", quark_bremsstrahlung_emissivity},
         {"eeBREMS", electron_bremsstrahlung_emissivity}}; // available processes

    auto Q_nu = [&](double r, double t, double T)
    {
        using namespace constants::scientific;
        double result = 0;
        for (auto it = processes.begin(); it != processes.end(); ++it)
            result += it->second(r, t, T);
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
    std::vector<std::vector<double>> others(3);
    time.reserve(cooling_n_points_estimate);
    surface_temp.reserve(cooling_n_points_estimate);
    others[0].reserve(cooling_n_points_estimate);
    others[1].reserve(cooling_n_points_estimate);
    others[2].reserve(cooling_n_points_estimate);

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
              << std::setw(indent) << "L^inf_nu[erg/s] "
              << std::setw(indent + 5) << "L^inf_nu_domin[erg/s]"
              << std::setw(indent) << "Dominant_processes " << '\n';

    while (t_curr < t_end)
    {
        double next_T;                             // predicted T
        bool reached_adaption_limit = false;       // control for adaptive solver
        bool reached_negative_temperature = false; // exclude NaN generation to "negative" temperature

        std::map<std::string, double> neutrino_luminosities; // placeholder for processes neutrino luminosity estimates

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

            for (auto it = processes.begin(); it != processes.end(); ++it)
            {
                neutrino_luminosities[it->first] = 0.0;
                for (size_t count = 0; count < radii.size() - 1; ++count)
                    neutrino_luminosities[it->first] += 4 * constants::scientific::Pi * radii[count] * radii[count] * exp_lambda(radii[count]) *
                                                        exp_phi(radii[count]) * exp_phi(radii[count]) * (radii[count + 1] - radii[count]) * it->second(radii[count], t_curr + t_step, profile[count]);
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
                               ss << std::scientific << std::setprecision(3) << "At t = " << 1.0E6 * t_curr / (constants::conversion::myr_over_s * constants::conversion::gev_s) << " [yr] target temperature difference exceeded (" << max_diff << " > " << cooling_max_diff_per_t_step << "), adapting (raise CoolingSolver.StepTolerance/NewtonTolerance if this halts progress)";
                               return ss.str();
                           },
                           "eq. cooling");
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

            for (auto it = processes.begin(); it != processes.end(); ++it)
            {
                neutrino_luminosities[it->first] = 0.0;
                for (size_t count = 0; count < radii.size() - 1; ++count)
                    neutrino_luminosities[it->first] += 4 * constants::scientific::Pi * radii[count] * radii[count] * exp_lambda(radii[count]) *
                                                        exp_phi(radii[count]) * exp_phi(radii[count]) * (radii[count + 1] - radii[count]) * it->second(radii[count], t_curr + t_step, next_T);
            }
        }
        temp_curr = next_T;
        t_curr += t_step;
        t_step *= exp_rate_estim;

        double neutrino_lum = 0.0,
               neutrino_lum_dominant = 0.0;
        std::vector<std::string> dominant_processes;
        for (auto it = neutrino_luminosities.begin(); it != neutrino_luminosities.end(); ++it)
        {
            neutrino_lum += it->second;
        }
        for (auto it = neutrino_luminosities.begin(); it != neutrino_luminosities.end(); ++it)
        {
            if (dominant_fraction <= it->second / neutrino_lum)
            {
                neutrino_lum_dominant += it->second;
                dominant_processes.push_back(it->first);
            }
        }

        // save in understandable units
        time.push_back(1.0E6 * t_curr / (constants::conversion::myr_over_s * constants::conversion::gev_s));
        surface_temp.push_back(auxiliaries::phys::te_tb_relation(temp_curr, r_ns, m_ns, crust_eta) * exp_phi_at_R * constants::conversion::gev_over_k);
        others[0].push_back(photon_luminosity(t_curr, temp_curr) * constants::conversion::gev_s / constants::conversion::erg_over_gev);
        others[1].push_back(neutrino_lum * constants::conversion::gev_s / constants::conversion::erg_over_gev);
        others[2].push_back(neutrino_lum_dominant * constants::conversion::gev_s / constants::conversion::erg_over_gev);
        
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
        std::cout << std::left << std::setw(indent) << time.back() << std::setw(indent) << surface_temp.back() << std::setw(indent) << others[0].back() << std::setw(indent) << others[1].back() << std::setw(indent) << others[2].back();
        for (auto &elem : dominant_processes)
        {
            std::cout << elem << ",";
        }
        std::cout << '\n';
    }
}