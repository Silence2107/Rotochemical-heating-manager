
#include "../../include/auxiliaries.h"
#include "../../include/cooling.h"
#include "../../include/constants.h"
#include "../../include/tov_solver.h"
#include "../../include/instantiator.hpp"

#include <vector>
#include <functional>
#include <cmath>
#include <fstream>
#include <sstream>
#include <fstream>

int main()
{
    using namespace instantiator;

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

    // photon luminosity
    auto photon_luminosity = cooling::predefined::photonic::surface_luminosity(r_ns, m_ns, crust_eta);

    // neutrino luminosity
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

    auto quark_murca_emissivity = cooling::predefined::neutrinic::quark_murca_emissivity(
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
        result += quark_ud_durca_emissivity(r, t, T) + quark_us_durca_emissivity(r, t, T) + quark_murca_emissivity(r, t, T) + quark_bremsstrahlung_emissivity(r, t, T);
        result += electron_bremsstrahlung_emissivity(r, t, T);
        return result * exp_phi(r) * exp_phi(r);
    };

    auto neutrino_luminosity = auxiliaries::math::integrate_volume<double, double>(
        std::function<double(double, double, double)>(Q_nu), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p, radius_step);

    // specific heat
    auto fermi_specific_heat_dens = auxiliaries::phys::fermi_specific_heat_density(
        k_fermi_of_nbar, m_stars_of_nbar, nbar, nbar_sf_shift, exp_phi, superfluid_p_temp, superfluid_n_temp, superconduct_q_gap);

    auto heat_capacity = auxiliaries::math::integrate_volume<double, double>(
        std::function<double(double, double, double)>(fermi_specific_heat_dens), 0, r_ns, exp_lambda, auxiliaries::math::IntegrationMode::kGaussLegendre_12p, radius_step);

    std::ofstream out("NormalizedEmissivities.txt");
    out << "r,km\t Q_{eDU}/T^6,erg/cm3s1K6\t Q_{muDU}/T^6,erg/cm3s1K6\t Q_{eMU}/T^8,erg/cm3s1K8\t Q_{muMU}/T^8,erg/cm3s1K8\t Q_{hnubr}/T^8,erg/cm3s1K8\t Q_{pPBF}/T^7,erg/cm3s1K7\t Q_{nPBF}/T^7,erg/cm3s1K7\n";
    for(double r = radius_step/2; r < r_ns; r += radius_step)
    {
        using namespace constants::conversion;
        using namespace constants::scientific;
        using namespace constants::species;
        double T_inf = 2E8 * exp_phi(r) / gev_over_k;
        out << r / km_gev << "\t" << hadron_durca_emissivity(r, electron, 0, T_inf) / erg_over_cm3_s_gev5 / pow(T_inf * gev_over_k / exp_phi(r), 6) << "\t"
            << hadron_durca_emissivity(r, muon, 0, T_inf) / erg_over_cm3_s_gev5 / pow(T_inf * gev_over_k / exp_phi(r), 6) << "\t"
            << hadron_murca_emissivity(r, electron, 0, T_inf) / erg_over_cm3_s_gev5 / pow(T_inf * gev_over_k / exp_phi(r), 8)
            << "\t" << hadron_murca_emissivity(r, muon, 0, T_inf) / erg_over_cm3_s_gev5 / pow(T_inf * gev_over_k / exp_phi(r), 8)
            << "\t" << hadron_bremsstrahlung_emissivity(r, 0, T_inf) / erg_over_cm3_s_gev5 / pow(T_inf * gev_over_k / exp_phi(r), 8)
            << "\t" << hadron_PBF_emissivity(r, proton, 0, T_inf) / erg_over_cm3_s_gev5 / pow(T_inf * gev_over_k / exp_phi(r), 7)
            << "\t" << hadron_PBF_emissivity(r, neutron, 0, T_inf) / erg_over_cm3_s_gev5 / pow(T_inf * gev_over_k / exp_phi(r), 7) << "\n";
    }
}