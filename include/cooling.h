#ifndef COOLING_H
#define COOLING_H

#include "../include/auxiliaries.h"

#include <vector>
#include <functional>
#include <map>
#include <string>

/// @brief Cooling related namespace
namespace cooling
{
    /// @brief Cooling equation solver
    namespace solver
    {
        /// @brief Solves the cooling equation dT^inf/dt = F(t,T^inf) for a rhs function and a given initial temperature (inf means in distant frame)
        /// @param cache cache support via CachedFunction wrapper; Contains of (t, T) pairs in exponential differencing scheme; If t > last cached time, gets reevaluated!
        /// @param t time at which the temperature is to be calculated [GeV^{-1}]
        /// @param cooling_rhs right hand side term F(t,T^inf) of the cooling equation, i.e. L^inf/Cv^inf [GeV^{2}]
        /// @param initial_temperature Temperature^inf at t=0 [GeV]
        /// @param base_time_step initial time step [GeV^{-1}] 
        /// @param exp_rate a constant multiplied applied to the time step at each iteration to cover multiple timescales
        /// @param interpolator interpolation function for the cache. It is safe to use cached version here. Remember to have sufficient amount of points for the interpolation
        /// @return Temperature^inf at given t [GeV]
        double stationary_cooling_cached(std::vector<std::vector<double>> &cache, double t, const std::function<double(double, double)> &cooling_rhs, double initial_temperature, double base_time_step, double exp_rate,
            const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator = 
            [](const std::vector<double> &x, const std::vector<double> &y, double val) { return auxiliaries::interpolate(x, y, auxiliaries::InterpolationMode::kLinear, val); });
    }
    
    /// @brief Predefined functionality, including cooling luminosities
    namespace predefined
    {
        /// @brief Specific heat calculator
        namespace specific_heat
        {
            /// @brief specific heat calculator based on Fermi gas model cv = sum (m_star * k_fermi)/3 * T
            /// @param cache cache support via CachedFunction wrapper; holds cache[0]=Cv/T^inf
            /// @param m_star_functions m_star functions for each species (GeV); this dictionary defines considered species by its keys
            /// @param k_fermi_functions k_fermi functions for each species (GeV); must have (at least) the same keys as m_star_functions
            /// @param nbar_of_r nbar(r) function in the NS (datafile units)
            /// @param exp_lambda_of_r exp(lambda(r)) function in the NS
            /// @param exp_phi_of_r exp(phi(r)) function in the NS
            /// @param r_ns NS radius [GeV^{-1}]
            /// @param radius_step radius step for the integration [GeV^{-1}]
            /// @return Cv(t, T^inf) [GeV^{0}] (integrated)
            std::function<double(double, double)> fermi_specific_heat_cached(std::vector<double> &cache, const std::map<std::string, std::function<double(double)>> &m_star_functions, const std::map<std::string, std::function<double(double)>> &k_fermi_functions, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_lambda_of_r, const std::function<double(double)> &exp_phi_of_r, double r_ns, double radius_step);
        }

        /// @brief Photon related cooling functionality
        namespace photonic
        {
            /// @brief Returns the cooling luminosity function of (t, T^inf) for a given temperature, radius, mass, light element share and phi function at radius R. The calculation is based on the formula Keisuke thesis.
            /// @param R NS radius [GeV^{-1}]
            /// @param M NS mass [GeV]
            /// @param eta g_14^2 * delta M / M , where g_14 is surface_gravity/(10^14 cm/s^2) and delta M is the light element mass on the surface
            /// @return cooling luminosity function of (t, T^inf) [GeV^{2}]
            std::function<double(double, double)> surface_luminosity(double R, double M, double eta);
        }

        /// @brief Neutrino related cooling functionality
        namespace neutrinic
        {
            /// @brief Returns the cooling luminosity function of (t, T^inf) for a given temperature, (n, p, l) m*'s and k_fermi's. The calculation is based on the formula Keisuke thesis.
            /// @param cache cache support via CachedFunction wrapper; holds cache[0]=L^inf/(T^inf)^6
            /// @param m_star_n neutron m* [GeV]
            /// @param m_star_p proton m* [GeV]
            /// @param m_star_l lepton m* [GeV]
            /// @param k_fermi_n neutron k_fermi [GeV]
            /// @param k_fermi_p proton k_fermi [GeV]
            /// @param k_fermi_l lepton k_fermi [GeV]
            /// @param nbar_of_r nbar(r) function in the NS (datafile units)
            /// @param exp_lambda_of_r exp(lambda(r)) function in the NS
            /// @param exp_phi_of_r exp(phi(r)) function in the NS
            /// @param r_ns NS radius [GeV^{-1}]
            /// @param radius_step radius step for the integration [GeV^{-1}]
            /// @return cooling luminosity function of (t, T^inf) [GeV^{2}]
            /// @note Covers processes: n -> p + l + anti-nu_l, p + l -> n + nu_l
            std::function<double(double, double)> hadron_durca_luminocity_cached(std::vector<double> &cache, const std::function<double(double)> &m_star_n, const std::function<double(double)> &m_star_p, const std::function<double(double)> &m_star_l, const std::function<double(double)> &k_fermi_n, const std::function<double(double)> &k_fermi_p, const std::function<double(double)> &k_fermi_l, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_lambda_of_r, const std::function<double(double)> &exp_phi_of_r, double r_ns, double radius_step);

            /// @brief Returns the cooling luminosity function of (t, T^inf) for a given temperature, (n, p, l) m*'s and k_fermi's. The calculation is based on the formula Keisuke thesis.
            /// @param cache cache support via CachedFunction wrapper; holds cache[0]=L^inf/(T^inf)^8
            /// @param m_star_n neutron m* [GeV]
            /// @param m_star_p proton m* [GeV]
            /// @param m_star_l lepton m* [GeV]
            /// @param k_fermi_n neutron k_fermi [GeV]
            /// @param k_fermi_p proton k_fermi [GeV]
            /// @param k_fermi_l lepton k_fermi [GeV]
            /// @param nbar_of_r nbar(r) function in the NS (datafile units)
            /// @param exp_lambda_of_r exp(lambda(r)) function in the NS
            /// @param exp_phi_of_r exp(phi(r)) function in the NS
            /// @param r_ns NS radius [GeV^{-1}]
            /// @param radius_step radius step for the integration [GeV^{-1}]
            /// @param nbar_conversion conversion factor from datafile nbar units to GeV^3
            /// @return cooling luminosity function of (t, T^inf) [GeV^{2}]
            /// @note Covers processes: n + N1 -> p + N2 + l + anti-nu_l, p + N2 + l -> n + N1 + nu_l
            std::function<double(double, double)> hadron_murca_luminocity_cached(std::vector<double> &cache, const std::function<double(double)> &m_star_n, const std::function<double(double)> &m_star_p, const std::function<double(double)> &m_star_l, const std::function<double(double)> &k_fermi_n, const std::function<double(double)> &k_fermi_p, const std::function<double(double)> &k_fermi_l, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_lambda_of_r, const std::function<double(double)> &exp_phi_of_r, double r_ns, double radius_step, double nbar_conversion);

            /// @brief Returns the cooling luminosity function of (t, T^inf) for a given temperature, (n, p, l) m*'s and k_fermi's. The calculation is based on https://arxiv.org/pdf/astro-ph/0012122.pdf.
            /// @param cache cache support via CachedFunction wrapper; holds cache[0]=L^inf/(T^inf)^8
            /// @param m_star_n neutron m* [GeV]
            /// @param m_star_p proton m* [GeV]
            /// @param k_fermi_n neutron k_fermi [GeV]
            /// @param k_fermi_p proton k_fermi [GeV]
            /// @param nbar_of_r nbar(r) function in the NS (datafile units)
            /// @param exp_lambda_of_r exp(lambda(r)) function in the NS
            /// @param exp_phi_of_r exp(phi(r)) function in the NS
            /// @param r_ns NS radius [GeV^{-1}]
            /// @param radius_step radius step for the integration [GeV^{-1}]
            /// @return cooling luminosity function of (t, T^inf) [GeV^{2}]
            /// @note Covers processes: n + n -> n + n + nu_l + anti-nu_l, p + p -> p + p + nu_l + anti-nu_l, n + p -> n + p + nu_l + anti-nu_l for (2?) neutrino sorts
            std::function<double(double, double)> hadron_bremsstrahlung_luminocity_cached(std::vector<double> &cache, const std::function<double(double)> &m_star_n, const std::function<double(double)> &m_star_p, const std::function<double(double)> &k_fermi_n, const std::function<double(double)> &k_fermi_p, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_lambda_of_r, const std::function<double(double)> &exp_phi_of_r, double r_ns, double radius_step);
        }
    }
}

#endif // COOLING_H