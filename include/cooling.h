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
        /// @param t time at which the temperature is to be calculated [GeV^{-1}] (initial time is assumed to be 0)
        /// @param cooling_rhs right hand side term F(t,T^inf) of the cooling equation, i.e. L^inf/Cv^inf [GeV^{2}]
        /// @param initial_temperature Temperature^inf at t=0 [GeV]
        /// @param base_time_step initial time step [GeV^{-1}]
        /// @param exp_rate a constant multiplying the time step at each iteration to cover multiple timescales
        /// @param interpolator interpolation function for the cache. It is safe to use cached version here. Remember to have sufficient amount of points for the interpolation
        /// @return Temperature^inf at given t [GeV]
        double stationary_cooling_cached(
            std::vector<std::vector<double>> &cache, double t, const std::function<double(double, double)> &cooling_rhs, double initial_temperature, double base_time_step, double exp_rate,
            const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator =
                [](const std::vector<double> &x, const std::vector<double> &y, double val)
            { return auxiliaries::interpolate(x, y, auxiliaries::InterpolationMode::kLinear, val); });
    }

    /// @brief Predefined functionality, including cooling luminosities
    namespace predefined
    {
        /// @brief Auxiliary physics related to cooling
        namespace auxiliary
        {
            /// @brief Te-Tb relation, based on Keisure thesis
            /// @param Tb internal temperature [GeV], measured by distant observer (inf)
            /// @param R NS radius [GeV^{-1}]
            /// @param M NS mass [GeV]
            /// @param eta NS light element share
            /// @return corresponding surface temperature [GeV]
            double te_tb_relation(double Tb, double R, double M, double eta);

            /// @brief Nucleon critical temperature parametrization from Keisuke thesis, (2.34)
            /// @param k_fermi fermi momentum [GeV], correspoding to the nucleon
            /// @param temp_ampl amplitude parameter [GeV]
            /// @param k_offs mean-like parameter [GeV]
            /// @param k_width variance-like parameter [GeV]
            /// @param quad_skew high order correction to the gaussian [dimensionless]
            /// @return T0 * exp[- ((k_fermi - k_offs) / (k_width))^2 - quad_skew * ((k_fermi - k_offs) / (k_width))^4]
            double critical_temperature_smeared_guassian(double k_fermi, double temp_ampl, double k_offs, double k_width, double quad_skew);

            /// @brief Enumerate choices for the nucleon critical temperature parametrization
            enum class CriticalTemperatureModel
            {
                kAO,
                kCCDK,
                kA,
                kB,
                kC,
                kA2
            };

            /// @brief Critical temperature in given model
            /// @param k_fermi fermi momentum [GeV], correspoding to the nucleon
            /// @param model choice of the model
            double critical_temperature(double k_fermi, CriticalTemperatureModel model);

            /// @brief Superfluid gap in 1S0 state (A)
            /// @param tau T/Tc
            double superfluid_gap_1s0(double tau);

            /// @brief Superfluid gap in 3P2 state (B)
            /// @param tau T/Tc
            double superfluid_gap_3p2(double tau);
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
            /// @brief Emissivity of neutrinos from hadronic DUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [data units]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [data units]
            /// @param nbar_of_r baryon density [data units] as a function of radius [GeV^{-1}]
            /// @param nbar_core_limit baryon density [data units] at the core
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_n_1s0 allow/forbid superfluidity in 1S0 state for neutrons
            /// @param superfluid_p_1s0 allow/forbid superfluidity in 1S0 state for protons
            /// @param superfluid_n_3p2 allow/forbid superfluidity in 3P2 state for neutrons
            /// @param superfluid_p_temp temperature of superfluid protons [GeV]
            /// @param superfluid_n_temp temperature of superfluid neutrons [GeV]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV]
            std::function<double(double, const std::string &, double, double)> hadron_durca_emissivity(
                const std::map<std::string, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<std::string, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from hadronic MUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [data units]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [data units]
            /// @param nbar_of_r baryon density [data units] as a function of radius [GeV^{-1}]
            /// @param nbar_core_limit baryon density [data units] at the core
            /// @param nbar_conversion conversion factor from baryon density [data units] to baryon density [GeV^{-3}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_n_1s0 allow/forbid superfluidity in 1S0 state for neutrons
            /// @param superfluid_p_1s0 allow/forbid superfluidity in 1S0 state for protons
            /// @param superfluid_n_3p2 allow/forbid superfluidity in 3P2 state for neutrons
            /// @param superfluid_p_temp temperature of superfluid protons [GeV]
            /// @param superfluid_n_temp temperature of superfluid neutrons [GeV]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV]
            std::function<double(double, const std::string &, double, double)> hadron_murca_emissivity(
                const std::map<std::string, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<std::string, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_core_limit, double nbar_conversion, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from hadronic bremsstrahlung reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [data units]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [data units]
            /// @param nbar_of_r baryon density [data units] as a function of radius [GeV^{-1}]
            /// @param ion_volume_frac volume fraction of ions as a function of baryon density [data units]
            /// @param nbar_core_limit baryon density [data units] at the core
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_n_1s0 allow/forbid superfluidity in 1S0 state for neutrons
            /// @param superfluid_p_1s0 allow/forbid superfluidity in 1S0 state for protons
            /// @param superfluid_n_3p2 allow/forbid superfluidity in 3P2 state for neutrons
            /// @param superfluid_p_temp temperature of superfluid protons [GeV]
            /// @param superfluid_n_temp temperature of superfluid neutrons [GeV]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV]
            std::function<double(double, double, double)> hadron_bremsstrahlung_emissivity(
                const std::map<std::string, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<std::string, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &ion_volume_frac, double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from hadronic pair breaking & formations
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [data units]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [data units]
            /// @param nbar_of_r baryon density [data units] as a function of radius [GeV^{-1}]
            /// @param nbar_core_limit baryon density [data units] at the core
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_n_1s0 allow/forbid superfluidity in 1S0 state for neutrons
            /// @param superfluid_p_1s0 allow/forbid superfluidity in 1S0 state for protons
            /// @param superfluid_n_3p2 allow/forbid superfluidity in 3P2 state for neutrons
            /// @param superfluid_p_temp temperature of superfluid protons [GeV]
            /// @param superfluid_n_temp temperature of superfluid neutrons [GeV]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV]
            std::function<double(double, const std::string &, double, double)> hadron_pbf_emissivity(
                const std::map<std::string, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<std::string, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);
        }
    }
}

#endif // COOLING_H