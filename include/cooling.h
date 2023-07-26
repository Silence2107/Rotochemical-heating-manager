#ifndef COOLING_H
#define COOLING_H

#include "../include/auxiliaries.h"

#include <vector>
#include <functional>
#include <map>

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
            { return auxiliaries::math::interpolate(x, y, auxiliaries::math::InterpolationMode::kLinear, val); });

        /// @brief Evolves the nonequilibrium cooling equation dT^inf/dt = -Qv^inf/cv + e^{-Lambda}/cv * 1/r^2 d/dr[r^2 e^{Phi - Lambda} lambda dT^inf/dr]
        /// @brief via 1 + 1D differencing scheme (backward exponential in time, forward linear in space). BCs are applied automatically, however, the initial profile must be set.
        /// @param t_curr time at which the initial temperature profile was calculated [GeV^{-1}]
        /// @param t_step time step at which the next profiles are to be calculated [GeV^{-1}]
        /// @param neutrino_rate density of neutrino energy loss Qv^inf [GeV^{2}] as a function of r, t, T^inf
        /// @param cv specific heat capacity cv [GeV^{-1}] as a function of r, t, T^inf
        /// @param lambda thermal conductivity [GeV^2] as a function of r, t, T^inf
        /// @param exp_lambda e^Lambda metric function of radius [GeV^{-1}]
        /// @param exp_phi e^Phi metric function of radius [GeV^{-1}]
        /// @param radii reper radii [GeV^{-1}] at which the profiles are calculated
        /// @param initial_profile initial temperature profile T^inf(r, t=t_curr) [GeV] as an array corresponding to radii
        /// @param te_tb Local surface temperature [GeV] as a function of undercrustal T_b^inf [GeV]
        /// @return Temperature profile T^inf(r, t) as an array corresponding to radii
        std::vector<double> nonequilibrium_cooling(
            double t_curr, double t_step, const std::function<double(double, double, double)> &neutrino_rate, const std::function<double(double, double, double)> &cv, const std::function<double(double, double, double)> &lambda,
            const std::function<double(double)> &exp_lambda, const std::function<double(double)> &exp_phi, const std::vector<double> &radii, const std::vector<double> &initial_profile,
            const std::function<double(double)> &te_tb);

        /// @brief Evolves the nonequilibrium cooling equation cv * dT^inf/dt = -Qv^inf + e^{-Lambda} / (4pir^2) d/dr Ld^inf; -lambda dT^inf/dr = Ld^inf/(4pir^2) e^{Lambda-Phi}
        /// @brief via 1 + 1D differencing scheme (backward exponential in time, forward linear in space). BCs are applied automatically, however, the initial profile must be set.
        /// @param t_curr time at which the initial temperature profile was calculated [GeV^{-1}]
        /// @param t_step time step at which the next profiles are to be calculated [GeV^{-1}]
        /// @param neutrino_rate density of neutrino energy loss Qv^inf [GeV^{2}] as a function of r, t, T^inf
        /// @param cv specific heat capacity cv [GeV^{-1}] as a function of r, t, T^inf
        /// @param lambda thermal conductivity [GeV^2] as a function of r, t, T^inf
        /// @param exp_lambda e^Lambda metric function of radius [GeV^{-1}]
        /// @param exp_phi e^Phi metric function of radius [GeV^{-1}]
        /// @param radii reper radii [GeV^{-1}] at which the profiles are calculated
        /// @param initial_profile initial temperature profile T^inf(r, t=t_curr) [GeV] as an array corresponding to radii
        /// @param te_tb Local surface temperature [GeV] as a function of undercrustal T_b^inf [GeV]
        /// @return Temperature, luminosity profiles [T^inf(r, t), Ld^inf(r, t)] as an [2][i_m] array corresponding to radii
        std::vector<double> nonequilibrium_cooling_2(
            double t_curr, double t_step, const std::function<double(double, double, double)> &neutrino_rate, const std::function<double(double, double, double)> &cv, const std::function<double(double, double, double)> &lambda,
            const std::function<double(double)> &exp_lambda, const std::function<double(double)> &exp_phi, const std::vector<double> &radii, const std::vector<double> &initial_profile,
            const std::function<double(double)> &te_tb);
    }
    /// @brief Predefined functionality, including cooling luminosities
    namespace predefined
    {
        /// @brief Photon related cooling functionality
        namespace photonic
        {
            /// @brief Returns the cooling luminosity function of (t, T^inf) for a given temperature, radius, mass, light element share and phi function at radius R.
            /// @param R NS radius [GeV^{-1}]
            /// @param M NS mass [GeV]
            /// @param eta g_14^2 * delta M / M , where g_14 is surface_gravity/(10^14 cm/s^2) and delta M is the light element mass on the surface
            /// @return cooling luminosity function of (t, T^inf) [GeV^{2}]
            /// @note Stefan-Boltzmann law
            std::function<double(double, double)> surface_luminosity(double R, double M, double eta);
        }

        /// @brief Neutrino related cooling functionality
        namespace neutrinic
        {
            /// @brief Emissivity of neutrinos from hadronic DUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param nbar_core_limit baryon density [GeV^3] at the core
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_n_1s0 allow/forbid superfluidity in 1S0 state for neutrons
            /// @param superfluid_p_1s0 allow/forbid superfluidity in 1S0 state for protons
            /// @param superfluid_n_3p2 allow/forbid superfluidity in 3P2 state for neutrons
            /// @param superfluid_p_temp temperature of superfluid protons [GeV]
            /// @param superfluid_n_temp temperature of superfluid neutrons [GeV]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV]
            /// @cite All - Yakovlev, Kaminker, 2000; AB joint superfluidity - Levenfish, Yakovlev, 1994
            std::function<double(double, const auxiliaries::phys::Species &, double, double)> hadron_durca_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from hadronic MUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param nbar_core_limit baryon density [GeV^3] at the core
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_n_1s0 allow/forbid superfluidity in 1S0 state for neutrons
            /// @param superfluid_p_1s0 allow/forbid superfluidity in 1S0 state for protons
            /// @param superfluid_n_3p2 allow/forbid superfluidity in 3P2 state for neutrons
            /// @param superfluid_p_temp temperature of superfluid protons [GeV]
            /// @param superfluid_n_temp temperature of superfluid neutrons [GeV]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV]
            /// @cite All - Yakovlev, Kaminker, 2000; Pion exchange coefficients - Yanagi, 2020
            std::function<double(double, const auxiliaries::phys::Species &, double, double)> hadron_murca_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from hadronic bremsstrahlung reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param ion_volume_frac volume fraction of ions as a function of baryon density [GeV^3]
            /// @param nbar_core_limit baryon density [GeV^3] at the core
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_n_1s0 allow/forbid superfluidity in 1S0 state for neutrons
            /// @param superfluid_p_1s0 allow/forbid superfluidity in 1S0 state for protons
            /// @param superfluid_n_3p2 allow/forbid superfluidity in 3P2 state for neutrons
            /// @param superfluid_p_temp temperature of superfluid protons [GeV]
            /// @param superfluid_n_temp temperature of superfluid neutrons [GeV]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV]
            /// @cite Yakovlev, Kaminker, 2000
            std::function<double(double, double, double)> hadron_bremsstrahlung_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &ion_volume_frac, double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from hadronic pair breaking & formations
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param nbar_core_limit baryon density [GeV^3] at the core
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_n_1s0 allow/forbid superfluidity in 1S0 state for neutrons
            /// @param superfluid_p_1s0 allow/forbid superfluidity in 1S0 state for protons
            /// @param superfluid_n_3p2 allow/forbid superfluidity in 3P2 state for neutrons
            /// @param superfluid_p_temp temperature of superfluid protons [GeV]
            /// @param superfluid_n_temp temperature of superfluid neutrons [GeV]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV]
            /// @cite Density expression - Page, 2009; superfluidity factors - Yakovlev, Kaminker, 1998
            std::function<double(double, const auxiliaries::phys::Species &, double, double)> hadron_pbf_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from quark DUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV]
            /// @cite Base density - Iwamoto, 1982; superconductivity effect - Blaschke, Grigorian 2001
            std::function<double(double, double, double)> quark_durca_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Emissivity of neutrinos from quark MUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV]
            /// @cite Base density - Iwamoto, 1982; superconductivity effect - Blaschke, Grigorian 2001
            std::function<double(double, double, double)> quark_murca_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Emissivity of neutrinos from quark Bremsstrahlung reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV]
            /// @cite Base density - Iwamoto, 1982; superconductivity effect - Blaschke, Grigorian 2001
            std::function<double(double, double, double)> quark_bremsstrahlung_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Emissivity of neutrinos from electron bremsstrahlung reaction (important in highly supressed quark phase)
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV]
            /// @cite All - Blaschke, Grigorian 2001
            std::function<double(double, double, double)> electron_bremsstrahlung_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi);
        }
    }
}

#endif // COOLING_H