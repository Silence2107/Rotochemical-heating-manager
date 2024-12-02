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
        /// @brief Solves the equilibrium cooling equation dT^inf/dt = F(t,T^inf) for a rhs function and a given initial temperature (inf means in distant frame)
        /// @param t_curr time at which the initial temperature profile was calculated [GeV^{-1}]
        /// @param t_step time step at which the next profiles are to be calculated [GeV^{-1}]
        /// @param cooling_rhs right hand side term F(t,T^inf) of the cooling equation, i.e. L^inf/Cv^inf [GeV^{2}]
        /// @param initial_temperature Temperature^inf at t=t_curr [GeV]
        /// @param newton_eps desirable relative accuracy of the solution
        /// @param newton_iter_max maximum number of iterations for the Newton-Raphson method
        /// @return [0] : Temperature^inf at t_curr + t_step [GeV]; [1] : boolean flag for reaching the adaption limit; [2] - boolean flag for negative temperature overshoot
        std::vector<double> equilibrium_cooling(
            double t_curr, double t_step, const std::function<double(double, double)> &cooling_rhs, double initial_temperature, double newton_eps, size_t newton_iter_max);

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
        /// @param newton_eps desirable relative accuracy of the solution
        /// @param newton_iter_max maximum number of iterations for the Newton-Raphson method
        /// @return Temperature, luminosity profiles [T^inf(r, t), Ld^inf(r, t)] as an [2][i_m + 1] array corresponding to radii; [2] control info : [0] - boolean flag for reaching the adaption limit, [1] - boolean flag for negative temperature overshoot
        std::vector<std::vector<double>> nonequilibrium_cooling(
            double t_curr, double t_step, const std::function<double(double, double, double)> &neutrino_rate, const std::function<double(double, double, double)> &cv, const std::function<double(double, double, double)> &lambda,
            const std::function<double(double)> &exp_lambda, const std::function<double(double)> &exp_phi, const std::vector<double> &radii, const std::vector<double> &initial_profile,
            const std::function<double(double)> &te_tb, double newton_eps, size_t newton_iter_max);

        /// @brief Solves the coupled cooling equation dXi(t)/dt = Fi(t,{Xi}). These depend on the context, e.g. they may reflect temperature evolution alongside chemical imbalances. 0th component is ALWAYS temperature (and ensured positive)
        /// @param t_curr time at which the initial values were calculated [GeV^{-1}]
        /// @param t_step time step at which the next values are to be calculated [GeV^{-1}]
        /// @param rhs right hand vector side term Fi(t, {Xi}) of the cooling equation
        /// @param initial_values value vector at t=t_curr [GeV]
        /// @param newton_eps desirable relative accuracy of the solution
        /// @param newton_iter_max maximum number of iterations for the Newton-Raphson method
        /// @return [0] : Evolved values at t_curr + t_step [GeV]; [1] control info : [0] - boolean flag for reaching the adaption limit, [1] - boolean flag for negative temperature overshoot
        std::vector<std::vector<double>> coupled_cooling(
            double t_curr, double t_step, const std::function<std::vector<double>(double, const std::vector<double> &)> &rhs,
            const std::vector<double> &initial_values, double newton_eps, size_t newton_iter_max);
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
            /// @param nbar_sf_shift lowest baryon density [GeV^3] with triplet pairing
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_p_temp proton SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @param superfluid_n_temp neutron SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV]
            /// @cite All - Yakovlev, Kaminker, 2000; AB joint superfluidity - Levenfish, Yakovlev, 1994
            std::function<double(double, const auxiliaries::phys::Species &, double, double)> hadron_durca_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_sf_shift, const std::function<double(double)> &exp_phi,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from hadronic MUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param nbar_sf_shift lowest baryon density [GeV^3] with triplet pairing
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_p_temp proton SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @param superfluid_n_temp neutron SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV]
            /// @cite All - Yakovlev, Kaminker, 2000; Pion exchange coefficients - Yanagi, 2020
            std::function<double(double, const auxiliaries::phys::Species &, double, double)> hadron_murca_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_sf_shift, const std::function<double(double)> &exp_phi,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from hadronic bremsstrahlung reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param ion_volume_frac volume fraction of ions as a function of baryon density [GeV^3]
            /// @param nbar_sf_shift lowest baryon density [GeV^3] with triplet pairing
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_p_temp proton SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @param superfluid_n_temp neutron SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV]
            /// @cite Yakovlev, Kaminker, 2000
            std::function<double(double, double, double)> hadron_bremsstrahlung_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &ion_volume_frac, double nbar_sf_shift, const std::function<double(double)> &exp_phi,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from hadronic pair breaking & formations
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param nbar_sf_shift lowest baryon density [GeV^3] with triplet pairing
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_p_temp proton SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @param superfluid_n_temp neutron SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV]
            /// @cite Density expression - Page, 2009; superfluidity factors - Yakovlev, Kaminker, 1998
            std::function<double(double, const auxiliaries::phys::Species &, double, double)> hadron_pbf_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_sf_shift, const std::function<double(double)> &exp_phi,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from up-down quark DUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV]
            /// @cite Base density - Iwamoto, 1982; superconductivity effect - Blaschke, Grigorian 2001
            std::function<double(double, double, double)> quark_ud_durca_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Emissivity of neutrinos from up-strange quark DUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV]
            /// @cite Base density - Iwamoto, 1982; superconductivity effect - Blaschke, Grigorian 2001
            std::function<double(double, double, double)> quark_us_durca_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Emissivity of neutrinos from up-down quark MUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV]
            /// @cite Base density - Iwamoto, 1982; superconductivity effect - Blaschke, Grigorian 2001
            std::function<double(double, double, double)> quark_ud_murca_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Emissivity of neutrinos from up-strange quark MUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV]
            /// @cite Base density - Iwamoto, 1982; superconductivity effect - Blaschke, Grigorian 2001
            std::function<double(double, double, double)> quark_us_murca_emissivity(
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

        /// @brief Rotochemical heating cooling related functionality
        namespace rotochemical
        {
            /// @brief Emissivity of neutrinos from hadronic DUrca reactions enhanced by rotochemical heating
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param nbar_sf_shift lowest baryon density [GeV^3] with triplet pairing
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_p_temp proton SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @param superfluid_n_temp neutron SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, const auxiliaries::phys::Species &, double, double, double)> hadron_durca_enhanced_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_sf_shift, const std::function<double(double)> &exp_phi,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from hadronic MUrca reactions enhanced by rotochemical heating
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param nbar_sf_shift lowest baryon density [GeV^3] with triplet pairing
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_p_temp proton SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @param superfluid_n_temp neutron SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, const auxiliaries::phys::Species &, double, double, double)> hadron_murca_enhanced_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_sf_shift, const std::function<double(double)> &exp_phi,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Emissivity of neutrinos from up-down quark DUrca reactions enhanced by rotochemical heating
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, double, double, double)> quark_ud_durca_enhanced_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Emissivity of neutrinos from up-strange quark DUrca reactions enhanced by rotochemical heating
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, double, double, double)> quark_us_durca_enhanced_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Emissivity of neutrinos from quark MUrca (up-down) reactions enhanced by rotochemical heating
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, double, double, double)> quark_ud_murca_enhanced_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Emissivity of neutrinos from quark MUrca (up-strange) reactions enhanced by rotochemical heating
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return emissivity [GeV^5] as a function of radius [GeV^{-1}], time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, double, double, double)> quark_us_murca_enhanced_emissivity(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Rate difference of neutrinos from forward and backward hadronic DUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param nbar_sf_shift lowest baryon density [GeV^3] with triplet pairing
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_p_temp proton SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @param superfluid_n_temp neutron SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @return Rate difference [GeV^4] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, const auxiliaries::phys::Species &, double, double, double)> hadron_durca_rate_difference(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_sf_shift, const std::function<double(double)> &exp_phi,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Rate difference of neutrinos from forward and backward hadronic MUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param nbar_sf_shift lowest baryon density [GeV^3] with triplet pairing
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superfluid_p_temp proton SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @param superfluid_n_temp neutron SF critical temperature [GeV] as a function of baryon density [GeV^3]
            /// @return Rate difference [GeV^4] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, const auxiliaries::phys::Species &, double, double, double)> hadron_murca_rate_difference(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                double nbar_sf_shift, const std::function<double(double)> &exp_phi,
                const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp);

            /// @brief Rate difference of neutrinos from forward and backward up-down quark DUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return Rate difference [GeV^4] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, double, double, double)> quark_ud_durca_rate_difference(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Rate difference of neutrinos from forward and backward up-strange quark DUrca reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return Rate difference [GeV^4] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, double, double, double)> quark_us_durca_rate_difference(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Rate difference of neutrinos from forward and backward quark MUrca (up-down) reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return Rate difference [GeV^4] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, double, double, double)> quark_ud_murca_rate_difference(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);

            /// @brief Rate difference of neutrinos from forward and backward quark MUrca (up-strange) reactions
            /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
            /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
            /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
            /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
            /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
            /// @return Rate difference [GeV^4] as a function of radius [GeV^{-1}], species, time [GeV], temperature [GeV] and chemical imbalance [GeV]
            /// @cite RH control function - Reisenegger 1994 (arXiv:astro-ph/9410035)
            std::function<double(double, double, double, double)> quark_us_murca_rate_difference(
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi, const std::function<double(double)> &superconduct_q_gap);
        }
    }
}

#endif // COOLING_H