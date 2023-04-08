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
        
            /// @brief Critical temperature in AO model
            /// @param k_fermi fermi momentum [GeV], correspoding to the nucleon
            double critical_temperature_ao(double k_fermi);

            /// @brief Critical temperature in CCDK model
            /// @param k_fermi fermi momentum [GeV], correspoding to the nucleon
            double critical_temperature_ccdk(double k_fermi);

            /// @brief Critical temperature in "a" model
            /// @param k_fermi fermi momentum [GeV], correspoding to the nucleon
            double critical_temperature_a(double k_fermi);

            /// @brief Critical temperature in "b" model
            /// @param k_fermi fermi momentum [GeV], correspoding to the nucleon
            double critical_temperature_b(double k_fermi);

            /// @brief Critical temperature in "c" model
            /// @param k_fermi fermi momentum [GeV], correspoding to the nucleon
            double critical_temperature_c(double k_fermi);

            /// @brief Critical temperature in "a2" model
            /// @param k_fermi fermi momentum [GeV], correspoding to the nucleon
            double critical_temperature_a2(double k_fermi);
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
        {}
    }
}

#endif // COOLING_H