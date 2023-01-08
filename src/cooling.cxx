
#include "../include/auxiliaries.h"
#include "../include/cooling.h"
#include "../include/constants.h"

#include <vector>
#include <functional>
#include <map>
#include <string>
#include <cmath>
#include <stdexcept>

double cooling::solver::stationary_cooling_cached(std::vector<std::vector<double>> &cache, double t, const std::function<double(double, double)> &cooling_rhs, double initial_temperature, double time_step,
                                                  const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator)
{
    // times must be positive
    if (t < 0 || time_step <= 0)
        throw std::invalid_argument("Evolution time must be positive");
    // inverse euler solver; I want this method to be stable for any time step, including huge ones
    auto back_euler_step = [&cooling_rhs, time_step](double t, double T)
    {
        // solve T_{n+1} - T_n - dt * F(t_{n+1}, T_{n+1}) = 0 with Newton's steps
        double eps = 1e-10;
        size_t max_iter = 100, iter = 0;
        double T_new = T;
        double F, F_shift;
        do
        {
            ++iter;
            F = cooling_rhs(t + time_step, T_new);
            auto temp_step = time_step / (t + time_step) * T_new;
            F_shift = cooling_rhs(t + time_step, T_new + temp_step) - F;
            T_new -= (T_new - T - time_step * F) / (1 - time_step * F_shift / temp_step);
            if (iter > max_iter)
                break;
        } while (std::abs(T_new - T - time_step * F) > eps * T);
        return T_new;
    };

    // if we conduct a first calculation, we need to fill the cache
    if (cache.empty())
    {
        size_t discr_size = std::ceil(t / time_step) + 2;
        cache = std::vector<std::vector<double>>(2, std::vector<double>(discr_size, 0));
        // fill the cache[0] with time linspaced values, and the cache[1] with the corresponding temperatures
        cache[0][0] = 0;
        cache[1][0] = initial_temperature;
        for (size_t i = 1; i < discr_size; ++i)
        {
            cache[0][i] = cache[0][i - 1] + time_step;
            cache[1][i] = back_euler_step(cache[0][i], cache[1][i - 1]);
        }
    }
    // if the maximum time exceeds the one in the cache, we need to extend the cache
    else if (t > cache[0].back())
    {
        size_t add_discr_size = std::ceil((t - cache[0].back()) / time_step) + 2;
        cache[0].resize(cache[0].size() + add_discr_size);
        cache[1].resize(cache[1].size() + add_discr_size);
        for (size_t i = cache[0].size() - add_discr_size; i < cache[0].size(); ++i)
        {
            cache[0][i] = cache[0][i - 1] + time_step;
            cache[1][i] = back_euler_step(cache[0][i], cache[1][i - 1]);
        }
    }
    // now we're sure the time is in the cache, we just interpolate
    return interpolator(cache[0], cache[1], t);
}

std::function<double(double, double)> cooling::predefined::photonic::surface_luminosity(double R, double M, double eta)
{
    return [=](double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::conversion;

        // calculate exp phi(R)
        double exp_phi_at_R = pow(1 - 2 * G * M / R, 0.5);

        // calculate g_14
        double g = G * M / (R * R * exp_phi_at_R);                    // gravitational acceleration at the surface in GeV
        double g_14 = g * gev_s * gev_s / (1.0E14 * km_gev * 1.0E-5); // normalized to 1E14 cm/s^2

        // local temperature on surface normalized to 1E9 K
        double T_local_9 = T * 1.0 / exp_phi_at_R * gev_over_k / 1.0E9;

        // now we just use the formulas from the paper
        double dzeta = T_local_9 - pow(7 * T_local_9 * pow(g_14, 0.5), 0.5) * 1.0E-3;
        double T_s6_Fe_to_4 = g_14 * (pow(7 * dzeta, 2.25) + pow(dzeta / 3, 1.25)),
               T_s6_a_to_4 = g_14 * pow(18.1 * T_local_9, 2.42);

        double a = (1.2 + pow(5.3 * 1.0E-6 / eta, 0.38)) * pow(T_local_9, 5.0 / 3);

        // surface temperature normalized to 1E6 K in 4th power
        double T_s6_to_4 = (a * T_s6_Fe_to_4 + T_s6_a_to_4) / (a + 1);

        return 4 * Pi * R * R * Sigma * T_s6_to_4 * pow(1.0E6/gev_over_k, 4) * pow(exp_phi_at_R, 2);
    };
}

std::function<double(double, double)> cooling::predefined::specific_heat::fermi_specific_heat_cached(std::vector<double> &cache, const std::map<std::string, std::function<double(double)>> &m_star_functions, const std::map<std::string, std::function<double(double)>> &k_fermi_functions, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_lambda_of_r, const std::function<double(double)> &exp_phi_of_r, double r_ns, double radius_step)
{
    if (cache.empty())
    {
        cache = std::vector<double>(1, 0);
        using namespace constants::scientific;

        // Cv/T^inf density
        auto cv_over_T = [=](double r)
        {
            double cv_over_T_dens = 0;
            for (auto it = m_star_functions.begin(); it != m_star_functions.end(); ++it)
            {
                auto key = it->first;
                double nbar = nbar_of_r(r);
                double m_star = m_star_functions.at(key)(nbar);
                double k_fermi = k_fermi_functions.at(key)(nbar);
                double exp_min_phi = 1.0 / exp_phi_of_r(r);
                cv_over_T_dens += m_star * k_fermi / 3.0 * exp_min_phi;
            }
            return cv_over_T_dens;
        };

        // calculate the integral
        for (double r = radius_step; r < r_ns; r += radius_step)
        {
            double jacob = 4 * Pi * r * r * exp_lambda_of_r(r);
            double cv_over_T_dens = cv_over_T(r);
            cache[0] += cv_over_T_dens * jacob * radius_step;
        }
    }
    return [&cache](double t, double T)
    {
        return cache[0] * T;
    };
}

std::function<double(double, double)> cooling::predefined::neutrinic::hadron_durca_luminocity_cached(std::vector<double> &cache, const std::function<double(double)> &m_star_n, const std::function<double(double)> &m_star_p, const std::function<double(double)> &m_star_l, const std::function<double(double)> &k_fermi_n, const std::function<double(double)> &k_fermi_p, const std::function<double(double)> &k_fermi_l, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_lambda_of_r, const std::function<double(double)> &exp_phi_of_r, double r_ns, double radius_step)
{
    if (cache.empty())
    {
        cache = std::vector<double>(1, 0);
        using namespace constants::scientific;
        using namespace constants::conversion;

        // Qv/T^6inf = f(r)
        auto qv_over_t6 = [=](double r)
        {
            double nbar = nbar_of_r(r);
            if (k_fermi_l(nbar) + k_fermi_p(nbar) - k_fermi_n(nbar) <= 0)
                return 0.0;
            double dens = (4.001E27 / 1.68E54) * (m_star_n(nbar) / M_N) * (m_star_p(nbar) / M_N) * m_star_l(nbar) *
                          pow(exp_phi_of_r(r), -6) * pow(gev_over_k, 6) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);
            return dens;
        };

        // calculate the integral
        for (double r = radius_step; r < r_ns; r += radius_step)
        {
            double jacob = 4 * Pi * r * r * exp_lambda_of_r(r) * exp_phi_of_r(r) * exp_phi_of_r(r);
            double cv_over_T_dens = qv_over_t6(r);
            cache[0] += cv_over_T_dens * jacob * radius_step;
        }
    }

    return [&cache](double t, double T)
    {
        return cache[0] * pow(T, 6);
    };
}