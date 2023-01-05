
#include "../include/auxiliaries.h"
#include "../include/cooling.h"
#include "../include/constants.h"

#include <vector>
#include <functional>
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
            F_shift = cooling_rhs(t + time_step, T_new + time_step * F) - F;
            T_new -= (T_new - T - time_step * F) / (1 - F_shift / F);
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