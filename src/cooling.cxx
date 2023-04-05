
#include "../include/auxiliaries.h"
#include "../include/cooling.h"
#include "../include/constants.h"

#include <vector>
#include <functional>
#include <map>
#include <string>
#include <cmath>
#include <stdexcept>

double cooling::solver::stationary_cooling_cached(std::vector<std::vector<double>> &cache, double t, const std::function<double(double, double)> &cooling_rhs, double initial_temperature, double base_time_step, double exp_rate,
                                                  const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator)
{
    // times must be positive
    if (t < 0 || base_time_step <= 0)
        throw std::invalid_argument("Evolution time must be positive");
    // inverse euler solver; I want this method to be stable for any time step, including huge ones
    auto back_euler_step = [&cooling_rhs](double t, double T, double time_step)
    {
        // solve T_{n+1} - T_n - dt * a^n * F(t_{n+1}, T_{n+1}) = 0 with Newton's steps
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

    // if we conduct a first calculation (or time exceeded existing cache), we need to fill the cache
    if (cache.empty() || cache[0].back() < t)
    {
        // create two empty vectors
        cache = std::vector<std::vector<double>>(2, std::vector<double>());
        // fill the cache[0] with time exp. growing values, and the cache[1] with the corresponding temperatures
        cache[0].push_back(0.0);
        cache[1].push_back(initial_temperature);
        double t_curr = 0.0, time_step = base_time_step;
        do
        {
            cache[1].push_back(back_euler_step(t_curr, cache[1].back(), time_step));
            t_curr += time_step;
            cache[0].push_back(t_curr);
            time_step *= exp_rate;
        } while (t > t_curr);
    }
    // now we're sure the time is in the cache, we just interpolate
    return interpolator(cache[0], cache[1], t);
}

double cooling::predefined::auxiliary::te_tb_relation(double Tb, double R, double M, double eta)
{
    using namespace constants::scientific;
    using namespace constants::conversion;

    // calculate exp phi(R)
    double exp_phi_at_R = pow(1 - 2 * G * M / R, 0.5);

    // calculate g_14
    double g = G * M / (R * R * exp_phi_at_R);                    // gravitational acceleration at the surface in GeV
    double g_14 = g * gev_s * gev_s / (1.0E14 * km_gev * 1.0E-5); // normalized to 1E14 cm/s^2

    // local temperature on surface normalized to 1E9 K
    double T_local_9 = Tb * 1.0 / exp_phi_at_R * gev_over_k / 1.0E9;

    // now we just use the formulas from the paper
    double dzeta = T_local_9 - pow(7 * T_local_9 * pow(g_14, 0.5), 0.5) * 1.0E-3;
    double T_s6_Fe_to_4 = g_14 * (pow(7 * dzeta, 2.25) + pow(dzeta / 3, 1.25)),
           T_s6_a_to_4 = g_14 * pow(18.1 * T_local_9, 2.42);

    double a = (1.2 + pow(5.3 * 1.0E-6 / eta, 0.38)) * pow(T_local_9, 5.0 / 3);

    // surface temperature normalized to 1E6 K in 4th power
    double T_s6_to_4 = (a * T_s6_Fe_to_4 + T_s6_a_to_4) / (a + 1);

    return pow(T_s6_to_4, 1.0 / 4) * 1.0E6 / gev_over_k;
}

std::function<double(double, double)> cooling::predefined::photonic::surface_luminosity(double R, double M, double eta)
{
    return [=](double t, double T)
    {
        using namespace constants::scientific;
        double exp_2phi_at_R = 1 - 2 * G * M / R;
        return 4 * Pi * R * R * Sigma * pow(cooling::predefined::auxiliary::te_tb_relation(T, R, M, eta), 4) * exp_2phi_at_R;
    };
}

std::function<double(double, double)> cooling::predefined::auxiliary::fermi_specific_heat_cached(std::vector<double> &cache, const std::map<std::string, std::function<double(double)>> &m_star_functions, const std::map<std::string, std::function<double(double)>> &k_fermi_functions, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_lambda_of_r, const std::function<double(double)> &exp_phi_of_r, double r_ns, double radius_step)
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
            double part_integrand = cv_over_T(r);
            cache[0] += part_integrand * jacob * radius_step;
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
            double pf_l = k_fermi_l(nbar),
                   pf_n = k_fermi_n(nbar),
                   pf_p = k_fermi_p(nbar);
            double mst_n = m_star_n(nbar),
                   mst_p = m_star_p(nbar),
                   mst_l = m_star_l(nbar);
            if (pf_l + pf_p - pf_n <= 0)
                return 0.0;
            double dens = (4.24E27 / 1.68E54) * (mst_n / M_N) * (mst_p / M_N) * mst_l *
                          pow(exp_phi_of_r(r), -6) * pow(gev_over_k, 6) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);
            return dens;
        };

        // calculate the integral
        for (double r = radius_step; r < r_ns; r += radius_step)
        {
            double jacob = 4 * Pi * r * r * exp_lambda_of_r(r) * exp_phi_of_r(r) * exp_phi_of_r(r);
            double part_integrand = qv_over_t6(r);
            cache[0] += part_integrand * jacob * radius_step;
        }
    }

    return [&cache](double t, double T)
    {
        return cache[0] * pow(T, 6);
    };
}

std::function<double(double, double)> cooling::predefined::neutrinic::hadron_murca_luminocity_cached(std::vector<double> &cache, const std::function<double(double)> &m_star_n, const std::function<double(double)> &m_star_p, const std::function<double(double)> &m_star_l, const std::function<double(double)> &k_fermi_n, const std::function<double(double)> &k_fermi_p, const std::function<double(double)> &k_fermi_l, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_lambda_of_r, const std::function<double(double)> &exp_phi_of_r, double r_ns, double radius_step, double nbar_conversion)
{
    if (cache.empty())
    {
        cache = std::vector<double>(1, 0);
        using namespace constants::scientific;
        using namespace constants::conversion;

        // Qv/T^8inf = f(r)
        auto qv_over_t8 = [=](double r)
        {
            double nbar = nbar_of_r(r);
            double pf_l = k_fermi_l(nbar),
                   pf_n = k_fermi_n(nbar),
                   pf_p = k_fermi_p(nbar);
            double mst_n = m_star_n(nbar),
                   mst_p = m_star_p(nbar),
                   mst_l = m_star_l(nbar);
            double alpha = 1.76 - 0.63 * pow(N_sat / (nbar * nbar_conversion), 2.0 / 3), beta = 0.68,
                   v_fl = pf_l / mst_l;
            double dens = (8.05E21 / 1.68E72) * v_fl * pow(mst_n / M_N, 3) * (mst_p / M_N) * pf_p *
                          pow(exp_phi_of_r(r), -8) * alpha * beta *
                          pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);
            if (pf_l + 3 * pf_p - pf_n > 0)
                dens += (8.05E21 / (8 * 1.68E72)) * (pow(pf_l + 3 * pf_p - pf_n, 2) / mst_l) * pow(mst_p / M_N, 3) * (mst_n / M_N) *
                        pow(exp_phi_of_r(r), -8) * alpha * beta *
                        pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);
            return dens;
        };

        // calculate the integral
        for (double r = radius_step; r < r_ns; r += radius_step)
        {
            double jacob = 4 * Pi * r * r * exp_lambda_of_r(r) * exp_phi_of_r(r) * exp_phi_of_r(r);
            double part_integrand = qv_over_t8(r);
            cache[0] += part_integrand * jacob * radius_step;
        }
    }

    return [&cache](double t, double T)
    {
        return cache[0] * pow(T, 8);
    };
}

std::function<double(double, double)> cooling::predefined::neutrinic::hadron_bremsstrahlung_luminocity_cached(std::vector<double> &cache, const std::function<double(double)> &m_star_n, const std::function<double(double)> &m_star_p, const std::function<double(double)> &k_fermi_n, const std::function<double(double)> &k_fermi_p, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_lambda_of_r, const std::function<double(double)> &exp_phi_of_r, double r_ns, double radius_step)
{
    if (cache.empty())
    {
        cache = std::vector<double>(1, 0);
        using namespace constants::scientific;
        using namespace constants::conversion;

        // Qv/T^8inf = f(r)
        auto qv_over_t8 = [=](double r)
        {
            double nbar = nbar_of_r(r);
            double pf_n = k_fermi_n(nbar),
                   pf_p = k_fermi_p(nbar);
            double mst_n = m_star_n(nbar),
                   mst_p = m_star_p(nbar);
            double alpha_nn = 0.59, alpha_np = 1.06, alpha_pp = 0.11,
                   beta_nn = 0.56, beta_np = 0.66, beta_pp = 0.7;
            int n_flavours = 3;
            double dens_nn = (7.5E19 / 1.68E72) * pow(mst_n / M_N, 4) * pf_n * n_flavours *
                             pow(exp_phi_of_r(r), -8) * alpha_nn * beta_nn *
                             pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3),
                   dens_pp = (7.5E19 / 1.68E72) * pow(mst_p / M_N, 4) * pf_p * n_flavours *
                             pow(exp_phi_of_r(r), -8) * alpha_pp * beta_pp *
                             pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3),
                   dens_np = (1.5E20 / 1.68E72) * pow(mst_p / M_N, 2) * pow(mst_n / M_N, 2) * pf_p * n_flavours *
                             pow(exp_phi_of_r(r), -8) * alpha_np * beta_np *
                             pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);

            return dens_nn + dens_np + dens_pp;
        };

        // calculate the integral
        for (double r = radius_step; r < r_ns; r += radius_step)
        {
            double jacob = 4 * Pi * r * r * exp_lambda_of_r(r) * exp_phi_of_r(r) * exp_phi_of_r(r);
            double part_integrand = qv_over_t8(r);
            cache[0] += part_integrand * jacob * radius_step;
        }
    }

    return [&cache](double t, double T)
    {
        return cache[0] * pow(T, 8);
    };
}

double cooling::predefined::auxiliary::critical_temperature_smeared_guassian(double k_fermi, double temp_ampl, double k_offs, double k_width, double quad_skew)
{
    return temp_ampl * exp(-pow((k_fermi - k_offs) / k_width, 2) - quad_skew * pow((k_fermi - k_offs) / k_width, 4));
}

double cooling::predefined::auxiliary::critical_temperature_ao(double k_fermi)
{
    using namespace constants::conversion;
    return cooling::predefined::auxiliary::critical_temperature_smeared_guassian(
        k_fermi, 2.35E9 / gev_over_k, 0.49 / (1.0E-18 * km_gev), 0.31 / (1.0E-18 * km_gev), 0.0);
}

double cooling::predefined::auxiliary::critical_temperature_ccdk(double k_fermi)
{
    using namespace constants::conversion;
    return cooling::predefined::auxiliary::critical_temperature_smeared_guassian(
        k_fermi, 6.6E9 / gev_over_k, 0.66 / (1.0E-18 * km_gev), 0.46 / (1.0E-18 * km_gev), 0.69);
}

double cooling::predefined::auxiliary::critical_temperature_a(double k_fermi)
{
    using namespace constants::conversion;
    return cooling::predefined::auxiliary::critical_temperature_smeared_guassian(
        k_fermi, 1.0E9 / gev_over_k, 1.8 / (1.0E-18 * km_gev), 0.5 / (1.0E-18 * km_gev), 0.0);
}

double cooling::predefined::auxiliary::critical_temperature_b(double k_fermi)
{
    using namespace constants::conversion;
    return cooling::predefined::auxiliary::critical_temperature_smeared_guassian(
        k_fermi, 3.0E9 / gev_over_k, 2.0 / (1.0E-18 * km_gev), 0.5 / (1.0E-18 * km_gev), 0.0);
}

double cooling::predefined::auxiliary::critical_temperature_c(double k_fermi)
{
    using namespace constants::conversion;
    return cooling::predefined::auxiliary::critical_temperature_smeared_guassian(
        k_fermi, 1.0E10 / gev_over_k, 2.5 / (1.0E-18 * km_gev), 0.7 / (1.0E-18 * km_gev), 0.0);
}

double cooling::predefined::auxiliary::critical_temperature_a2(double k_fermi)
{
    using namespace constants::conversion;
    return cooling::predefined::auxiliary::critical_temperature_smeared_guassian(
        k_fermi, 5.5E9 / gev_over_k, 2.3 / (1.0E-18 * km_gev), 0.9 / (1.0E-18 * km_gev), 0.0);
}