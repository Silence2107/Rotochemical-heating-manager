
#include "../include/auxiliaries.h"
#include "../include/cooling.h"
#include "../include/constants.h"

#include <vector>
#include <functional>
#include <map>
#include <string>
#include <cmath>
#include <stdexcept>
// #include <iostream>

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
        double eps = 1e-5;
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
            if (T_new < 0)
                throw std::runtime_error("Reached negative temperature with current method; Encountered in stationary_cooling_cached");
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
            // std::cout << cache[0].back() << " " << cache[1].back() << '\n';
        } while (t > t_curr);
    }
    // now we're sure the time is in the cache, we just interpolate
    return interpolator(cache[0], cache[1], t);
}

std::function<double(double, double, double)> cooling::predefined::auxiliary::fermi_specific_heat_density(
    const std::map<auxiliaries::Species, std::function<double(double)>> &k_fermi_of_nbar,
    const std::map<auxiliaries::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
    double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
    bool superconduct_u, bool superconduct_d, bool superconduct_s, const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp,
    const std::function<double(double)> &superconduct_u_temp, const std::function<double(double)> &superconduct_d_temp, const std::function<double(double)> &superconduct_s_temp)
{
    return [=](double r, double t, double T)
    {
        using namespace constants::scientific;

        double cv_dens = 0;
        for (auto it = m_stars_of_nbar.begin(); it != m_stars_of_nbar.end(); ++it)
        {
            auto key = it->first;
            double nbar_val = nbar_of_r(r);
            double m_star = m_stars_of_nbar.at(key)(nbar_val);
            double k_fermi = k_fermi_of_nbar.at(key)(nbar_val);
            double T_loc = T / exp_phi(r);
            double diff = m_star * k_fermi / 3.0 * T_loc;

            // superfluid factors
            auto r_A = [&](double v)
            {
                return pow(0.4186 + sqrt(pow(1.007, 2.0) + pow(0.5010 * v, 2.0)), 2.5) *
                       exp(1.456 - sqrt(pow(1.456, 2.0) + pow(v, 2.0)));
            };
            auto r_B = [&](double v)
            {
                return pow(0.6893 + sqrt(pow(0.790, 2.0) + pow(0.2824 * v, 2.0)), 2.0) *
                       exp(1.934 - sqrt(pow(1.934, 2.0) + pow(v, 2.0)));
            };

            // proton superfluidity?
            if (key == proton && superfluid_p_1s0)
            {
                double tau = T_loc / superfluid_p_temp(k_fermi);
                if (tau < 1.0)
                {
                    using namespace cooling::predefined::auxiliary;
                    diff *= r_A(superfluid_gap_1s0(tau));
                }
            }
            // neutron superfluidity?
            else if (key == neutron && (superfluid_n_3p2 || superfluid_n_1s0))
            {
                double tau = T_loc / superfluid_p_temp(k_fermi);
                if (tau < 1.0)
                {
                    using namespace cooling::predefined::auxiliary;
                    // 1S0/3P2 division at core entrance
                    if (superfluid_n_1s0 && superfluid_n_3p2)
                        diff *= (nbar_val > nbar_core_limit) ? r_B(superfluid_gap_3p2(tau)) : r_A(superfluid_gap_1s0(tau));
                    // 3P2 only
                    else if (superfluid_n_3p2)
                        diff *= r_B(superfluid_gap_3p2(tau));
                    // 1S0 only
                    else
                        diff *= r_A(superfluid_gap_1s0(tau));
                }
            }
            // u quark superconductivity?
            else if (key == uquark && superconduct_u)
            {
                double tau = T_loc / superconduct_u_temp(k_fermi);
                if (tau < 1.0)
                {
                    using namespace cooling::predefined::auxiliary;
                    diff *= 3.1 / pow(tau, 2.5) * exp(-1.76 / tau * sqrt(1.0 - tau));
                }
            }
            // d quark superconductivity?
            else if (key == dquark && superconduct_d)
            {
                double tau = T_loc / superconduct_d_temp(k_fermi);
                if (tau < 1.0)
                {
                    using namespace cooling::predefined::auxiliary;
                    diff *= 3.1 / pow(tau, 2.5) * exp(-1.76 / tau * sqrt(1.0 - tau));
                }
            }
            // s quark superconductivity?
            else if (key == squark && superconduct_s)
            {
                double tau = T_loc / superconduct_s_temp(k_fermi);
                if (tau < 1.0)
                {
                    using namespace cooling::predefined::auxiliary;
                    diff *= 3.1 / pow(tau, 2.5) * exp(-1.76 / tau * sqrt(1.0 - tau));
                }
            }
            cv_dens += diff;
        }
        return cv_dens;
    };
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
    // dzeta > 0 -- from rotochemical notes
    double T_s6_Fe_to_4 = (dzeta > 0 ? g_14 * (pow(7 * dzeta, 2.25) + pow(dzeta / 3, 1.25)) : 0.0),
           T_s6_a_to_4 = g_14 * pow(18.1 * T_local_9, 2.42);

    double a = (1.2 + pow(5.3 * 1.0E-6 / eta, 0.38)) * pow(T_local_9, 5.0 / 3);

    // surface temperature normalized to 1E6 K in 4th power
    double T_s6_to_4 = (a * T_s6_Fe_to_4 + T_s6_a_to_4) / (a + 1);
    double Ts = pow(T_s6_to_4, 1.0 / 4) * 1.0E6 / gev_over_k;
    return (Ts < Tb / exp_phi_at_R) ? Ts : Tb / exp_phi_at_R;
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

double cooling::predefined::auxiliary::critical_temperature_smeared_guassian(double k_fermi, double temp_ampl, double k_offs, double k_width, double quad_skew)
{
    return temp_ampl * exp(-pow((k_fermi - k_offs) / k_width, 2) - quad_skew * pow((k_fermi - k_offs) / k_width, 4));
}

double cooling::predefined::auxiliary::superfluid_gap_1s0(double tau)
{
    return std::sqrt(1 - tau) * (1.456 - 0.157 / std::sqrt(tau) + 1.764 / tau);
}

double cooling::predefined::auxiliary::superfluid_gap_3p2(double tau)
{
    return std::sqrt(1 - tau) * (0.7893 + 1.188 / tau);
}

double cooling::predefined::auxiliary::critical_temperature(double k_fermi, cooling::predefined::auxiliary::CriticalTemperatureModel model)
{
    using namespace constants::conversion;
    using namespace cooling::predefined::auxiliary;
    switch (model)
    {
    case CriticalTemperatureModel::kA:
        return critical_temperature_smeared_guassian(
            k_fermi, 1.0E9 / gev_over_k, 1.8 / (1.0E-18 * km_gev), 0.5 / (1.0E-18 * km_gev), 0.0);
    case CriticalTemperatureModel::kB:
        return critical_temperature_smeared_guassian(
            k_fermi, 3.0E9 / gev_over_k, 2.0 / (1.0E-18 * km_gev), 0.5 / (1.0E-18 * km_gev), 0.0);
    case CriticalTemperatureModel::kC:
        return critical_temperature_smeared_guassian(
            k_fermi, 1.0E10 / gev_over_k, 2.5 / (1.0E-18 * km_gev), 0.7 / (1.0E-18 * km_gev), 0.0);
    case CriticalTemperatureModel::kA2:
        return critical_temperature_smeared_guassian(
            k_fermi, 5.5E9 / gev_over_k, 2.3 / (1.0E-18 * km_gev), 0.9 / (1.0E-18 * km_gev), 0.0);
    case CriticalTemperatureModel::kCCDK:
        return critical_temperature_smeared_guassian(
            k_fermi, 6.6E9 / gev_over_k, 0.66 / (1.0E-18 * km_gev), 0.46 / (1.0E-18 * km_gev), 0.69);
    case CriticalTemperatureModel::kAO:
        return critical_temperature_smeared_guassian(
            k_fermi, 2.35E9 / gev_over_k, 0.49 / (1.0E-18 * km_gev), 0.31 / (1.0E-18 * km_gev), 0.0);
    default:
        throw std::runtime_error("Unknown critical temperature model");
    }
}

std::function<double(double, const auxiliaries::Species &, double, double)> cooling::predefined::neutrinic::hadron_durca_emissivity(
    const std::map<auxiliaries::Species, std::function<double(double)>> &k_fermi_of_nbar,
    const std::map<auxiliaries::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
    double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
    const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp)
{
    return [=](double r, const auxiliaries::Species &lepton_flavour, double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::conversion;

        double nbar_val = nbar_of_r(r);
        double pf_l = k_fermi_of_nbar.at(lepton_flavour)(nbar_val),
               pf_n = k_fermi_of_nbar.at(neutron)(nbar_val),
               pf_p = k_fermi_of_nbar.at(proton)(nbar_val);
        double mst_n = m_stars_of_nbar.at(neutron)(nbar_val),
               mst_p = m_stars_of_nbar.at(proton)(nbar_val),
               mst_l = m_stars_of_nbar.at(lepton_flavour)(nbar_val);
        double T_loc = T / exp_phi(r);
        if (pf_l + pf_p - pf_n <= 0)
            return 0.0;
        double dens = (4.001E27 / 1.68E54) * (mst_n / M_N) * (mst_p / M_N) * mst_l *
                      pow(T_loc, 6) * pow(gev_over_k, 6) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);

        // Superfluid factors

        auto r_A = [](double v)
        {
            return pow(0.2312 + sqrt(pow(0.7688, 2.0) + pow(0.1438 * v, 2.0)), 5.5) *
                   exp(3.427 - sqrt(pow(3.427, 2.0) + pow(v, 2.0)));
        };
        auto r_B = [](double v)
        {
            return pow(0.2546 + sqrt(pow(0.7454, 2.0) + pow(0.1284 * v, 2.0)), 5) *
                   exp(2.701 - sqrt(pow(2.701, 2.0) + pow(v, 2.0)));
        };

        // I calculate inverse taus so that I do not hit 1/0, and also it signals about whether superfluidity is present
        double tau_p_inv = superfluid_p_temp(pf_p) / T_loc,
               tau_n_inv = superfluid_n_temp(pf_n) / T_loc;

        // normal fluidity
        if (tau_p_inv <= 1.0 && tau_n_inv <= 1.0)
        {
            return dens;
        }
        // pure proton superfluidity?
        else if (tau_p_inv > 1.0 && tau_n_inv <= 1.0)
        {
            using namespace cooling::predefined::auxiliary;
            return dens * r_A(superfluid_gap_1s0(1 / tau_p_inv));
        }
        // pure neutron superfluidity?
        else if (tau_p_inv <= 1.0 && tau_n_inv > 1.0)
        {
            using namespace cooling::predefined::auxiliary;
            //  1S0/3P2 division at core entrance
            if (superfluid_n_3p2 && superfluid_n_1s0)
            {
                return dens * (nbar_val > nbar_core_limit) ? r_B(superfluid_gap_3p2(1 / tau_n_inv)) : r_A(superfluid_gap_1s0(1 / tau_n_inv));
            }
            // 1S0 only
            else if (superfluid_n_1s0)
            {
                return dens * r_A(superfluid_gap_1s0(1 / tau_n_inv));
            }
            // 3P2 only
            else
            {
                return dens * r_B(superfluid_gap_3p2(1 / tau_n_inv));
            }
        }
        // combined superfluidity
        else
        {
            // Here we have quite cumbersome factors, so I do not define them in outer scope

            // for n1S0 neutron superfluidity
            auto r_AA = [](double v1, double v2)
            {
                using constants::scientific::Pi;
                double u = pow(v1, 2.0) + pow(v2, 2.0),
                       w = pow(v2, 2.0) - pow(v1, 2.0);
                double u1 = 1.8091 + sqrt(pow(v1, 2.0) + pow(2.2476, 2.0)),
                       u2 = 1.8091 + sqrt(pow(v2, 2.0) + pow(2.2476, 2.0));
                double D = 1.52 * pow(u1 * u2, 1.5) * (pow(u1, 2.0) + pow(u2, 2.0)) * exp(-u1 - u2);
                double pe = 0.5 * (u + 0.43847 + sqrt(pow(w, 2.0) + 8.3680 * u + 491.32)),
                       ps = 0.5 * (u + sqrt(pow(w, 2.0) + 5524.8 * u + 6.7737)),
                       q = 0.5 * (u + 12.421 - sqrt(pow(w, 2.0) + 16.350 * u + 45.171)),
                       p = 0.5 * (u + 12.421 + sqrt(pow(w, 2.0) + 16.350 * u + 45.171));
                double K2 = 7.0 * pow(Pi, 4.0) / 60 * sqrt(p - q),
                       K1 = pow(Pi, 2.0) * (1.0 / 6 * sqrt(p - q) * (p + 2 * q) - 1.0 / 2 * q * sqrt(p) * log((sqrt(p) + sqrt(p - q)) / sqrt(q))),
                       K0 = sqrt(p - q) / 120 * (6 * p * p + 83 * p * q + 16 * q * q) - sqrt(p) * q / 8 * (4 * p + 3 * q) * log((sqrt(p) + sqrt(p - q)) / sqrt(q));
                double S = 5040.0 / 457 * pow(Pi, 6.0) * (K0 + K1 + 0.42232 * K2) * sqrt(Pi * sqrt(ps) / 2) * exp(-sqrt(pe));
                return D + S * u / (u + 0.9163);
            };

            // for n3P2 neutron superfluidity
            auto r_AB = [](double tau_n, double tau_p)
            {
                size_t index_1 = ceil(-10 * log10(tau_p) - 1),
                       index_2 = ceil(-10 * log10(tau_n) - 1);
                if (index_1 > 15 || index_2 > 15)
                    return 0.0;

                // No analytic formula for this case, so I use a (lg tau_p, lg tau_n) \in
                // [-0.1, -0.2 ,..., -1.6]^2 -> lg R_AB table
                // TO BE CHECKED
                std::vector<std::vector<double>> table = {
                    {-0.1480, -0.2358, -0.3454, -0.4843, -0.6628, -0.8931, -1.1920, -1.5829,
                     -2.0996, -2.7905, -3.7203, -4.9696, -6.6362, -8.8387, -11.7229, -15.4698},
                    {-0.2497, -0.3467, -0.4664, -0.6159, -0.8028, -1.0372, -1.3340, -1.7157,
                     -2.2164, -2.8870, -3.7959, -5.0266, -6.6779, -8.8687, -11.7441, -15.4846},
                    {-0.3762, -0.4838, -0.6161, -0.7786, -0.9770, -1.2183, -1.5140, -1.8842,
                     -2.3638, -3.0076, -3.8895, -5.0966, -6.7291, -8.9055, -11.7703, -15.5031},
                    {-0.5392, -0.6563, -0.8008, -0.9784, -1.1931, -1.4475, -1.7463, -2.1042,
                     -2.5562, -3.1641, -4.0100, -5.1864, -6.7946, -8.9527, -11.8041, -15.5271},
                    {-0.7536, -0.8760, -1.0284, -1.2187, -1.4516, -1.7275, -2.0416, -2.3939,
                     -2.8138, -3.3736, -4.1706, -5.3054, -6.8813, -9.0153, -11.8489, -15.5591},
                    {-1.0400, -1.1637, -1.3149, -1.5071, -1.7492, -2.0456, -2.3905, -2.7635,
                     -3.1622, -3.6632, -4.3927, -5.4691, -7.0001, -9.1008, -11.9103, -15.6029},
                    {-1.4429, -1.5528, -1.6915, -1.8697, -2.0996, -2.3929, -2.7568, -3.1806,
                     -3.6166, -4.0735, -4.7127, -5.7034, -7.1680, -9.2206, -11.9958, -15.6639},
                    {-1.9984, -2.0919, -2.2095, -2.3609, -2.5584, -2.8168, -3.1538, -3.5870,
                     -4.1108, -4.6351, -5.1916, -6.0538, -7.4142, -9.3931, -12.1173, -15.7500},
                    {-2.7712, -2.8458, -2.9393, -3.0594, -3.2166, -3.4248, -3.7030, -4.0764,
                     -4.5818, -5.2366, -5.8918, -6.6043, -7.7919, -9.6493, -12.2935, -15.8731},
                    {-3.8375, -3.8944, -3.9653, -4.0564, -4.1758, -4.3350, -4.5502, -4.8440,
                     -5.2498, -5.8270, -6.6393, -7.4719, -8.4097, -10.0489, -12.5574, -16.0523},
                    {-5.2898, -5.3317, -5.3840, -5.4512, -5.5395, -5.6577, -5.8185, -6.0401,
                     -6.3474, -6.7806, -7.4218, -8.4047, -9.4662, -10.7245, -12.9735, -16.3217},
                    {-7.2398, -7.2700, -7.3078, -7.3565, -7.4208, -7.5072, -7.6253, -7.7883,
                     -8.0154, -8.3335, -8.7876, -9.4801, -10.6314, -11.9835, -13.6956, -16.7499},
                    {-9.8243, -9.8457, -9.8726, -9.9075, -9.9539, -10.0166, -10.1023, -10.2209,
                     -10.3860, -10.6174, -10.9431, -11.4119, -12.1425, -13.4464, -15.1608, -17.5078},
                    {-13.2126, -13.2276, -13.2466, -13.2715, -13.3047, -13.3498, -13.4117, -13.4972,
                     -13.6167, -13.7836, -14.0179, -14.3487, -14.8273, -15.5846, -17.0157, -19.1725},
                    {-17.6168, -17.6272, -17.6405, -17.6580, -17.6816, -17.7139, -17.7584, -17.8201,
                     -17.9060, -18.0259, -18.1939, -18.4302, -18.7643, -19.2494, -20.0245, -21.5537},
                    {-23.3023, -23.3095, -23.3188, -23.3311, -23.3477, -23.3707, -23.4025, -23.4467,
                     -23.5084, -23.5945, -23.7149, -23.8836, -24.1212, -24.4575, -24.9466, -25.7333}};

                double base_exp = table[index_1][index_2];
                // linear corrections
                double tau_p_corr = (index_1 > 14) ? 0.0 : (table[index_1 + 1][index_2] - base_exp) / 0.1 * (-10 * log10(tau_p) - 1 - index_1);
                double tau_n_corr = (index_2 > 14) ? 0.0 : (table[index_1][index_2 + 1] - base_exp) / 0.1 * (-10 * log10(tau_n) - 1 - index_2);

                double r_comp = pow(10.0, base_exp + tau_p_corr + tau_n_corr);

                return (pow(tau_n, 2.0) + pow(tau_p, 2.0) < 3 * 3) ? r_comp : r_comp * exp(-sqrt(pow(tau_n, 2.0) + pow(tau_p, 2.0)) / 3.0);
            };
        }
        return 0.0;
    };
}

std::function<double(double, const auxiliaries::Species &, double, double)> cooling::predefined::neutrinic::hadron_murca_emissivity(
    const std::map<auxiliaries::Species, std::function<double(double)>> &k_fermi_of_nbar,
    const std::map<auxiliaries::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
    double nbar_core_limit, double nbar_conversion, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
    const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp)
{
    return [=](double r, const auxiliaries::Species &lepton_flavour, double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::conversion;

        double nbar_val = nbar_of_r(r);
        double pf_l = k_fermi_of_nbar.at(lepton_flavour)(nbar_val),
               pf_n = k_fermi_of_nbar.at(neutron)(nbar_val),
               pf_p = k_fermi_of_nbar.at(proton)(nbar_val);
        double mst_n = m_stars_of_nbar.at(neutron)(nbar_val),
               mst_p = m_stars_of_nbar.at(proton)(nbar_val),
               mst_l = m_stars_of_nbar.at(lepton_flavour)(nbar_val);
        double T_loc = T / exp_phi(r);
        double alpha = 1.76 - 0.63 * pow(N_sat / (nbar_val * nbar_conversion), 2.0 / 3), beta = 0.68,
               v_fl = pf_l / mst_l;
        double dens_n = (8.05E21 / 1.68E72) * v_fl * pow(mst_n / M_N, 3) * (mst_p / M_N) * pf_p *
                        pow(T_loc, 8) * alpha * beta *
                        pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);
        double dens_p = (pf_l + 3 * pf_p - pf_n > 0) ? (8.05E21 / (8 * 1.68E72)) * (pow(pf_l + 3 * pf_p - pf_n, 2) / mst_l) * pow(mst_p / M_N, 3) * (mst_n / M_N) *
                                                           pow(T_loc, 8) * alpha * beta *
                                                           pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3)
                                                     : 0.0;

        // Superfluid factors
        double r_Mn_n = 1.0, r_Mn_p = 1.0, r_Mp_n = 1.0, r_Mp_p = 1.0;
        auto r_Mn_n_1S0 = [&](double v)
        {
            return pow(0.2414 + sqrt(pow(0.7586, 2.0) + pow(0.1318 * v, 2.0)), 7.0) * exp(5.339 - sqrt(pow(5.339, 2.0) + pow(2 * v, 2.0)));
        };
        auto r_Mn_p_1S0 = [&](double v)
        {
            double a = 0.1477 + sqrt(pow(0.8523, 2.0) + pow(0.1175 * v, 2.0)),
                   b = 0.1477 + sqrt(pow(0.8523, 2.0) + pow(0.1297 * v, 2.0));
            return (pow(a, 7.5) + pow(b, 5.5)) / 2.0 * exp(3.4370 - sqrt(pow(3.4370, 2.0) + pow(v, 2.0)));
        };
        auto r_Mp_p_1S0 = r_Mn_n_1S0;
        auto r_Mp_n_1S0 = r_Mn_p_1S0;
        auto r_Mp_n_3P2 = [&](double v)
        {
            double a = 0.1612 + sqrt(pow(0.8388, 2.0) + pow(0.1117 * v, 2.0)),
                   b = 0.1612 + sqrt(pow(0.8388, 2.0) + pow(0.1274 * v, 2.0));
            return (pow(a, 7.0) + pow(b, 5.0)) / 2.0 * exp(2.398 - sqrt(pow(2.398, 2.0) + pow(v, 2.0)));
        };
        auto r_Mn_n_3P2 = r_Mp_p_1S0; // Following Yakovlev similarity criteria

        // proton superfluidity?
        if (superfluid_p_1s0)
        {
            double tau = T_loc / superfluid_p_temp(pf_p);
            if (tau < 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                r_Mn_p = r_Mn_p_1S0(superfluid_gap_1s0(tau));
                r_Mp_p = r_Mp_p_1S0(superfluid_gap_1s0(tau));
            }
        }
        // neutron superfluidity?
        if (superfluid_n_3p2 || superfluid_n_1s0)
        {
            double tau = T_loc / superfluid_n_temp(pf_n);
            if (tau < 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                //  1S0/3P2 division at core entrance
                if (superfluid_n_3p2 && superfluid_n_1s0)
                {
                    r_Mn_n = (nbar_val > nbar_core_limit) ? r_Mn_n_3P2(superfluid_gap_3p2(tau)) : r_Mn_n_1S0(superfluid_gap_1s0(tau));
                    r_Mp_n = (nbar_val > nbar_core_limit) ? r_Mp_n_3P2(superfluid_gap_3p2(tau)) : r_Mp_n_1S0(superfluid_gap_1s0(tau));
                }
                // 1S0 only
                else if (superfluid_n_1s0)
                {
                    r_Mn_n = r_Mn_n_1S0(superfluid_gap_1s0(tau));
                    r_Mp_n = r_Mp_n_1S0(superfluid_gap_1s0(tau));
                }
                // 3P2 only
                else
                {
                    r_Mn_n = r_Mn_n_3P2(superfluid_gap_3p2(tau));
                    r_Mp_n = r_Mp_n_3P2(superfluid_gap_3p2(tau));
                }
            }
        }
        return dens_n * std::min(r_Mn_n, r_Mp_n) + dens_p * std::min(r_Mn_p, r_Mp_p);
    };
}

std::function<double(double, double, double)> cooling::predefined::neutrinic::hadron_bremsstrahlung_emissivity(
    const std::map<auxiliaries::Species, std::function<double(double)>> &k_fermi_of_nbar,
    const std::map<auxiliaries::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
    const std::function<double(double)> &ion_volume_frac, double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
    const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp)
{
    return [=](double r, double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::conversion;

        double nbar_val = nbar_of_r(r);
        double pf_n = k_fermi_of_nbar.at(neutron)(nbar_val),
               pf_p = k_fermi_of_nbar.at(proton)(nbar_val);
        double mst_n = m_stars_of_nbar.at(neutron)(nbar_val),
               mst_p = m_stars_of_nbar.at(proton)(nbar_val);
        double alpha_nn = 0.59, alpha_np = 1.06, alpha_pp = 0.11,
               beta_nn = 0.56, beta_np = 0.66, beta_pp = 0.7;
        double T_loc = T / exp_phi(r);
        int n_flavours = 3;
        double dens_nn = (7.5E19 / 1.68E72) * pow(mst_n / M_N, 4) * pf_n * n_flavours *
                         pow(T_loc, 8) * alpha_nn * beta_nn *
                         pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3),
               dens_pp = (7.5E19 / 1.68E72) * pow(mst_p / M_N, 4) * pf_p * n_flavours *
                         pow(T_loc, 8) * alpha_pp * beta_pp *
                         pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3),
               dens_np = (1.5E20 / 1.68E72) * pow(mst_p / M_N, 2) * pow(mst_n / M_N, 2) * pf_p * n_flavours *
                         pow(T_loc, 8) * alpha_np * beta_np *
                         pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);

        // Suppression due to ion excluded volume
        double eta_ion = 0.0; // if p1s0 is not enabled, we should account for this
        if (!superfluid_p_1s0)
            eta_ion = ion_volume_frac(nbar_val);

        // Superfluid factors
        double r_nn_n = 1.0, r_nn_p = 1.0,
               r_np_n = 1.0, r_np_p = 1.0,
               r_pp_n = 1.0, r_pp_p = 1.0;
        // r_nn_p, r_pp_n, are not affected by superfluidity, they are left as 1.0

        auto r_nn_n_1S0 = [&](double v)
        {
            double c = 0.1747 + sqrt(pow(0.8253, 2.0) + pow(0.07933 * v, 2.0)),
                   d = 0.7333 + sqrt(pow(0.2667, 2.0) + pow(0.1678 * v, 2.0));
            return 0.5 * (pow(c, 2.0) * exp(4.228 - sqrt(pow(4.228, 2.0) + pow(2 * v, 2.0))) +
                          pow(d, 7.5) * exp(7.762 - sqrt(pow(7.762, 2.0) + pow(3 * v, 2.0))));
        };
        auto r_pp_p_1S0 = r_nn_n_1S0;
        auto r_np_p_1S0 = [&](double v)
        {
            double a = 0.9982 + sqrt(pow(0.0018, 2.0) + pow(0.3815 * v, 2.0)),
                   b = 0.3949 + sqrt(pow(0.6051, 2.0) + pow(0.2666 * v, 2.0));
            return 1.0 / 2.732 * (a * exp(1.306 - sqrt(pow(1.306, 2.0) + pow(v, 2.0))) + 1.732 * pow(b, 7.0) * exp(3.303 - sqrt(pow(3.303, 2.0) + pow(4 * v, 2.0))));
        };
        auto r_np_n_1S0 = r_np_p_1S0;
        auto r_np_n_3P2 = r_np_n_1S0;
        auto r_nn_n_3P2 = r_nn_n_1S0;

        // proton superfluidity?
        if (superfluid_p_1s0)
        {
            double tau = T_loc / superfluid_p_temp(pf_p);
            if (tau < 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                r_pp_p = r_pp_p_1S0(superfluid_gap_1s0(tau));
                r_np_p = r_np_p_1S0(superfluid_gap_1s0(tau));
            }
        }

        // neutron superfluidity?
        if (superfluid_n_3p2 || superfluid_n_1s0)
        {
            double tau = T_loc / superfluid_n_temp(pf_n);
            if (tau < 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                // 1S0/3P2 division at core entrance
                if (superfluid_n_3p2 && superfluid_n_1s0)
                {
                    r_nn_n = (nbar_val > nbar_core_limit) ? r_nn_n_3P2(superfluid_gap_3p2(tau)) : r_nn_n_1S0(superfluid_gap_1s0(tau));
                    r_np_n = (nbar_val > nbar_core_limit) ? r_np_n_3P2(superfluid_gap_3p2(tau)) : r_np_n_1S0(superfluid_gap_1s0(tau));
                }
                // 1S0 only
                else if (superfluid_n_1s0)
                {
                    r_nn_n = r_nn_n_1S0(superfluid_gap_1s0(tau));
                    r_np_n = r_np_n_1S0(superfluid_gap_1s0(tau));
                }
                // 3P2 only
                else
                {
                    r_nn_n = r_nn_n_3P2(superfluid_gap_3p2(tau));
                    r_np_n = r_np_n_3P2(superfluid_gap_3p2(tau));
                }
            }
        }
        double r_nn = std::min(r_nn_n, r_nn_p),
               r_np = std::min(r_np_n, r_np_p),
               r_pp = std::min(r_pp_n, r_pp_p);

        return dens_nn * r_nn * (1 - eta_ion) + dens_np * r_np + dens_pp * r_pp;
    };
}

std::function<double(double, const auxiliaries::Species &, double, double)> cooling::predefined::neutrinic::hadron_pbf_emissivity(
    const std::map<auxiliaries::Species, std::function<double(double)>> &k_fermi_of_nbar,
    const std::map<auxiliaries::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
    double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
    const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp)
{
    return [=](double r, const auxiliaries::Species &hadron, double t, double T)
    {
        // the process is not allowed in normal matter
        if (!(superfluid_n_3p2 || superfluid_n_1s0 || superfluid_p_1s0))
            return 0.0;
        using namespace constants::scientific;
        using namespace constants::conversion;

        double nbar_val = nbar_of_r(r);
        double pf = k_fermi_of_nbar.at(hadron)(nbar_val);
        double mst = m_stars_of_nbar.at(hadron)(nbar_val);
        double a_s, a_t;
        if (hadron == neutron)
        {
            a_s = 1.0 + 1.588 * pow(pf / M_N, 2.0) * (1.0 + 0.262 * pow(mst / M_N, -2.0));
            a_t = 4.17;
        }
        else if (hadron == proton)
        {
            a_s = 0.0064 + 1.588 * pow(pf / M_N, 2.0) * (1.0 + 0.262 * pow(mst / M_N, -2.0));
            a_t = 0.0;
        }
        else
            throw std::runtime_error("Unexpected species: " + std::to_string(static_cast<int>(hadron.identify())) + "; Encountered in hadron_PBF_emissivity");
        double T_loc = T / exp_phi(r);
        int n_flavours = 3;
        double base_dens = (1.17E21 / 1.0E63) * (mst / M_N) * (pf / M_N) * n_flavours *
                           pow(T_loc, 7) *
                           pow(gev_over_k, 7) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);

        // Superfluid factors
        auto f_s = [&](double v)
        {
            return (0.602 * pow(v, 2.0) + 0.5942 * pow(v, 4.0) + 0.288 * pow(v, 6.0)) *
                   pow(0.5547 + sqrt(pow(0.4453, 2.0) + 0.01130 * pow(v, 2.0)), 0.5) *
                   exp(2.245 - sqrt(pow(2.245, 2.0) + pow(2 * v, 2.0)));
        };
        auto f_t = [&](double v)
        {
            return (1.204 * pow(v, 2.0) + 3.733 * pow(v, 4.0) + 0.3191 * pow(v, 6.0)) / (1 + 0.3155 * pow(v, 2.0)) *
                   pow(0.7591 + sqrt(pow(0.2409, 2.0) + 0.3145 * pow(v, 2.0)), 2.0) *
                   exp(0.4616 - sqrt(pow(0.4616, 2.0) + pow(2 * v, 2.0)));
        };

        // proton superfluidity?
        if (superfluid_p_1s0 && (hadron == proton))
        {
            double tau = T_loc / superfluid_p_temp(pf);
            if (tau < 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                return base_dens * a_s * f_s(superfluid_gap_1s0(tau));
            }
        }

        // neutron superfluidity?
        if ((superfluid_n_3p2 || superfluid_n_1s0) && (hadron == neutron))
        {
            double tau = T_loc / superfluid_n_temp(pf);
            if (tau < 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                // 1S0/3P2 division at core entrance
                if (superfluid_n_3p2 && superfluid_n_1s0)
                {
                    return base_dens * ((nbar_val > nbar_core_limit) ? a_t * f_t(superfluid_gap_3p2(tau)) : a_s * f_s(superfluid_gap_1s0(tau)));
                }
                // 1S0 only
                else if (superfluid_n_1s0)
                {
                    return base_dens * a_s * f_s(superfluid_gap_1s0(tau));
                }
                // 3P2 only
                else
                {
                    return base_dens * a_t * f_t(superfluid_gap_3p2(tau));
                }
            }
        }
        return 0.0;
    };
}

std::function<double(double, double, double)> cooling::predefined::neutrinic::quark_durca_emissivity(
    const std::map<auxiliaries::Species, std::function<double(double)>> &k_fermi_of_nbar,
    const std::map<auxiliaries::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
    const std::function<double(double)> &exp_phi, bool superconduct_u, bool superconduct_d, bool superconduct_s,
    const std::function<double(double)> &superconduct_u_temp, const std::function<double(double)> &superconduct_d_temp, const std::function<double(double)> &superconduct_s_temp)
{
    return [=](double r, double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::conversion;

        double nbar_val = nbar_of_r(r);
        double pf_u = k_fermi_of_nbar.at(uquark)(nbar_val),
               pf_d = k_fermi_of_nbar.at(dquark)(nbar_val),
               pf_s = k_fermi_of_nbar.at(squark)(nbar_val);
        double pf_l = 0.0,
               pf_l_mult_onemincos = 0.0;
        double alpha_c = 1.0, m_s = 0.01;
        double T_loc = T / exp_phi(r);

        double dens_ud = 3.74E-10 * alpha_c * pf_u * pf_d * pow(T_loc, 6);
        for (auto it = k_fermi_of_nbar.begin(); it != k_fermi_of_nbar.end(); ++it)
        {
            if (it->first.classify() == auxiliaries::Species::ParticleClassification::kLepton)
                pf_l += it->second(nbar_val);
        }
        dens_ud *= pf_l;

        double x = pf_s / m_s, eta = sqrt(1 + x * x);
        double mu_s = (x > 0.001) ? (eta / x + 8 * alpha_c / (3 * Pi) * (1 - 3 / (eta * x) * log(x + eta))) * pf_s : m_s;
        double dens_us = 1.208E-11 * mu_s * pf_u * pow(T_loc, 6);
        for (auto it = k_fermi_of_nbar.begin(); it != k_fermi_of_nbar.end(); ++it)
        {
            if (it->first.classify() == auxiliaries::Species::ParticleClassification::kLepton)
            {
                double pf_l = it->second(nbar_val);
                if (fabs(pf_s - pf_u) > pf_l)
                    continue;
                if (pf_s * pf_u * pf_l == 0.0)
                    continue;
                // Hopefully this is an enough contraint

                double theta_14 = acos((pf_s * pf_s + pf_l * pf_l - pf_u * pf_u) / (2.0 * pf_s * pf_l)),
                       theta_13 = acos((pf_s * pf_s + pf_u * pf_u - pf_l * pf_l) / (2.0 * pf_s * pf_u));
                pf_l_mult_onemincos += pf_l * (1.0 - cos(theta_14 + theta_13));
            }
        }
        dens_us *= pf_l_mult_onemincos;

        double dens = dens_ud + dens_us;

        // superconductivity?
        if (superconduct_u)
        {
            double tau = T_loc / superconduct_u_temp(pf_u);
            if (tau < 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                dens *= exp(-1.76 / tau * sqrt(1.0 - tau));
            }
        }
        return dens;
    };
}

std::function<double(double, double, double)> cooling::predefined::neutrinic::quark_murca_emissivity(
    const std::map<auxiliaries::Species, std::function<double(double)>> &k_fermi_of_nbar,
    const std::map<auxiliaries::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
    const std::function<double(double)> &exp_phi, bool superconduct_u, bool superconduct_d, bool superconduct_s,
    const std::function<double(double)> &superconduct_u_temp, const std::function<double(double)> &superconduct_d_temp, const std::function<double(double)> &superconduct_s_temp)
{
    return [=](double r, double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::conversion;

        double nbar_val = nbar_of_r(r);
        double pf_u = k_fermi_of_nbar.at(uquark)(nbar_val),
               pf_d = k_fermi_of_nbar.at(dquark)(nbar_val);
        double alpha_c = 1.0;
        double T_loc = T / exp_phi(r);

        // Estimates
        double dens = alpha_c * alpha_c * 1.29E-10 * (pf_u + pf_d) * pow(T_loc, 8);

        // superconductivity?
        if (superconduct_u || superconduct_d)
        {
            double tau_u_inv = superconduct_u_temp(pf_u)/T_loc,
                     tau_d_inv = superconduct_d_temp(pf_d)/T_loc;
            if (tau_u_inv > 1.0 && tau_d_inv < 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                auto tau = 1 / tau_u_inv;
                dens *= exp(- 2.0 * 1.76 / tau * sqrt(1.0 - tau));
            }
            else if (tau_u_inv < 1.0 && tau_d_inv > 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                auto tau = 1 / tau_d_inv;
                dens *= exp(- 2.0 * 1.76 / tau * sqrt(1.0 - tau));
            }
            else if (tau_u_inv > 1.0 && tau_d_inv > 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                auto tau = 1 / tau_u_inv;
                dens *= exp(- 1.76 / tau * sqrt(1.0 - tau));
                tau = 1 / tau_d_inv;
                dens *= exp(- 1.76 / tau * sqrt(1.0 - tau));
            }
        }
        return dens;
    };
}

std::function<double(double, double, double)> cooling::predefined::neutrinic::quark_bremsstrahlung_emissivity(
    const std::map<auxiliaries::Species, std::function<double(double)>> &k_fermi_of_nbar,
    const std::map<auxiliaries::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
    const std::function<double(double)> &exp_phi, bool superconduct_u, bool superconduct_d, bool superconduct_s,
    const std::function<double(double)> &superconduct_u_temp, const std::function<double(double)> &superconduct_d_temp, const std::function<double(double)> &superconduct_s_temp)
{
    return [=](double r, double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::conversion;

        double nbar_val = nbar_of_r(r);
        double pf_u = k_fermi_of_nbar.at(uquark)(nbar_val),
               pf_d = k_fermi_of_nbar.at(dquark)(nbar_val);
        double alpha_c = 1.0;
        double T_loc = T / exp_phi(r);

        // Estimates
        double dens_u = 1.36E-10 * pf_u * pow(T_loc, 8),
                dens_d = 1.36E-10 * pf_d * pow(T_loc, 8);

        // superconductivity?
        if (superconduct_u)
        {
            auto tau = T_loc / superconduct_u_temp(pf_u);
            if (tau < 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                dens_u *= exp(- 2.0 * 1.76 / tau * sqrt(1.0 - tau));
            }
        }
        if (superconduct_d)
        {
            auto tau = T_loc / superconduct_d_temp(pf_d);
            if (tau < 1.0)
            {
                using namespace cooling::predefined::auxiliary;
                dens_d *= exp(- 2.0 * 1.76 / tau * sqrt(1.0 - tau));
            }
        }
        return dens_u + dens_d;
    };
}

std::function<double(double, double, double)> cooling::predefined::neutrinic::electron_bremsstrahlung_emissivity(
                const std::map<auxiliaries::Species, std::function<double(double)>> &k_fermi_of_nbar,
                const std::map<auxiliaries::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
                const std::function<double(double)> &exp_phi)
{
    return [=](double r, double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::conversion;

        double nbar_val = nbar_of_r(r);
        double pf_e = k_fermi_of_nbar.at(electron)(nbar_val);
        double T_loc = T / exp_phi(r);

        // Ignore superfluidity, since the reaction is not relevant in hadronic phase
        double dens = (2.427E16 / 1.0E72) * pf_e * pow(T_loc, 8) *
                        pow(gev_over_k, 8) * erg_over_gev / gev_s * (km_gev * 1.0E-18) / pow(km_gev * 1.0E-5, 3);

        return dens;
    };
}