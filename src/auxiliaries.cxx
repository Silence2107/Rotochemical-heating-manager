
#include "../include/auxiliaries.h"
#include "../include/constants.h"

#include <cmath>
#include <string>
#include <algorithm>
#include <vector>
#include <functional>
#include <fstream>
#include <sstream>
#include <map>

std::vector<std::vector<double>> auxiliaries::io::read_tabulated_file(const std::string &path, std::pair<size_t, size_t> columns, std::pair<size_t, size_t> rows, double empty_value)
{
    std::ifstream fstr(path);
    if (!fstr.is_open())
        throw std::runtime_error("Cannot open file " + path + ". Encountered in auxiliaries::io::read_tabulated_file");
    std::vector<std::vector<double>> table;
    std::vector<std::string> lines;
    std::string line;
    while (std::getline(fstr, line))
        lines.push_back(line);
    if (rows.second == 0)
        rows.second = lines.size();
    if (columns.second == 0)
    {
        std::string cleared_line = auxiliaries::io::retrieve_cleared_line(lines[rows.first]);
        std::stringstream ss(cleared_line);
        std::string str;
        while (ss >> str)
            ++columns.second;
    }
    if (rows.second <= rows.first || columns.second <= columns.first)
        throw std::runtime_error("Invalid rows or columns count extracted from input file. Encountered in auxiliaries::io::read_tabulated_file");
    std::vector<std::vector<std::string>> str_data(rows.second - rows.first, std::vector<std::string>(columns.second - columns.first));
    for (size_t i = rows.first; i < rows.second; ++i)
    {
        std::string cleared_line = auxiliaries::io::retrieve_cleared_line(lines[i]);
        std::stringstream ss(cleared_line);
        for (size_t j = columns.first; j < columns.second; ++j)
        {
            std::string str;
            if (!(ss >> str))
                str = std::to_string(empty_value);
            // skip columns before the first one
            if (j < columns.first)
                continue;
            str_data[i - rows.first][j - columns.first] = str;
        }
    }
    for (size_t i = columns.first; i < columns.second; ++i)
    {
        std::vector<double> column(rows.second - rows.first);
        for (size_t j = rows.first; j < rows.second; ++j)
        {
            try
            {
                column[j - rows.first] = std::stod(str_data[j - rows.first][i - columns.first]);
            }
            catch (const std::invalid_argument &e)
            {
                throw std::runtime_error("Intractable argument: " + std::string(e.what()) + ". Encountered in auxiliaries::io::read_tabulated_file");
            }
        }
        table.push_back(column);
    }
    return table;
}

std::string auxiliaries::io::retrieve_cleared_line(const std::string &line)
{
    if (line.empty())
        return line;
    std::string cleared_line{line};
    for (int i = 0; i < cleared_line.length() - 1; ++i)
    {
        if (cleared_line[i] <= 32 && cleared_line[i + 1] <= 32)
        {                             // if two characters in a row are auxiliary symbols or whitespaces
            cleared_line.erase(i, 1); // erase one
            cleared_line[i] = ' ';    // replace with whitespace
            --i;                      // and get one position back
        }
    } // strips additional whitespaces
    if (cleared_line[0] == ' ')
        cleared_line.erase(0, 1); // erase first whitespace if present
    if (cleared_line.back() == ' ')
        cleared_line.pop_back(); // erase last whitespace   if present
    return cleared_line;
}

double auxiliaries::math::interpolate(const std::vector<double> &input, const std::vector<double> &output,
                                auxiliaries::math::InterpolationMode mode, double x, bool extrapolate, bool enable_checks)
{
    if (enable_checks)
    {
        if (input.size() != output.size())
            throw std::runtime_error("Input and output arrays have different sizes. Encountered in auxiliaries::math::interpolate");
    }
    // determine if input is increasing or decreasing
    bool decreasingly_sorted = (input[0] > input[1]) ? true : false;
    if (!extrapolate && enable_checks)
        if ((!decreasingly_sorted && (x < input.front() || x > input.back())) || (decreasingly_sorted && (x > input.front() || x < input.back())))
            throw std::runtime_error("Searched value is out of range. Encountered in auxiliaries::math::interpolate");

    // find index of x in input sorted array in reasonable time
    size_t low_pos;
    if (!decreasingly_sorted)
        low_pos = std::upper_bound(input.begin(), input.end(), x) -
                  input.begin();
    else
        low_pos = std::upper_bound(input.begin(), input.end(), x, std::greater<double>()) -
                  input.begin();
    // extrapolation: in case x is out of range, we extrapolate from the nearest polinomial
    if(extrapolate)
    {
        if (low_pos == input.size())
            --low_pos;
        if (low_pos == 0)
            ++low_pos;
    }
    // Conventional definition of low_pos : x is within [input[low_pos], input[low_pos+1]), if not extrapolating
    --low_pos;
    // if x is equal to input[i], return output[i]
    if (x == input[low_pos])
        return output[low_pos];

    // interpolate
    switch (mode)
    {
    case auxiliaries::math::InterpolationMode::kLinear:
    {
        if (enable_checks)
            if (input.size() < 2)
                throw std::runtime_error("Cannot perform linear interpolation with less than 2 points. Encountered in auxiliaries::math::interpolate");
        return output[low_pos] + (output[low_pos + 1] - output[low_pos]) * (x - input[low_pos]) / (input[low_pos + 1] - input[low_pos]);
    }
    break;
    case auxiliaries::math::InterpolationMode::kCubic:
    {
        if (enable_checks)
            if (input.size() < 5)
                throw std::runtime_error("Cannot perform cubic interpolation with less than 5 points. Encountered in auxiliaries::math::interpolate");
        auto tridiagonal_solve = [](const std::vector<double> &subdiag, const std::vector<double> &diag, const std::vector<double> &superdiag, const std::vector<double> &rhs)
        {
            size_t n = diag.size();
            std::vector<double> v(n);     // The solution vector
            std::vector<double> c(n - 1); // new superdiagonal
            std::vector<double> g(n);     // new right hand side
            c[0] = superdiag[0] / diag[0];
            g[0] = rhs[0] / diag[0];
            for (size_t i = 1; i < n - 1; ++i)
            {
                c[i] = superdiag[i] / (diag[i] - subdiag[i - 1] * c[i - 1]);
                g[i] = (rhs[i] - subdiag[i - 1] * g[i - 1]) / (diag[i] - subdiag[i - 1] * c[i - 1]);
            }
            g[n - 1] = (rhs[n - 1] - subdiag[n - 2] * g[n - 2]) / (diag[n - 1] - subdiag[n - 2] * c[n - 2]);
            v[n - 1] = g[n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                v[i] = g[i] - c[i] * v[i + 1];
            }
            return v;
        };
        // Solve for quadratic coefficients
        std::vector<double> subdiag(input.size() - 4);
        std::vector<double> diag(input.size() - 3);
        std::vector<double> superdiag(input.size() - 4);
        std::vector<double> rhs(input.size() - 3);
        for (size_t i = 0; i < input.size() - 4; ++i)
        {
            subdiag[i] = (input[i + 2] - input[i + 1]);
            superdiag[i] = (input[i + 2] - input[i + 1]);
        }
        for (size_t i = 0; i < input.size() - 3; ++i)
        {
            diag[i] = 2 * (input[i + 2] - input[i]);
            rhs[i] = 3 * ((output[i + 2] - output[i + 1]) / (input[i + 2] - input[i + 1]) - (output[i + 1] - output[i]) / (input[i + 1] - input[i]));
        }
        std::vector<double> quadratic_coeffs = tridiagonal_solve(subdiag, diag, superdiag, rhs);
        quadratic_coeffs.insert(quadratic_coeffs.begin(), 0);
        quadratic_coeffs.push_back(0);
        // Solve for cubic and linear coefficients
        std::vector<double> cubic_coeffs(input.size() - 1), linear_coeffs(input.size() - 1);
        for (size_t i = 0; i < input.size() - 2; ++i)
        {
            cubic_coeffs[i] = (quadratic_coeffs[i + 1] - quadratic_coeffs[i]) / (3 * (input[i + 1] - input[i]));
            linear_coeffs[i] = (output[i + 1] - output[i]) / (input[i + 1] - input[i]) - (input[i + 1] - input[i]) * (2 * quadratic_coeffs[i] + quadratic_coeffs[i + 1]) / 3;
        }
        linear_coeffs[input.size() - 2] = linear_coeffs[input.size() - 3] + 2 * (input[input.size() - 2] - input[input.size() - 3]) * quadratic_coeffs[input.size() - 3] + 3 * (input[input.size() - 2] - input[input.size() - 3]) * (input[input.size() - 2] - input[input.size() - 3]) * cubic_coeffs[input.size() - 3];
        cubic_coeffs[input.size() - 2] = (output[input.size() - 1] - output[input.size() - 2] - linear_coeffs[input.size() - 2] * (input[input.size() - 1] - input[input.size() - 2])) / ((input[input.size() - 1] - input[input.size() - 2]) * (input[input.size() - 1] - input[input.size() - 2]) * (input[input.size() - 1] - input[input.size() - 2]));
        return cubic_coeffs[low_pos] * (x - input[low_pos]) * (x - input[low_pos]) * (x - input[low_pos]) + quadratic_coeffs[low_pos] * (x - input[low_pos]) * (x - input[low_pos]) + linear_coeffs[low_pos] * (x - input[low_pos]) + output[low_pos];
    }
    break;
    default:
        throw std::runtime_error("Unknown interpolation mode. Encountered in auxiliaries::math::interpolate");
    }
}

double auxiliaries::math::interpolate_cached(std::function<double(double)> &cache, const std::vector<double> &input, const std::vector<double> &output, auxiliaries::math::InterpolationMode mode, double x, bool extrapolate, bool enable_checks)
{
    // determine if input is increasing or decreasing
    bool decreasingly_sorted = (input[0] > input[1]) ? true : false;
    if (!extrapolate && enable_checks)
        if ((!decreasingly_sorted && (x < input.front() || x > input.back())) || (decreasingly_sorted && (x > input.front() || x < input.back())))
            throw std::runtime_error("Searched value is out of range. Encountered in auxiliaries::math::interpolate_cached");
    if (!cache) // if a callable is not stored, cache one
    {
        if (enable_checks)
        {
            if (input.size() != output.size())
                throw std::runtime_error("Input and output arrays have different sizes. Encountered in auxiliaries::math::interpolate_cached");
        }
        
        // interpolate
        switch (mode)
        {
        case auxiliaries::math::InterpolationMode::kLinear:
        {
            if (enable_checks)
                if (input.size() < 2)
                    throw std::runtime_error("Cannot perform linear interpolation with less than 2 points. Encountered in auxiliaries::math::interpolate_cached");
            cache = [=](double x)
            {
                // find index of x in input sorted array in reasonable time
                size_t low_pos;
                if (!decreasingly_sorted)
                    low_pos = std::upper_bound(input.begin(), input.end(), x) -
                              input.begin();
                else
                    low_pos = std::upper_bound(input.begin(), input.end(), x, std::greater<double>()) -
                              input.begin();
                // extrapolation: in case x is out of range, we extrapolate from the nearest polinomial
                if (extrapolate)
                {
                    if (low_pos == input.size())
                        --low_pos;
                    if (low_pos == 0)
                        ++low_pos;
                }
                // Conventional definition of low_pos : x is within [input[low_pos], input[low_pos+1]), if not extrapolating
                --low_pos;
                // if x is equal to input[i], return output[i]
                if (x == input[low_pos])
                    return output[low_pos];

                return output[low_pos] + (output[low_pos + 1] - output[low_pos]) * (x - input[low_pos]) / (input[low_pos + 1] - input[low_pos]);
            };
        }
        break;
        case auxiliaries::math::InterpolationMode::kCubic:
        {
            if (enable_checks)
                if (input.size() < 5)
                    throw std::runtime_error("Cannot perform cubic interpolation with less than 5 points. Encountered in auxiliaries::math::interpolate_cached");
            auto tridiagonal_solve = [](const std::vector<double> &subdiag, const std::vector<double> &diag, const std::vector<double> &superdiag, const std::vector<double> &rhs)
            {
                size_t n = diag.size();
                std::vector<double> v(n);     // The solution vector
                std::vector<double> c(n - 1); // new superdiagonal
                std::vector<double> g(n);     // new right hand side
                c[0] = superdiag[0] / diag[0];
                g[0] = rhs[0] / diag[0];
                for (size_t i = 1; i < n - 1; ++i)
                {
                    c[i] = superdiag[i] / (diag[i] - subdiag[i - 1] * c[i - 1]);
                    g[i] = (rhs[i] - subdiag[i - 1] * g[i - 1]) / (diag[i] - subdiag[i - 1] * c[i - 1]);
                }
                g[n - 1] = (rhs[n - 1] - subdiag[n - 2] * g[n - 2]) / (diag[n - 1] - subdiag[n - 2] * c[n - 2]);
                v[n - 1] = g[n - 1];
                for (int i = n - 2; i >= 0; --i)
                {
                    v[i] = g[i] - c[i] * v[i + 1];
                }
                return v;
            };
            // Solve for quadratic coefficients
            std::vector<double> subdiag(input.size() - 4);
            std::vector<double> diag(input.size() - 3);
            std::vector<double> superdiag(input.size() - 4);
            std::vector<double> rhs(input.size() - 3);
            for (size_t i = 0; i < input.size() - 4; ++i)
            {
                subdiag[i] = (input[i + 2] - input[i + 1]);
                superdiag[i] = (input[i + 2] - input[i + 1]);
            }
            for (size_t i = 0; i < input.size() - 3; ++i)
            {
                diag[i] = 2 * (input[i + 2] - input[i]);
                rhs[i] = 3 * ((output[i + 2] - output[i + 1]) / (input[i + 2] - input[i + 1]) - (output[i + 1] - output[i]) / (input[i + 1] - input[i]));
            }
            std::vector<double> quadratic_coeffs = tridiagonal_solve(subdiag, diag, superdiag, rhs);
            quadratic_coeffs.insert(quadratic_coeffs.begin(), 0);
            quadratic_coeffs.push_back(0);
            // Solve for cubic and linear coefficients
            std::vector<double> cubic_coeffs(input.size() - 1), linear_coeffs(input.size() - 1);
            for (size_t i = 0; i < input.size() - 2; ++i)
            {
                cubic_coeffs[i] = (quadratic_coeffs[i + 1] - quadratic_coeffs[i]) / (3 * (input[i + 1] - input[i]));
                linear_coeffs[i] = (output[i + 1] - output[i]) / (input[i + 1] - input[i]) - (input[i + 1] - input[i]) * (2 * quadratic_coeffs[i] + quadratic_coeffs[i + 1]) / 3;
            }
            linear_coeffs[input.size() - 2] = linear_coeffs[input.size() - 3] + 2 * (input[input.size() - 2] - input[input.size() - 3]) * quadratic_coeffs[input.size() - 3] +
                                              3 * (input[input.size() - 2] - input[input.size() - 3]) * (input[input.size() - 2] - input[input.size() - 3]) * cubic_coeffs[input.size() - 3];
            cubic_coeffs[input.size() - 2] = (output[input.size() - 1] - output[input.size() - 2] -
                                              linear_coeffs[input.size() - 2] * (input[input.size() - 1] - input[input.size() - 2])) /
                                             ((input[input.size() - 1] - input[input.size() - 2]) * (input[input.size() - 1] - input[input.size() - 2]) * (input[input.size() - 1] - input[input.size() - 2]));
            cache = [=](double x)
            {
                // find index of x in input sorted array in reasonable time
                size_t low_pos;
                if (!decreasingly_sorted)
                    low_pos = std::upper_bound(input.begin(), input.end(), x) -
                              input.begin();
                else
                    low_pos = std::upper_bound(input.begin(), input.end(), x, std::greater<double>()) -
                              input.begin();
                // extrapolation: in case x is out of range, we extrapolate from the nearest polinomial
                if (extrapolate)
                {
                    if (low_pos == input.size())
                        --low_pos;
                    if (low_pos == 0)
                        ++low_pos;
                }
                // Conventional definition of low_pos : x is within [input[low_pos], input[low_pos+1]), if not extrapolating
                --low_pos;
                // if x is equal to input[i], return output[i]
                if (x == input[low_pos])
                    return output[low_pos];

                return cubic_coeffs[low_pos] * (x - input[low_pos]) * (x - input[low_pos]) * (x - input[low_pos]) + quadratic_coeffs[low_pos] * (x - input[low_pos]) * (x - input[low_pos]) +
                       linear_coeffs[low_pos] * (x - input[low_pos]) + output[low_pos];
            };
        }
        break;
        default:
            throw std::runtime_error("Unknown interpolation mode. Encountered in auxiliaries::math::interpolate_cached");
        }
    }
    return cache(x);
}

auxiliaries::phys::Species::Species(auxiliaries::phys::Species::ParticleType type, auxiliaries::phys::Species::ParticleClassification classification)
    : m_type(type), m_classification(classification)
{
}

bool auxiliaries::phys::Species::operator==(const auxiliaries::phys::Species &other) const
{
    return m_type == other.m_type;
}

bool auxiliaries::phys::Species::operator<(const auxiliaries::phys::Species &other) const
{
    return m_type < other.m_type;
}

std::function<double(double, double, double)> auxiliaries::phys::fermi_specific_heat_density(
    const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
    const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
    double nbar_core_limit, const std::function<double(double)> &exp_phi, bool superfluid_n_1s0, bool superfluid_p_1s0, bool superfluid_n_3p2,
    const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp, const std::function<double(double)> &superconduct_q_gap)
{
    return [=](double r, double t, double T)
    {
        using namespace constants::scientific;
        using namespace constants::species;

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
                    using namespace auxiliaries::phys;
                    diff *= r_A(superfluid_gap_1s0(tau));
                }
            }
            // neutron superfluidity?
            else if (key == neutron && (superfluid_n_3p2 || superfluid_n_1s0))
            {
                double tau = T_loc / superfluid_p_temp(k_fermi);
                if (tau < 1.0)
                {
                    using namespace auxiliaries::phys;
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
            // quark superconductivity?
            else if (key.classify() == auxiliaries::phys::Species::ParticleClassification::kQuark)
            {
                // following Blaschke except the gap is provided externally
                auto exp_factor = superconduct_q_gap(nbar_val) / T_loc;
                // estimate critical temperature as 0.15 GeV (consider later if we need Tc(nbar))
                if (exp_factor > 1.0)
                    diff *= 3.1 * pow(critical_temperature(0.0, CriticalTemperatureModel::kHadronToQGP) / T_loc, 2.5) * exp(-exp_factor);
            }
            cv_dens += diff;
        }
        return cv_dens;
    };
}

double auxiliaries::phys::te_tb_relation(double Tb, double R, double M, double eta)
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

double auxiliaries::phys::critical_temperature_smeared_guassian(double k_fermi, double temp_ampl, double k_offs, double k_width, double quad_skew)
{
    return temp_ampl * exp(-pow((k_fermi - k_offs) / k_width, 2) - quad_skew * pow((k_fermi - k_offs) / k_width, 4));
}

double auxiliaries::phys::superfluid_gap_1s0(double tau)
{
    return std::sqrt(1 - tau) * (1.456 - 0.157 / std::sqrt(tau) + 1.764 / tau);
}

double auxiliaries::phys::superfluid_gap_3p2(double tau)
{
    return std::sqrt(1 - tau) * (0.7893 + 1.188 / tau);
}

double auxiliaries::phys::critical_temperature(double k_fermi, auxiliaries::phys::CriticalTemperatureModel model)
{
    using namespace constants::conversion;
    using namespace auxiliaries::phys;
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
    case CriticalTemperatureModel::kHadronToQGP:
        return 0.15;
    default:
        throw std::runtime_error("Unknown critical temperature model");
    }
}