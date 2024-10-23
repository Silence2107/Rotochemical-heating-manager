
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
        RHM_THROW(std::runtime_error, "Cannot open file " + path + ". ");
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
        RHM_THROW(std::runtime_error, "Invalid rows or columns count extracted from input file at " + path + ".");
    std::vector<std::vector<std::string>> str_data(rows.second - rows.first, std::vector<std::string>(columns.second - columns.first));
    for (size_t i = rows.first; i < rows.second; ++i)
    {
        std::string cleared_line = auxiliaries::io::retrieve_cleared_line(lines[i]);
        std::stringstream ss(cleared_line);
        for (size_t j = 0; j < columns.second; ++j)
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
                RHM_THROW(std::runtime_error, "Intractable argument: " + std::string(e.what()) + ".");
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
            RHM_THROW(std::runtime_error, "Input and output arrays have different sizes.");
    }
    // determine if input is increasing or decreasing
    bool decreasingly_sorted = (input[0] > input[1]) ? true : false;
    if (!extrapolate && enable_checks)
        if ((!decreasingly_sorted && (x < input.front() || x > input.back())) || (decreasingly_sorted && (x > input.front() || x < input.back())))
            RHM_THROW(std::runtime_error, "Searched value is out of range.");

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

    // interpolate
    switch (mode)
    {
    case auxiliaries::math::InterpolationMode::kLinear:
    {
        if (enable_checks)
            if (input.size() < 2)
                RHM_THROW(std::runtime_error, "Cannot perform linear interpolation with less than 2 points.");
        return output[low_pos] + (output[low_pos + 1] - output[low_pos]) * (x - input[low_pos]) / (input[low_pos + 1] - input[low_pos]);
    }
    break;
    case auxiliaries::math::InterpolationMode::kCubic:
    {
        if (enable_checks)
            if (input.size() < 5)
                RHM_THROW(std::runtime_error, "Cannot perform cubic interpolation with less than 5 points.");
        auto tridiagonal_solve = [](const std::vector<double> &subdiag, const std::vector<double> &diag, const std::vector<double> &superdiag, const std::vector<double> &rhs)
        {
            auxiliaries::math::MatrixD A(diag.size(), diag.size(), 0.0);
            for (size_t i = 0; i < diag.size(); ++i)
            {
                A.at(i, i) = diag[i];
                if (i < diag.size() - 1)
                    A.at(i, i + 1) = superdiag[i];
                if (i > 0)
                    A.at(i, i - 1) = subdiag[i - 1];
            }
            return A.tridiagonal_solve(rhs);
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
        RHM_THROW(std::runtime_error, "Unknown interpolation mode.");
    }
}

double auxiliaries::math::interpolate_cached(std::function<double(double)> &cache, const std::vector<double> &input, const std::vector<double> &output, auxiliaries::math::InterpolationMode mode, double x, bool extrapolate, bool enable_checks)
{
    // determine if input is increasing or decreasing
    bool decreasingly_sorted = (input[0] > input[1]) ? true : false;
    if (!extrapolate && enable_checks)
        if ((!decreasingly_sorted && (x < input.front() || x > input.back())) || (decreasingly_sorted && (x > input.front() || x < input.back())))
            RHM_THROW(std::runtime_error, "Searched value is out of range.");
    if (!cache) // if a callable is not stored, cache one
    {
        if (enable_checks)
        {
            if (input.size() != output.size())
                RHM_THROW(std::runtime_error, "Input and output arrays have different sizes.");
        }

        // interpolate
        switch (mode)
        {
        case auxiliaries::math::InterpolationMode::kLinear:
        {
            if (enable_checks)
                if (input.size() < 2)
                    RHM_THROW(std::runtime_error, "Cannot perform linear interpolation with less than 2 points.");
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
        // This interpolation mode became possible only thanks to https://signalsmith-audio.co.uk/writing/2021/monotonic-smooth-interpolation/
        case auxiliaries::math::InterpolationMode::kCubic:
        {
            if (enable_checks)
                if (input.size() < 2)
                    RHM_THROW(std::runtime_error, "Cannot perform cubic interpolation with less than 2 points.");
            std::vector<double> linear_coeffs(input.size()),
                quadratic_coeffs(input.size() - 1),
                cubic_coeffs(input.size() - 1),
                slopes(input.size() - 1);
            for (size_t i = 0; i < input.size() - 1; ++i)
            {
                slopes[i] = (output[i + 1] - output[i]) / (input[i + 1] - input[i]);
            }
            // deal with linear coefficients, which correspond to the slopes at the nodes
            for (size_t i = 0; i < input.size(); ++i)
            {
                // edge cases shall have natural slopes
                if (i == 0)
                {
                    linear_coeffs[i] = slopes.front();
                    continue;
                }
                else if (i == input.size() - 1)
                {
                    linear_coeffs[i] = slopes.back();
                    continue;
                }
                // we enforce monotonicity by imposing extremum at discrete slope sign change
                else if (slopes[i - 1] * slopes[i] <= 0)
                {
                    linear_coeffs[i] = 0;
                    continue;
                }
                // averaging slopes with weights
                else
                    linear_coeffs[i] = (slopes[i - 1] * (input[i + 1] - input[i]) + slopes[i] * (input[i] - input[i - 1])) / (input[i + 1] - input[i - 1]);
                
                // In the latter case, make sure linear coefficients are not too large to avoid overshooting.
                // Note that the signs of the slopes and linear_coeffs are at this point the same.
                if (linear_coeffs[i] / slopes[i] > 3 || linear_coeffs[i] / slopes[i - 1] > 3)
                    linear_coeffs[i] = 3 * (linear_coeffs[i] > 0) * std::min(std::abs(slopes[i]), std::abs(slopes[i - 1]));
            }
            // deal with the rest
            for (size_t i = 0; i < input.size() - 1; ++i)
            {
                cubic_coeffs[i] = (linear_coeffs[i] + linear_coeffs[i + 1] - 2 * slopes[i]) / pow(input[i + 1] - input[i], 2);
                quadratic_coeffs[i] = (3 * slopes[i] - linear_coeffs[i + 1] - 2 * linear_coeffs[i]) / (input[i + 1] - input[i]);
            }
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

                double s = x - input[low_pos];

                return cubic_coeffs[low_pos] * s * s * s + quadratic_coeffs[low_pos] * s * s + linear_coeffs[low_pos] * s + output[low_pos];
            };
        }
        break;
        default:
            RHM_THROW(std::runtime_error, "Unknown interpolation mode.");
        }
    }
    return cache(x);
}

auxiliaries::phys::Species::Species(auxiliaries::phys::Species::ParticleType type, auxiliaries::phys::Species::ParticleClassification classification, const std::string &name, double mass, double qcharge, double bcharge)
    : m_type(type), m_classification(classification), m_name(name), m_mass(mass), m_qcharge(qcharge), m_bcharge(bcharge)
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
    double nbar_sf_shift, const std::function<double(double)> &exp_phi, const std::function<double(double)> &superfluid_p_temp, 
    const std::function<double(double)> &superfluid_n_temp, const std::function<double(double)> &superconduct_q_gap)
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
            if (key == proton)
            {
                double T_c = superfluid_p_temp(nbar_val);
                if (T_loc < T_c)
                {
                    double tau = T_loc / T_c;
                    using namespace auxiliaries::phys;
                    diff *= r_A(superfluid_gap_1s0(tau));
                }
            }
            // neutron superfluidity?
            else if (key == neutron)
            {
                double T_c = superfluid_n_temp(nbar_val);
                if (T_loc < T_c)
                {
                    double tau = T_loc / T_c;
                    using namespace auxiliaries::phys;
                    diff *= (nbar_val > nbar_sf_shift ? r_B(superfluid_gap_3p2(tau)) : r_A(superfluid_gap_1s0(tau)));
                }
            }
            // quark superconductivity?
            else if (key.classify() == auxiliaries::phys::Species::ParticleClassification::kQuark)
            {
                // following Blaschke except the gap is provided externally
                auto exp_factor = superconduct_q_gap(nbar_val) / T_loc;
                // estimate Tc = 0.4 Delta
                if (exp_factor > 1.0)
                    diff *= 3.1 * pow(0.4 * exp_factor, 2.5) * exp(-exp_factor);
            }
            cv_dens += diff;
        }
        return cv_dens;
    };
}

std::function<double(double, double, double)> auxiliaries::phys::thermal_conductivity_FI(const std::function<double(double)> &rho, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_phi)
{
    return [=](double r, double t, double T)
    {
        using namespace constants::conversion;
        using namespace constants::scientific;
        double rho_cm3_over_g = rho(nbar_of_r(r)) / g_over_cm3_gev4,
               T_loc_8 = T / exp_phi(r) * (gev_over_k / 1E8);

        // scales, present in the calculations
        double rho_scale_cm3_over_g = 2.0E14,
               T_melt_8 = 2.4E5 * pow(rho_cm3_over_g, 1.0 / 3) * 1E-8;

        // units: erg / (cm * s * K) -> GeV^2
        double erg_over_cm_s_k_gev2 = erg_over_gev * gev_over_k / (gev_s * 1E-5 * km_gev);
        // In what follows we mostly refer to F & I, with the exception of some simplifications
        // introduced specifically for typical NS conditions.

        // (1) Quantum liquid region

        // usually requires to check if T < T_fp \approx mst_p \ge 0.1 M_N \sim 100 MeV, which as way above NS core temperatures

        // density must usually not exceed the one where nucleon matter become subdominant,
        // but due to consideration of either npemu or npeq matter, we shall omit this condition
        if (rho_scale_cm3_over_g <= rho_cm3_over_g)
        {
            return 1E23 * (rho_cm3_over_g / 1E14) / T_loc_8 * erg_over_cm_s_k_gev2;
        }

        // (2) Liquid metal region

        // usually requires to check if T < T_ep = p_fe = (ne/nsat)^{1/3} 1.7 fm^{-1} =
        // = [for most EoS Tmelt < T will be under crust, so ne \ge 1E-3 nB] \gapprox
        // (Ye)^{1/3} 1.68/5 GeV \gapprox 10^{0:2} MeV <- only wrong for extremely small Ye
        // which may occur at crust/exotic ph. However, my opinion that it's safe to
        // extend in general.

        // auxiliary variables; setup for taylor expansion
        double x = 0.08594 * log(rho_cm3_over_g) - 1.7949;

        std::vector<double> a = {-43.8644, -11.6758, 1.2698, 0.2798, 2.5652, -0.3879};
        double y = 0.0;

        for (size_t index = 0; index < a.size(); ++index)
        {
            y += a[index] * pow(x, 1.0 * index);
        }

        // rho_scale_cm3_over_g > rho_cm3_over_g now holds
        if (T_melt_8 < T_loc_8)
        {
            return 1.0 / (1.0 / (1E14 * pow(rho_cm3_over_g, 1.0 / 3) * T_loc_8) + exp(y) * T_loc_8) * erg_over_cm_s_k_gev2;
        }

        // (3) Solid region

        // auxiliary variables; setup for taylor expansions

        std::vector<double> b = {-41.1677, -7.8991, 3.4603, -0.8061};
        double z = 0.0;

        for (size_t index = 0; index < b.size(); ++index)
        {
            z += b[index] * pow(x, 1.0 * index);
        }

        double u = 3.6E-6 * pow(rho_cm3_over_g, 1.0 / 2) / T_loc_8;
        double s;
        if (u > 5)
        {
            s = Pi * Pi / (6 * u);
        }
        else
        {
            s = 1.0 - u / 4 + pow(u, 2.0) / 36 - pow(u, 4.0) / 3600 +
                pow(u, 6.0) / 211680 - pow(u, 8.0) / 10886400 + pow(u, 10.0) / 526901760;
        }

        // rho_scale_cm3_over_g > rho_cm3_over_g && T_melt_8 >= T_loc_8 now hold
        return 1.0 / (exp(z) * s + exp(y) * T_loc_8) * erg_over_cm_s_k_gev2;
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
    return std::min(Ts, Tb / exp_phi_at_R);
}

double auxiliaries::phys::critical_temperature_smeared_guassian(double k_fermi, double temp_ampl, double k_offs, double k_width, double quad_skew)
{
    return temp_ampl * exp(-pow((k_fermi - k_offs) / k_width, 2) - quad_skew * pow((k_fermi - k_offs) / k_width, 4));
}

double auxiliaries::phys::critical_temperature_double_lorenzian(double k_fermi, double t0, double k0, double k1, double k2, double k3)
{
    if (k_fermi <= k0 || k_fermi >= k2)
        return 0.0;
    return t0 * pow(k_fermi - k0, 2) / (pow(k_fermi - k0, 2) + pow(k1, 2)) * pow(k_fermi - k2, 2) / (pow(k_fermi - k2, 2) + pow(k3, 2));
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
    double fm_gev = 1.0E-18 * km_gev;
    switch (model)
    {
    case CriticalTemperatureModel::kGIPSF_NS:
        return critical_temperature_double_lorenzian(
            k_fermi, 8.8E-3 / 1.76, 0.18 / fm_gev, sqrt(0.1) / fm_gev, 1.2 / fm_gev, sqrt(0.6) / fm_gev);
    case CriticalTemperatureModel::kMSH_NS:
        return critical_temperature_double_lorenzian(
            k_fermi, 2.45E-3 / 1.76, 0.18 / fm_gev, sqrt(0.05) / fm_gev, 1.4 / fm_gev, sqrt(0.1) / fm_gev);
    case CriticalTemperatureModel::kAWP2_NS:
        return critical_temperature_double_lorenzian(
            k_fermi, 28E-3 / 1.76, 0.2 / fm_gev, sqrt(1.5) / fm_gev, 1.7 / fm_gev, sqrt(2.5) / fm_gev);
    case CriticalTemperatureModel::kSFB_NS:
        return critical_temperature_double_lorenzian(
            k_fermi, 45E-3 / 1.76, 0.1 / fm_gev, sqrt(4.5) / fm_gev, 1.55 / fm_gev, sqrt(2.5) / fm_gev);
    case CriticalTemperatureModel::kAO_NT:
        return critical_temperature_double_lorenzian(
            k_fermi, 4E-3 / 8.42, 1.2 / fm_gev, sqrt(0.45) / fm_gev, 3.3 / fm_gev, sqrt(5.0) / fm_gev);
    case CriticalTemperatureModel::kTTOA_NT:
        return critical_temperature_double_lorenzian(
            k_fermi, 2.1E-3 / 8.42, 1.1 / fm_gev, sqrt(0.6) / fm_gev, 3.2 / fm_gev, sqrt(2.4) / fm_gev);
    case CriticalTemperatureModel::kBEEHS_NT:
        return critical_temperature_double_lorenzian(
            k_fermi, 0.45E-3 / 8.42, 1.0 / fm_gev, sqrt(0.4) / fm_gev, 3.2 / fm_gev, sqrt(0.25) / fm_gev);
    case CriticalTemperatureModel::kTTAV_NT:
        return critical_temperature_double_lorenzian(
            k_fermi, 3.0E-3 / 8.42, 1.1 / fm_gev, sqrt(0.6) / fm_gev, 2.92 / fm_gev, sqrt(3.0) / fm_gev);
    case CriticalTemperatureModel::kA_NT:
        return critical_temperature_smeared_guassian(
            k_fermi, 1E9 / gev_over_k, 1.8 / fm_gev, 0.5 / fm_gev, 0.0);
    case CriticalTemperatureModel::kB_NT:
        return critical_temperature_smeared_guassian(
            k_fermi, 3E9 / gev_over_k, 2.0 / fm_gev, 0.5 / fm_gev, 0.0);
    case CriticalTemperatureModel::kC_NT:
        return critical_temperature_smeared_guassian(
            k_fermi, 1E10 / gev_over_k, 2.5 / fm_gev, 0.7 / fm_gev, 0.0);
    case CriticalTemperatureModel::kCCDK_PS:
        return critical_temperature_double_lorenzian(
            k_fermi, 102E-3 / 1.76, 0 / fm_gev, sqrt(9.0) / fm_gev, 1.3 / fm_gev, sqrt(1.5) / fm_gev);
    case CriticalTemperatureModel::kAO_PS:
        return critical_temperature_double_lorenzian(
            k_fermi, 14E-3 / 1.76, 0.15 / fm_gev, sqrt(0.22) / fm_gev, 1.05 / fm_gev, sqrt(3.8) / fm_gev);
    case CriticalTemperatureModel::kBS_PS:
        return critical_temperature_double_lorenzian(
            k_fermi, 17E-3 / 1.76, 0 / fm_gev, sqrt(2.9) / fm_gev, 0.8 / fm_gev, sqrt(0.08) / fm_gev);
    case CriticalTemperatureModel::kBCLL_PS:
        return critical_temperature_double_lorenzian(
            k_fermi, 1.69E-3 / 1.76, 0.05 / fm_gev, sqrt(0.07) / fm_gev, 1.05 / fm_gev, sqrt(0.16) / fm_gev);
    default:
        RHM_THROW(std::runtime_error, "Unknown critical temperature model.");
    }
}

auxiliaries::math::MatrixD::MatrixD(size_t rows, size_t cols, double aloc_value) : m_matrix(rows, std::vector<double>(cols, aloc_value)) {}

auxiliaries::math::MatrixD::MatrixD(const std::vector<std::vector<double>> &matrix) : m_matrix(matrix) {}

auxiliaries::math::MatrixD::MatrixD(const std::vector<double> &column)
{
    m_matrix = std::vector<std::vector<double>>(column.size(), std::vector<double>(1, 0.0));
    for (size_t i = 0; i < column.size(); i++)
    {
        m_matrix.at(i).at(0) = column.at(i);
    }
}

double &auxiliaries::math::MatrixD::at(size_t i, size_t j)
{
    return m_matrix.at(i).at(j);
}

const double &auxiliaries::math::MatrixD::at(size_t i, size_t j) const
{
    return m_matrix.at(i).at(j);
}

const std::vector<double> &auxiliaries::math::MatrixD::row(size_t i) const
{
    return m_matrix.at(i);
}

std::vector<double> auxiliaries::math::MatrixD::column(size_t j) const
{
    std::vector<double> result;
    for (size_t i = 0; i < this->rows(); i++)
    {
        result.push_back(this->at(i, j));
    }
    return result;
}

void auxiliaries::math::MatrixD::set_row(size_t i, const std::vector<double> &row)
{
    m_matrix.at(i) = row;
}

void auxiliaries::math::MatrixD::set_column(size_t j, const std::vector<double> &column)
{
    for (size_t i = 0; i < this->rows(); i++)
    {
        this->at(i, j) = column.at(i);
    }
}

size_t auxiliaries::math::MatrixD::rows() const
{
    return m_matrix.size();
}

size_t auxiliaries::math::MatrixD::columns() const
{
    return m_matrix.at(0).size();
}

double auxiliaries::math::MatrixD::det() const
{
    if (this->rows() != this->columns())
    {
        RHM_THROW(std::runtime_error, "Non-square matrix has no determinant.");
    }
    if (this->rows() == 1)
    {
        return at(0, 0);
    }
    double result = 1.0;
    // Gaussian elimination
    MatrixD top_diag = *this;
    for (size_t row = 0; row < this->rows(); ++row)
    {
        double max_entry = std::abs(top_diag.at(row, row));
        size_t max_entry_row = row;
        // Swap considered row with the one with abs max entry in this column
        for (size_t row2 = row + 1; row2 < this->rows(); ++row2)
        {
            if (std::abs(top_diag.at(row2, row)) > max_entry)
            {
                max_entry = std::abs(top_diag.at(row2, row));
                max_entry_row = row2;
            }
        }
        if (max_entry == 0.0)
        {
            return 0.0;
        }
        if (max_entry_row != row)
        {
            auto temp = top_diag.row(row);
            top_diag.set_row(row, top_diag.row(max_entry_row));
            top_diag.set_row(max_entry_row, temp);
            result *= -1.0;
        }
        for (size_t row2 = row + 1; row2 < this->rows(); ++row2)
        {
            double factor = top_diag.at(row2, row) / top_diag.at(row, row);
            for (size_t col = row; col < this->columns(); ++col)
            {
                top_diag.at(row2, col) -= factor * top_diag.at(row, col);
            }
        }
    }
    for (size_t row = 0; row < this->rows(); ++row)
    {
        result *= top_diag.at(row, row);
    }
    return result;
}

auxiliaries::math::MatrixD auxiliaries::math::MatrixD::transpose() const
{
    MatrixD result(this->columns(), this->rows());
    for (size_t i = 0; i < this->rows(); ++i)
    {
        for (size_t j = 0; j < this->columns(); ++j)
        {
            result.at(j, i) = this->at(i, j);
        }
    }
    return result;
}

auxiliaries::math::MatrixD auxiliaries::math::MatrixD::inverse() const
{
    if (this->rows() != this->columns())
    {
        RHM_THROW(std::runtime_error, "Non-square matrix has no inverse.");
    }
    MatrixD result(this->rows(), this->columns());
    // Gaussian elimination of augmented matrix
    MatrixD augmented = MatrixD(this->rows(), this->columns() * 2, 0.0);
    for (size_t row = 0; row < this->rows(); ++row)
    {
        for (size_t col = 0; col < this->columns(); ++col)
        {
            augmented.at(row, col) = this->at(row, col);
        }
        augmented.at(row, row + this->columns()) = 1.0;
    }
    for (size_t row = 0; row < this->rows(); ++row)
    {
        double max_entry = std::abs(augmented.at(row, row));
        size_t max_entry_row = row;
        // Swap considered row with the one with abs max entry in this column
        for (size_t row2 = row + 1; row2 < augmented.rows(); ++row2)
        {
            if (std::abs(augmented.at(row2, row)) > max_entry)
            {
                max_entry = std::abs(augmented.at(row2, row));
                max_entry_row = row2;
            }
        }
        if (max_entry == 0.0)
        {
            RHM_THROW(std::runtime_error, "Matrix is not invertible.");
        }
        if (max_entry_row != row)
        {
            auto temp = augmented.row(row);
            augmented.set_row(row, augmented.row(max_entry_row));
            augmented.set_row(max_entry_row, temp);
        }
        // Forward elimination
        for (size_t row2 = row + 1; row2 < augmented.rows(); ++row2)
        {
            double factor = augmented.at(row2, row) / augmented.at(row, row);
            for (size_t col = row; col < augmented.columns(); ++col)
            {
                augmented.at(row2, col) -= factor * augmented.at(row, col);
            }
        }
    }
    // Backward elimination
    for (size_t row = augmented.rows() - 1; row > 0; --row)
    {
        for (size_t row2 = row; row2 > 0; --row2)
        {
            // use row2 - 1 as a targeted row further to ensure unsigned underflow does not occur
            double factor = augmented.at(row2 - 1, row) / augmented.at(row, row);
            for (size_t col = row; col < augmented.columns(); ++col)
            {
                augmented.at(row2 - 1, col) -= factor * augmented.at(row, col);
            }
        }
    }
    // Obtain inverse
    for (size_t row = 0; row < this->rows(); ++row)
    {
        double factor = 1.0 / augmented.at(row, row);
        for (size_t col = 0; col < this->columns(); ++col)
        {
            result.at(row, col) = factor * augmented.at(row, col + this->columns());
        }
    }
    return result;
}

std::vector<double> auxiliaries::math::MatrixD::solve(const std::vector<double> &rhs) const
{
    if (this->rows() != this->columns())
    {
        RHM_THROW(std::runtime_error, "Non-square matrix cannot be subjected to solver.");
    }
    if (this->rows() != rhs.size())
    {
        RHM_THROW(std::runtime_error, "Matrix and right-hand side vector dimensions do not match.");
    }

    // Gaussian elimination of augmented matrix
    MatrixD augmented = MatrixD(this->rows(), this->columns() + 1, 0.0);
    for (size_t row = 0; row < this->rows(); ++row)
    {
        for (size_t col = 0; col < this->columns(); ++col)
        {
            augmented.at(row, col) = this->at(row, col);
        }
        augmented.at(row, this->columns()) = rhs[row];
    }
    for (size_t row = 0; row < this->rows(); ++row)
    {
        double max_entry = std::abs(augmented.at(row, row));
        size_t max_entry_row = row;
        // Swap considered row with the one with abs max entry in this column
        for (size_t row2 = row + 1; row2 < augmented.rows(); ++row2)
        {
            if (std::abs(augmented.at(row2, row)) > max_entry)
            {
                max_entry = std::abs(augmented.at(row2, row));
                max_entry_row = row2;
            }
        }
        if (max_entry == 0.0)
        {
            RHM_THROW(std::runtime_error, "Matrix is not invertible.");
        }
        if (max_entry_row != row)
        {
            auto temp = augmented.row(row);
            augmented.set_row(row, augmented.row(max_entry_row));
            augmented.set_row(max_entry_row, temp);
        }
        // Forward elimination
        for (size_t row2 = row + 1; row2 < augmented.rows(); ++row2)
        {
            double factor = augmented.at(row2, row) / augmented.at(row, row);
            for (size_t col = row; col < augmented.columns(); ++col)
            {
                augmented.at(row2, col) -= factor * augmented.at(row, col);
            }
        }
    }
    // Backsustitution
    std::vector<double> result(this->rows(), 0.0);
    for (size_t row = augmented.rows(); row > 0; --row)
    {
        // here row is forwarded one space ahead to not cause size_t overflow
        double row_sum = 0.0;
        for (size_t col = augmented.columns() - 2; col > row - 1; --col)
        {
            row_sum += augmented.at(row - 1, col) * result[col];
        }
        result[row - 1] = (augmented.at(row - 1, augmented.columns() - 1) - row_sum) / augmented.at(row - 1, row - 1);
    }
    return result;
}

auxiliaries::math::MatrixD auxiliaries::math::MatrixD::tridiagonal_inverse() const
{
    if (this->rows() != this->columns())
    {
        RHM_THROW(std::runtime_error, "Non-square matrix has no inverse.");
    }
    MatrixD result(this->rows(), this->columns());
    // Gaussian elimination of augmented matrix
    MatrixD augmented = MatrixD(this->rows(), this->columns() * 2, 0.0);
    for (size_t row = 0; row < this->rows(); ++row)
    {
        for (size_t col = 0; col < this->columns(); ++col)
        {
            augmented.at(row, col) = this->at(row, col);
        }
        if (augmented.at(row, row) == 0.0)
        {
            RHM_THROW(std::runtime_error, "Tridiagonal inverse cannot be applied with zeros on the diagonal.");
        }
        augmented.at(row, row + this->columns()) = 1.0;
    }
    for (size_t row = 0; row < this->rows() - 1; ++row)
    {
        // Forward elimination (only one row below to eliminate)
        double factor = augmented.at(row + 1, row) / augmented.at(row, row);
        for (size_t col = row; col < augmented.columns(); ++col)
        {
            augmented.at(row + 1, col) -= factor * augmented.at(row, col);
        }
    }
    for (size_t row = augmented.rows() - 1; row > 0; --row)
    {
        // Backward elimination (only one row above to eliminate)
        double factor = augmented.at(row - 1, row) / augmented.at(row, row);
        for (size_t col = row; col < augmented.columns(); ++col)
        {
            augmented.at(row - 1, col) -= factor * augmented.at(row, col);
        }
    }
    // Obtain inverse
    for (size_t row = 0; row < this->rows(); ++row)
    {
        double factor = 1.0 / augmented.at(row, row);
        for (size_t col = 0; col < this->columns(); ++col)
        {
            result.at(row, col) = factor * augmented.at(row, col + this->columns());
        }
    }
    return result;
}

std::vector<double> auxiliaries::math::MatrixD::tridiagonal_solve(const std::vector<double> &rhs) const
{
    if (this->rows() != this->columns())
    {
        RHM_THROW(std::runtime_error, "Non-square matrix cannot be subjected to tridiagonal solver.");
    }
    if (this->rows() != rhs.size())
    {
        RHM_THROW(std::runtime_error, "Matrix and right-hand side vector dimensions do not match.");
    }

    size_t n = this->rows();
    std::vector<double> v(n);     // The solution vector
    std::vector<double> c(n - 1); // new superdiagonal
    std::vector<double> g(n);     // new right hand side
    c[0] = this->at(0, 1) / this->at(0, 0);
    g[0] = rhs[0] / this->at(0, 0);
    for (size_t i = 1; i < n - 1; ++i)
    {
        c[i] = this->at(i, i + 1) / (this->at(i, i) - this->at(i, i - 1) * c[i - 1]);
        g[i] = (rhs[i] - this->at(i, i - 1) * g[i - 1]) / (this->at(i, i) - this->at(i, i - 1) * c[i - 1]);
    }
    g[n - 1] = (rhs[n - 1] - this->at(n - 1, n - 2) * g[n - 2]) / (this->at(n - 1, n - 1) - this->at(n - 1, n - 2) * c[n - 2]);
    v[n - 1] = g[n - 1];
    for (int i = n - 2; i >= 0; --i)
    {
        v[i] = g[i] - c[i] * v[i + 1];
    }
    return v;
}

auxiliaries::math::MatrixD auxiliaries::math::MatrixD::operator+(const MatrixD &other) const
{
    if (this->rows() != other.rows() || this->columns() != other.columns())
    {
        RHM_THROW(std::runtime_error, "Matrix dimensions do not match.");
    }
    MatrixD result(this->rows(), this->columns(), 0.0);
    for (size_t row = 0; row < this->rows(); ++row)
    {
        for (size_t col = 0; col < this->columns(); ++col)
        {
            result.at(row, col) = this->at(row, col) + other.at(row, col);
        }
    }
    return result;
}

auxiliaries::math::MatrixD auxiliaries::math::MatrixD::operator-(const MatrixD &other) const
{
    if (this->rows() != other.rows() || this->columns() != other.columns())
    {
        RHM_THROW(std::runtime_error, "Matrix dimensions do not match.");
    }
    MatrixD result(this->rows(), this->columns(), 0.0);
    for (size_t row = 0; row < this->rows(); ++row)
    {
        for (size_t col = 0; col < this->columns(); ++col)
        {
            result.at(row, col) = this->at(row, col) - other.at(row, col);
        }
    }
    return result;
}

auxiliaries::math::MatrixD auxiliaries::math::MatrixD::operator*(const MatrixD &other) const
{
    if (this->columns() != other.rows())
    {
        RHM_THROW(std::runtime_error, "Matrix dimensions do not match.");
    }
    MatrixD result(this->rows(), other.columns(), 0.0);
    for (size_t row = 0; row < this->rows(); ++row)
    {
        for (size_t col = 0; col < other.columns(); ++col)
        {
            for (size_t i = 0; i < this->columns(); ++i)
            {
                result.at(row, col) += this->at(row, i) * other.at(i, col);
            }
        }
    }
    return result;
}

auxiliaries::math::MatrixD auxiliaries::math::operator*(const MatrixD &matrix, double scalar)
{
    MatrixD result(matrix.rows(), matrix.columns(), 0.0);
    for (size_t row = 0; row < matrix.rows(); ++row)
    {
        for (size_t col = 0; col < matrix.columns(); ++col)
        {
            result.at(row, col) = scalar * matrix.at(row, col);
        }
    }
    return result;
}

auxiliaries::math::MatrixD auxiliaries::math::operator*(double scalar, const MatrixD &matrix)
{
    return matrix * scalar;
}

std::ostream &auxiliaries::math::operator<<(std::ostream &os, const MatrixD &matrix)
{
    for (size_t row = 0; row < matrix.rows(); ++row)
    {
        for (size_t col = 0; col < matrix.columns(); ++col)
        {
            os << matrix.at(row, col) << " ";
        }
        os << '\n';
    }
    return os;
}