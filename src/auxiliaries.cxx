
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
#include <ctime>
#include <iomanip>
#include <iostream>

std::vector<std::vector<double>> auxiliaries::io::read_tabulated_file(const std::string &path, std::pair<size_t, size_t> columns, std::pair<size_t, size_t> rows, double empty_value)
{
    std::ifstream fstr(path);
    if (!fstr.is_open())
        RHM_ERROR("Cannot open file " + path + ". ");
    std::vector<std::vector<double>> table;
    std::vector<std::string> relevant_lines;
    std::string line;
    size_t line_count = 0;
    while (std::getline(fstr, line))
    {
        // if the last line is reached, stop reading
        if (line_count >= rows.second && rows.second != 0)
            break;
        // if the first line is reached, while the last is not yet reached, memorize
        if (line_count >= rows.first)
            relevant_lines.push_back(line);

        ++line_count;
    }
    // handle row number special cases
    if (rows.second == 0)
        rows.second = relevant_lines.size() + rows.first;
    if (rows.first >= rows.second)
        RHM_ERROR("Invalid rows range (" + std::to_string(rows.first) + ", " + std::to_string(rows.second) + ") extracted from input file at " + path + ".");

    // handle column number special cases
    if (columns.second == 0)
    {
        std::string cleared_line = auxiliaries::io::retrieve_cleared_line(relevant_lines[0]);
        std::stringstream ss(cleared_line);
        std::string str;
        while (ss >> str)
            ++columns.second;
    }
    if (columns.first >= columns.second)
        RHM_ERROR("Invalid columns range (" + std::to_string(columns.first) + ", " + std::to_string(columns.second) + ") extracted from input file at " + path + ".");

    std::vector<std::vector<std::string>> str_data(rows.second - rows.first, std::vector<std::string>(columns.second - columns.first));
    for (size_t i = rows.first; i < rows.second; ++i)
    {
        std::string cleared_line = auxiliaries::io::retrieve_cleared_line(relevant_lines[i - rows.first]);
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
                RHM_ERROR("Intractable argument: " + std::string(e.what()) + ".");
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

auxiliaries::io::Logger::LogLevel auxiliaries::io::Logger::g_log_level;
std::ostream *auxiliaries::io::Logger::g_stream_ptr;

void auxiliaries::io::Logger::log(std::function<bool()> &&lazy_condition, LogLevel level, std::function<std::string()> &&lazy_message, std::string appendix) const
{
    static const std::map<auxiliaries::io::Logger::LogLevel, std::string> log_level_map = {
        {LogLevel::kTrace, "TRACE"},
        {LogLevel::kDebug, "DEBUG"},
        {LogLevel::kInfo, " INFO"},
        {LogLevel::kError, "ERROR"}};
    std::stringstream ss;
    if (level >= auxiliaries::io::Logger::g_log_level)
    {
        if (!lazy_condition())
            return;
        if (appendix != "")
            appendix.insert(0, ", ");
        // current time
        std::time_t t = std::time(nullptr);
        std::tm tm = *std::localtime(&t);

        ss << std::setfill('0') << std::setw(2) << tm.tm_hour << ":" << std::setfill('0') << std::setw(2) << tm.tm_min << ":" << std::setfill('0') << std::setw(2) << tm.tm_sec;
        ss << " [" << log_level_map.at(level) << "] <" << header << appendix << "> : " << lazy_message();
        *auxiliaries::io::Logger::g_stream_ptr << ss.str() << std::endl;
        // stream is flushed immediately. Otherwise stuck programs might not be able to communicate their issues
    }
    if (level == LogLevel::kError)
    {
        // write to cerr as well, if stream is not cerr
        if (auxiliaries::io::Logger::g_stream_ptr != &std::cerr)
            std::cerr << ss.str() << std::endl;
    }
}

void auxiliaries::math::Interpolator::instantiate(const std::vector<double> &input, const std::vector<double> &output, auxiliaries::math::Interpolator::InterpolationMode mode, bool extrapolate, bool enable_checks)
{
    if (enable_checks)
    {
        if (input.size() != output.size())
            RHM_ERROR("Input and output arrays have different sizes. Object name: " + m_obj_name +
                      ", input size: " + std::to_string(input.size()) + ", output size: " + std::to_string(output.size()));
    }
    m_input = input;
    m_mode = mode;
    m_extrapolate = extrapolate;
    m_enable_checks = enable_checks;
    // applicable to all modes, m_weights[0] shall be output values
    m_weights.push_back(output);
    switch (mode)
    {
    case auxiliaries::math::Interpolator::InterpolationMode::kLinear:
    {
        if (enable_checks)
            if (input.size() < 2)
                RHM_ERROR("Cannot perform linear interpolation with less than 2 points. Object name: " + m_obj_name +
                          ", input size: " + std::to_string(input.size()));
        std::vector<double> linear_coeffs(input.size() - 1);
        for (size_t i = 0; i < input.size() - 1; ++i)
        {
            linear_coeffs[i] = (output[i + 1] - output[i]) / (input[i + 1] - input[i]);
        }
        // 1+1 coefficients for linear interpolation
        m_weights.push_back(linear_coeffs);
        break;
    }
    // This interpolation mode became possible only thanks to https://signalsmith-audio.co.uk/writing/2021/monotonic-smooth-interpolation/
    case auxiliaries::math::Interpolator::InterpolationMode::kCubic:
    {
        if (enable_checks)
            if (input.size() < 2)
                RHM_ERROR("Cannot perform cubic interpolation with less than 2 points. Object name: " + m_obj_name +
                          ", input size: " + std::to_string(input.size()));
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
        // 1+3 coefficients for cubic interpolation
        m_weights.push_back(linear_coeffs);
        m_weights.push_back(quadratic_coeffs);
        m_weights.push_back(cubic_coeffs);
        break;
    }
    default:
        RHM_ERROR("Unknown interpolation mode.");
    }
    m_instantiated = true;
}

double auxiliaries::math::Interpolator::operator()(double arg) const
{
    if (m_enable_checks && !m_instantiated)
        RHM_ERROR("Uninitialized interpolator invoked. This error signifies a serious bug in the production. Object name: " + m_obj_name);
    bool decreasingly_sorted = (m_input.front() > m_input.back());
    if (m_enable_checks && !m_extrapolate)
        if ((!decreasingly_sorted && (arg < m_input.front() || arg > m_input.back())) || (decreasingly_sorted && (arg > m_input.front() || arg < m_input.back())))
        {
            std::stringstream ss;
            ss << std::scientific << std::setprecision(3) << "Searched value is out of the interpolator's range. Object name: " << m_obj_name << ", input range: [" << (decreasingly_sorted ? m_input.back() : m_input.front()) << ", " << (decreasingly_sorted ? m_input.front() : m_input.back()) << "], passed arg: " << arg;
            std::string msg = ss.str();
            RHM_ERROR(msg);
        }
    size_t low_pos;
    if (!decreasingly_sorted)
        low_pos = std::upper_bound(m_input.begin(), m_input.end(), arg) -
                  m_input.begin();
    else
        low_pos = std::upper_bound(m_input.begin(), m_input.end(), arg, std::greater<double>()) -
                  m_input.begin();
    // extrapolation: in case x is out of range, we extrapolate from the nearest polinomial
    if (m_extrapolate)
    {
        if (low_pos == m_input.size())
            --low_pos;
        if (low_pos == 0)
            ++low_pos;
    }
    // Conventional definition of low_pos : x is within [input[low_pos], input[low_pos+1]), if not extrapolating
    --low_pos;
    // if x is equal to input[i], return output[i]
    if (arg == m_input[low_pos])
        return m_weights[0][low_pos];

    // resolve shift from the passed arg
    double s = arg - m_input[low_pos];
    switch (m_mode)
    {
    case auxiliaries::math::Interpolator::InterpolationMode::kLinear:
    {
        return m_weights[0][low_pos] + m_weights[1][low_pos] * s;
    }
    case auxiliaries::math::Interpolator::InterpolationMode::kCubic:
    {
        return m_weights[0][low_pos] + m_weights[1][low_pos] * s +
               m_weights[2][low_pos] * s * s + m_weights[3][low_pos] * s * s * s;
    }
    default:
        RHM_ERROR("Unknown interpolation mode.");
    }
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

std::function<double(double, double, double)> auxiliaries::phys::thermal_conductivity_crust_Flowers_Itoh(const std::function<double(double)> &rho, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_phi)
{
    return [=](double r, double t, double T)
    {
        using namespace constants::conversion;
        using namespace constants::scientific;
        double rho_cm3_over_g = rho(nbar_of_r(r)) / g_over_cm3_gev4,
               T_loc_8 = T / exp_phi(r) * (gev_over_k / 1E8);

        // scales, present in the calculations
        double T_melt_8 = 2.4E5 * pow(rho_cm3_over_g, 1.0 / 3) * 1E-8;

        // units: erg / (cm * s * K) -> GeV^2
        double erg_over_cm_s_k_gev2 = erg_over_gev * gev_over_k / (gev_s * 1E-5 * km_gev);
        // Derivations of F & I, low density regions

        // (2) Liquid metal region

        // usually requires to check if T < T_ep = p_fe = (ne/nsat)^{1/3} 1.7 fm^{-1} =
        // = [for most EoS Tmelt < T will be under crust, so ne \ge 1E-3 nB] \gapprox
        // (Ye)^{1/3} 1.68/5 GeV \gapprox 10^{0:2} MeV <- only wrong for extremely small Ye

        // auxiliary variables; setup for taylor expansion
        double x = 0.08594 * log(rho_cm3_over_g) - 1.7949;

        std::vector<double> a = {-43.8644, -11.6758, 1.2698, 0.2798, 2.5652, -0.3879};
        double y = 0.0;

        for (size_t index = 0; index < a.size(); ++index)
        {
            y += a[index] * pow(x, 1.0 * index);
        }

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

        return 1.0 / (exp(z) * s + exp(y) * T_loc_8) * erg_over_cm_s_k_gev2;
    };
}

std::function<double(double, double, double)> auxiliaries::phys::thermal_conductivity_core_Flowers_Itoh(const std::function<double(double)> &rho, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_phi)
{
    return [=](double r, double t, double T)
    {
        using namespace constants::conversion;
        using namespace constants::scientific;
        double rho_cm3_over_g = rho(nbar_of_r(r)) / g_over_cm3_gev4,
               T_loc_8 = T / exp_phi(r) * (gev_over_k / 1E8);

        // units: erg / (cm * s * K) -> GeV^2
        double erg_over_cm_s_k_gev2 = erg_over_gev * gev_over_k / (gev_s * 1E-5 * km_gev);
        // Derivations of F & I, high density regions

        // (1) Quantum liquid region

        // usually requires to check if T < T_fp \approx mst_p \ge 0.1 M_N \sim 100 MeV, which is way above NS core temperatures

        // density must usually not exceed the one where nucleon matter become subdominant,
        // but due to consideration of either npemu or npeq matter, we shall omit this condition
        return 1E23 * (rho_cm3_over_g / 1E14) / T_loc_8 * erg_over_cm_s_k_gev2;
    };
}

std::function<double(double, double, double)> auxiliaries::phys::thermal_conductivity_core_Shternin_Yakovlev(
    const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
    const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
    const std::function<double(double)> &exp_phi, const std::function<double(double)> &superfluid_p_temp)
{
    // following Shternin & Yakovlev 2007, lepton carriers dominate thermal conductivity
    // Two considerations are possible: (1) only electrons, (2) electrons and muons

    // Electrons must be defined in the map
    if (!k_fermi_of_nbar.count(constants::species::electron))
    {
        return [](double, double, double)
        {
            return 0.0;
        };
    }

    return [=](double r, double t, double T)
    {
        using namespace constants::conversion;
        using namespace constants::scientific;
        using namespace constants::species;

        using auxiliaries::phys::Species;

        double nbar_val = nbar_of_r(r);
        double kf_e = k_fermi_of_nbar.at(electron)(nbar_val);
        double mst_e = m_stars_of_nbar.at(electron)(nbar_val);
        double T_loc = T / exp_phi(r);
        
        // if electrons are absent
        if (kf_e == 0)
            return 0.0;

        // decide if muons are present
        double kf_mu = 0;
        double mst_mu = muon.mass();
        if (k_fermi_of_nbar.count(muon))
        {
            kf_mu = k_fermi_of_nbar.at(muon)(nbar_val);
            mst_mu = m_stars_of_nbar.at(muon)(nbar_val);
        }

        // screening wave numbers
        double screening_pref = 4.0 / (137 * Pi);
        double q_t2 = 0, q_l2 = 0;
        for (auto it = m_stars_of_nbar.begin(); it != m_stars_of_nbar.end(); ++it)
        {
            auto key = it->first;
            if (key.qcharge() == 0 || key.classify() == Species::ParticleClassification::kMeson)
                continue;
            double kf = k_fermi_of_nbar.at(key)(nbar_val);
            double mst = m_stars_of_nbar.at(key)(nbar_val);
            q_t2 += screening_pref * kf * kf;
            q_l2 += screening_pref * kf * mst;
        }

        // way to calculate perpendicular frequencies
        auto perp_freq = [&](const Species &carrier, const Species &scatterer)
        {
            if (!k_fermi_of_nbar.count(carrier) || !k_fermi_of_nbar.count(scatterer))
                return 0.0;
            double kf_c = k_fermi_of_nbar.at(carrier)(nbar_val);
            double mst_c = m_stars_of_nbar.at(carrier)(nbar_val);
            double kf_s = k_fermi_of_nbar.at(scatterer)(nbar_val);
            double mst_s = m_stars_of_nbar.at(scatterer)(nbar_val);
            if (kf_c == 0 || kf_s == 0)
                return 0.0;
            return 24.0 * 1.202 / (137 * 137 * Pi * Pi * Pi) *
                   T_loc * pow(kf_s, 2.0) * kf_c / (mst_c * q_t2);
        };

        // way to calculate parallel frequencies
        auto paral_freq = [&](const Species &carrier, const Species &scatterer)
        {
            if (!k_fermi_of_nbar.count(carrier) || !k_fermi_of_nbar.count(scatterer))
                return 0.0;
            double kf_c = k_fermi_of_nbar.at(carrier)(nbar_val);
            double mst_c = m_stars_of_nbar.at(carrier)(nbar_val);
            double kf_s = k_fermi_of_nbar.at(scatterer)(nbar_val);
            double mst_s = m_stars_of_nbar.at(scatterer)(nbar_val);
            if (kf_c == 0 || kf_s == 0)
                return 0.0;
            return 4.0 * Pi * Pi / (137 * 137 * 5) *
                   pow(T_loc, 2) * pow(mst_s, 2.0) * mst_c / (kf_c * pow(q_l2, 1.5));
        };

        // way to calculate prime frequencies
        auto prime_freq = [&](const Species &carrier, const Species &scatterer)
        {
            if (!k_fermi_of_nbar.count(carrier) || !k_fermi_of_nbar.count(scatterer))
                return 0.0;
            double kf_c = k_fermi_of_nbar.at(carrier)(nbar_val);
            double mst_c = m_stars_of_nbar.at(carrier)(nbar_val);
            double kf_s = k_fermi_of_nbar.at(scatterer)(nbar_val);
            double mst_s = m_stars_of_nbar.at(scatterer)(nbar_val);
            if (kf_c == 0 || kf_s == 0)
                return 0.0;
            return 12.0 * 18.52 / (137 * 137 * Pi * Pi * Pi) *
                   pow(T_loc, 5.0 / 3) * pow(kf_s, 2.0) * mst_c / (pow(q_t2, 1.0 / 3) * kf_c * q_l2);
        };

        // superfluid factors

        // suppression for paralel electron-proton frequency
        auto R_paral_ep = [](double v)
        {
            // formula from the cited paper contained errors. This is the corrected fit
            return (0.998 + (2.04 + 0.68 * sqrt(v) + 5.7 * pow(v, 2.0) + 1.71 * pow(v, 4.0)) * exp(-1.04 * v)) *
                   exp(-sqrt(1.23 + v * v));
        };

        // suppression for perpendicular frequencies. r is (pfe2 + pfm2) / pfp2
        auto R_perp = [](double v, double r)
        {
            double p1 = 0.48 - 0.17 * r;
            double p3 = pow((1 - p1) * 45 * 1.202 / (4 * Pi * Pi * r), 2.0);
            return p1 * exp(-0.14 * v * v) + (1 - p1) / sqrt(1 + p3 * v * v);
        };

        // suppression for prime frequencies. r is (pfe2 + pfm2) / pfp2
        auto R_prime = [](double v, double r)
        {
            return pow(r + 1, 1.0 / 3) / pow(pow(r + 1, 2.0) - 0.757 * v + pow(0.50651 * v, 2.0), 1.0 / 6);
        };

        // get to calculating frequencies of collisions

        // electron-electron collisions
        double nu_ee_perp = perp_freq(electron, electron),
               nu_ee_paral = paral_freq(electron, electron);
        // electron-muon collisions
        double nu_em_perp = perp_freq(electron, muon),
               nu_em_paral = paral_freq(electron, muon),
               nu_em_prime = prime_freq(electron, muon);
        // electron-proton collisions
        double nu_ep_perp = perp_freq(electron, proton),
               nu_ep_paral = paral_freq(electron, proton);

        // muon-muon collisions
        double nu_mm_perp = perp_freq(muon, muon),
               nu_mm_paral = paral_freq(muon, muon);
        // muon-electron collisions
        double nu_me_perp = perp_freq(muon, electron),
               nu_me_paral = paral_freq(muon, electron),
               nu_me_prime = prime_freq(muon, electron);
        // muon-proton collisions
        double nu_mp_perp = perp_freq(muon, proton),
               nu_mp_paral = paral_freq(muon, proton);

        // resolve superfluidity
        double r_perp_total = 1,
               r_paral_ep = 1,
               r_prime = 1;
        double tau_p_inv = superfluid_p_temp(nbar_val) / T_loc;
        if (tau_p_inv > 1)
        {
            double proton_gap_scaled = superfluid_gap_1s0(1 / tau_p_inv);
            // pfp2 / (pfe2 + pfm2). Exercise extra care in quark matter, where it zeroes out
            double proton_ratio_inv = 0;
            double kf_p = 0;
            if (k_fermi_of_nbar.count(proton) && k_fermi_of_nbar.at(proton)(nbar_val) != 0)
            {
                kf_p = k_fermi_of_nbar.at(proton)(nbar_val);
                proton_ratio_inv = pow(kf_p, 2.0) / (pow(kf_e, 2.0) + pow(kf_mu, 2.0));
            }
            // If protons are absent (e.g. in quark matter), leave unmodified
            if (proton_ratio_inv != 0)
            {
                // if proton_ratio_inv is anomalously low (e.g.mixed phase?), r_perp_total grows linearly. We restrict this growth
                r_perp_total = std::min(R_perp(proton_gap_scaled, 1.0 / proton_ratio_inv), 1.0);
                // r_prime is restricted already even for low proton_ratio_inv
                r_prime = R_prime(proton_gap_scaled, 1.0 / proton_ratio_inv);
                r_paral_ep = R_paral_ep(proton_gap_scaled);
            }
        }

        double nu_ee = nu_ee_perp * r_perp_total + nu_ee_paral,
               nu_em = nu_em_perp * r_perp_total + nu_em_paral,
               nu_ep = nu_ep_perp * r_perp_total + nu_ep_paral * r_paral_ep,
               nu_mm = nu_mm_perp * r_perp_total + nu_mm_paral,
               nu_me = nu_me_perp * r_perp_total + nu_me_paral,
               nu_mp = nu_mp_perp * r_perp_total + nu_mp_paral;

        nu_em_prime *= r_prime;
        nu_me_prime *= r_prime;

        double nu_e = nu_ee + nu_em + nu_ep,
               nu_m = nu_mm + nu_me + nu_mp;

        // resolve relaxation times
        double tau_e = 1 / nu_e, tau_m = 0;
        double freq_denominator = nu_e * nu_m - nu_em_prime * nu_me_prime;
        if (kf_mu != 0 && freq_denominator != 0)
        {
            tau_e = (nu_m - nu_em_prime) / freq_denominator;
            tau_m = (nu_e - nu_me_prime) / freq_denominator;
        }

        return T_loc / 9 * (pow(kf_e, 3.0) * tau_e / mst_e + pow(kf_mu, 3.0) * tau_m / mst_mu);
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
        RHM_ERROR("Unknown critical temperature model.");
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
        RHM_ERROR("Non-square matrix has no determinant.");
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
        RHM_ERROR("Non-square matrix has no inverse.");
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
            RHM_ERROR("Matrix is not invertible.");
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
        RHM_ERROR("Non-square matrix cannot be subjected to solver.");
    }
    if (this->rows() != rhs.size())
    {
        RHM_ERROR("Matrix and right-hand side vector dimensions do not match.");
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
            RHM_ERROR("Matrix is not invertible.");
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
        RHM_ERROR("Non-square matrix has no inverse.");
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
            RHM_ERROR("Tridiagonal inverse cannot be applied with zeros on the diagonal.");
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
        RHM_ERROR("Non-square matrix cannot be subjected to tridiagonal solver.");
    }
    if (this->rows() != rhs.size())
    {
        RHM_ERROR("Matrix and right-hand side vector dimensions do not match.");
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
        RHM_ERROR("Matrix dimensions do not match.");
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
        RHM_ERROR("Matrix dimensions do not match.");
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
        RHM_ERROR("Matrix dimensions do not match.");
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