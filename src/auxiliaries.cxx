
#include "../include/auxiliaries.h"

#include <string>
#include <algorithm>
#include <vector>
#include <functional>
#include <fstream>

std::vector<double> auxiliaries::eos_data(const std::vector<double> &input, const std::function<std::vector<double>(const std::vector<double> &, std::ifstream &)> &eos, std::ifstream &fstr)
{
	std::vector<double> result; // empty placeholder for output array

	result = eos(input, fstr);

	fstr.clear();
	fstr.seekg(0, std::ios::beg); // reset fstream to initial state for further usage
	return result;				  // output
}

std::string auxiliaries::retrieve_cleared_line(const std::string &line)
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

double auxiliaries::interpolate(const std::vector<double> &input, const std::vector<double> &output,
                                auxiliaries::InterpolationMode mode, double x, bool disable_checks)
{
    if (!disable_checks)
    {
        if (input.size() != output.size())
            throw std::runtime_error("Input and output arrays have different sizes. Encountered in auxiliaries::interpolate");
    }
    // determine if input is increasing or decreasing
    bool decreasingly_sorted = (input[0] > input[1]) ? true : false;
    if (!disable_checks)
        if ((!decreasingly_sorted && (x < input.front() || x > input.back())) || (decreasingly_sorted && (x > input.front() || x < input.back())))
            throw std::runtime_error("Searched value is out of range. Encountered in auxiliaries::interpolate");

    // find index of x in input sorted array in reasonable time
    size_t low_pos;
    if (!decreasingly_sorted)
        low_pos = std::upper_bound(input.begin(), input.end(), x) -
                  input.begin();
    else
        low_pos = std::upper_bound(input.begin(), input.end(), x, std::greater<double>()) -
                  input.begin();
    --low_pos;
    // if x is equal to input[i], return output[i]
    if (x == input[low_pos])
        return output[low_pos];

    // interpolate
    switch (mode)
    {
    case auxiliaries::InterpolationMode::kLinear:
    {
        if (!disable_checks)
            if (input.size() < 2)
                throw std::runtime_error("Cannot perform linear interpolation with less than 2 points. Encountered in auxiliaries::interpolate");
        return output[low_pos] + (output[low_pos + 1] - output[low_pos]) * (x - input[low_pos]) / (input[low_pos + 1] - input[low_pos]);
    }
    break;
    case auxiliaries::InterpolationMode::kCubic:
    {
        if (!disable_checks)
            if (input.size() < 5)
                throw std::runtime_error("Cannot perform cubic interpolation with less than 5 points. Encountered in auxiliaries::interpolate");
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
        throw std::runtime_error("Unknown interpolation mode. Encountered in auxiliaries::interpolate");
    }
}

double auxiliaries::interpolate_cached(std::function<double(double)> &cache, const std::vector<double> &input, const std::vector<double> &output, auxiliaries::InterpolationMode mode, double x, bool disable_checks)
{
    if (!cache) // if a callable is not stored, cache one
    {
        if (!disable_checks)
        {
            if (input.size() != output.size())
                throw std::runtime_error("Input and output arrays have different sizes. Encountered in auxiliaries::interpolate");
        }
        // determine if input is increasing or decreasing
        bool decreasingly_sorted = (input[0] > input[1]) ? true : false;
        if (!disable_checks)
            if ((!decreasingly_sorted && (x < input.front() || x > input.back())) || (decreasingly_sorted && (x > input.front() || x < input.back())))
                throw std::runtime_error("Searched value is out of range. Encountered in auxiliaries::interpolate");

        // interpolate
        switch (mode)
        {
        case auxiliaries::InterpolationMode::kLinear:
        {
            if (!disable_checks)
                if (input.size() < 2)
                    throw std::runtime_error("Cannot perform linear interpolation with less than 2 points. Encountered in auxiliaries::interpolate");
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
                --low_pos;
                // if x is equal to input[i], return output[i]
                if (x == input[low_pos])
                    return output[low_pos];

                return output[low_pos] + (output[low_pos + 1] - output[low_pos]) * (x - input[low_pos]) / (input[low_pos + 1] - input[low_pos]);
            };
        }
        break;
        case auxiliaries::InterpolationMode::kCubic:
        {
            if (!disable_checks)
                if (input.size() < 5)
                    throw std::runtime_error("Cannot perform cubic interpolation with less than 5 points. Encountered in auxiliaries::interpolate");
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
            throw std::runtime_error("Unknown interpolation mode. Encountered in auxiliaries::interpolate");
        }
    }
    // determine if input is increasing or decreasing
    bool decreasingly_sorted = (input[0] > input[1]) ? true : false;
    if (!disable_checks)
        if ((!decreasingly_sorted && (x < input.front() || x > input.back())) || (decreasingly_sorted && (x > input.front() || x < input.back())))
            throw std::runtime_error("Searched value is out of range. Encountered in auxiliaries::interpolate");
    return cache(x);
}

auxiliaries::Species::Species(auxiliaries::Species::ParticleType type, auxiliaries::Species::ParticleClassification classification)
    : m_type(type), m_classification(classification)
{
}

bool auxiliaries::Species::operator==(const auxiliaries::Species& other) const
{
    return m_type == other.m_type;
}

bool auxiliaries::Species::operator<(const auxiliaries::Species &other) const
{
    return m_type < other.m_type;
}