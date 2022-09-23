
#include "../include/auxiliaries.h"

#include <string>
#include <algorithm>
#include <vector>
#include <functional>

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
        if (input.empty())
            throw std::runtime_error("Input array is empty. Encountered in auxiliaries::interpolate");
    }
    // determine if input is increasing or decreasing
    bool decreasingly_sorted = (input[0] > input[1]) ? true : false;
    if (!disable_checks)
        if ( (!decreasingly_sorted && (x < input.front() || x > input.back()) ) || (decreasingly_sorted && (x > input.front() || x < input.back()) ) )
            throw std::runtime_error("Searched value is out of range. Encountered in auxiliaries::interpolate");

    // find index of x in input sorted array in reasonable time
    size_t i;
    if (!decreasingly_sorted)
        i = std::upper_bound(input.begin(), input.end(), x) -
            input.begin();
    else
        i = std::upper_bound(input.begin(), input.end(), x, std::greater<double>()) -
            input.begin();
    --i;
    // if x is equal to input[i], return output[i]
    if (x == input[i])
        return output[i];

    // interpolate
    switch (mode)
    {
    case auxiliaries::InterpolationMode::kLinear:
        return output[i] + (output[i + 1] - output[i]) * (x - input[i]) / (input[i + 1] - input[i]);
        break;
    case auxiliaries::InterpolationMode::kCubic:
        // throw non-implemented exception
        throw std::runtime_error("Cubic interpolation is not implemented yet. Encountered in auxiliaries::interpolate");
        /*if (i == 0)
            return output[i] + (output[i + 1] - output[i]) * (x - input[i]) / (input[i + 1] - input[i]);
        else if (i == input.size() - 1)
            return output[i] + (output[i] - output[i - 1]) * (x - input[i]) / (input[i] - input[i - 1]);
        else
            return output[i] + (output[i + 1] - output[i]) * (x - input[i]) / (input[i + 1] - input[i]) +
                   (output[i] - output[i - 1]) * (x - input[i]) / (input[i] - input[i - 1]);*/
        break;
    default:
        throw std::runtime_error("Unknown interpolation mode. Encountered in auxiliaries::interpolate");
    }
}