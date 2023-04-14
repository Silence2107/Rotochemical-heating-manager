#ifndef AUXILIARIES_H
#define AUXILIARIES_H

#include <string>
#include <vector>
#include <functional>

#include "../include/constants.h"

/// @brief various auxiliary functionality;
/// contains:
/// function retrieve_cleared_line;
/// class CachedFunc;
namespace auxiliaries
{
    /// @brief function that cleares line from stream
    /// @param line line to clear out
    /// @return line, cleared of auxiliary symbols, preceeding, trailing or excessive whitespaces; separation between words is done with whitespaces
    std::string retrieve_cleared_line(const std::string &line);

    /// @brief Allows to wrap std::function of general form so that to support caching.
    /// @tparam Cache  Cache type;
    /// Cache is supposed to be a first parameter of a function (non-const reference!);
    /// User is supposed to write cache usage inside the function manually;
    /// After wrapping, cache is fully handled by the class, with possible changes at each call
    /// @tparam Foutput Function's output type
    /// @tparam ...Args Function arguments besides cache
    /// You are supposed to put them after Cache&amp; parameter
    template <class Cache, class Foutput, class... Args>
    class CachedFunc
    {
    private:
        std::function<Foutput(Cache &, Args...)> func; // callable that may use caching
        Cache cache;                                   // data which should be saved
    public:
        /// @brief public constructor of CachedFunc
        CachedFunc(const std::function<Foutput(Cache &, Args...)> &func, Cache defaultcache = Cache()) : func{func} 
        {
            erase(defaultcache);
        }
        /// @brief Invokes wrapped function based on cache and provided arguments
        /// @param ...args Function arguments besides cache
        /// @return Foutput from this invocation
        Foutput operator()(Args... args)
        {
            return func(cache, args...);
        }
        /// @brief Returns wrapped function as std::function, with cache running under the hood
        /// @return std::function<Foutput(Args...)>
        operator std::function<Foutput(Args...)>()
        {
            return [this](Args... args) { return func(this->cache, args...); };
        }
        /// @brief Erases cache
        /// @param defaultval Value that is substituted inside the cache; equals 'Cache()' by default, pass your version if needed/required
        void erase(Cache defaultval = Cache())
        {
            cache = defaultval;
        }
    };

    /// @brief Interpolation modes for arrays
    enum class InterpolationMode
    {
        /// @brief Linear interpolation; requires at least two points to interpolate
        kLinear,
        /// @brief Cubic interpolation (spline); requires at least five points to interpolate
        kCubic
    };

    /// @brief function that interpolates array based on mode provided;
    /// Make sure to have sufficient amount of points when choosing mode;
    /// Also make sure to have arrays sorted beforehand. Sorting will not be done here;
    /// Note : both increasing and decreasing input arrays are supported
    /// @param input input array to interpolate
    /// @param output output array to interpolate
    /// @param mode interpolation mode
    /// @param x x-coordinate of point to interpolate
    /// @param disable_checks false by default. If true, no checks are performed. For performance reasons.
    /// @return interpolated value
    double interpolate(const std::vector<double> &input, const std::vector<double> &output, InterpolationMode mode, double x, bool disable_checks = false);

    /// @brief function that interpolates array based on mode provided;
    /// Make sure to have sufficient amount of points when choosing mode;
    /// Also make sure to have arrays sorted beforehand. Sorting will not be done here;
    /// Note : both increasing and decreasing input arrays are supported
    /// @param cache Cache support. Wrap it with CachedFunc to use it. 
    /// @param input input array to interpolate
    /// @param output output array to interpolate
    /// @param mode interpolation mode
    /// @param x x-coordinate of point to interpolate
    /// @param disable_checks false by default. If true, no checks are performed. For performance reasons.
    /// @return interpolated value
    double interpolate_cached(std::function<double(double)> &cache, const std::vector<double> &input, const std::vector<double> &output, InterpolationMode mode, double x, bool disable_checks = false);

    /// @brief integrator over NS volume
    /// @param function to integrate with parameters double r, args...
    /// @param rmin lower bound
    /// @param rmax upper bound
    /// @param r_step integration step
    /// @param exp_lambda exponent of lambda function(to be used in jacobian)
    /// @return integral value as a function of args...
    template <class... Args>
    std::function<double(Args... args)> integrate_volume(const std::function<double(double, Args...)> &function, double rmin, double rmax, double r_step, const std::function<double(double)> &exp_lambda)
    {
        return [=](Args... args) -> double
        {
            /*double result = 0.0;
            for (double r = rmin + r_step / 2; r < rmax - r_step / 2; r += r_step)
            {
                result += function(r, args...) * 4 * constants::scientific::Pi * r * r * r_step * exp_lambda(r);
            }
            return result;*/
            // 6-point Gauss-Legendre quadrature
            std::vector<double> weights = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704},
                                 points = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521};
            double result = 0.0;
            auto integrand = [=](double r) -> double
            {
                return function(r, args...) * 4 * constants::scientific::Pi * r * r * exp_lambda(r);
            };
            for (int i = 0; i < weights.size(); i++)
            {
                result += (rmax - rmin) / 2 * weights[i] * integrand((rmax - rmin) / 2 * points[i] + (rmax + rmin) / 2);
            }
            return result;
        };
    }
}

#endif