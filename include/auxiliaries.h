#ifndef AUXILIARIES_H
#define AUXILIARIES_H

#include <string>
#include <vector>
#include <functional>

/// @brief various auxiliary functionality
namespace auxiliaries
{
    /// @brief function that reads tabulated file and organizes it into matrix data
    /// @param path path to file to be read
    /// @param columns column slice of signature tuple[first, last) to be read; leave default to self-deduce from the first line
    /// @param rows row slice of signature tuple[first, last) to be read; leave default to read all rows
    /// @param empty_value value to be assigned to empty cells; defaults to 0
    /// @return a vector of vectors, where each vector is a column of the table; absent values are represented by empty_value
    std::vector<std::vector<double>> read_tabulated_file(const std::string &path, std::pair<size_t, size_t> columns = {0, 0}, std::pair<size_t, size_t> rows = {0, 0}, double empty_value = 0);

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
            return [this](Args... args)
            { return func(this->cache, args...); };
        }
        /// @brief Erases cache
        /// @param defaultval Value that is substituted inside the cache; equals 'Cache()' by default, pass your version if needed/required
        void erase(Cache defaultval = Cache())
        {
            cache = defaultval;
        }
    };

    /// @brief Deduction guide for CachedFunc; couldn't put it to work so far
    // template <class Cache, class Foutput, class... Args>
    // CachedFunc(const std::function<Foutput(Cache &, Args...)> &func) -> CachedFunc<Cache, Foutput, Args...>;

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
    /// @param extrapolate false by default; if true, extrapolation beyond input bounds is allowed
    /// @param enable_checks true by default. If false, no checks are performed. Performance gain is not significant, but present.
    /// @return interpolated value
    double interpolate(const std::vector<double> &input, const std::vector<double> &output, InterpolationMode mode, double x, bool extrapolate = false, bool enable_checks = true);

    /// @brief function that interpolates array based on mode provided;
    /// Make sure to have sufficient amount of points when choosing mode;
    /// Also make sure to have arrays sorted beforehand. Sorting will not be done here;
    /// Note : both increasing and decreasing input arrays are supported
    /// @param cache Cache support. Wrap it with CachedFunc to use it.
    /// @param input input array to interpolate
    /// @param output output array to interpolate
    /// @param mode interpolation mode
    /// @param x x-coordinate of point to interpolate
    /// @param extrapolate false by default; if true, extrapolation beyond input bounds is allowed
    /// @param enable_checks true by default. If false, no checks are performed. Performance gain is not significant, but present.
    /// @return interpolated value
    double interpolate_cached(std::function<double(double)> &cache, const std::vector<double> &input, const std::vector<double> &output, InterpolationMode mode, double x, bool extrapolate = false, bool enable_checks = true);

    /// @brief Integration modes
    enum class IntegrationMode
    {
        /// @brief Rectangular integration
        kRectangular,
        /// @brief Gauss-Legendre integration with 6 points
        kGaussLegendre_6p,
        /// @brief Gauss-Legendre integration with 12 points
        kGaussLegendre_12p
    };

    /// @brief integrator over NS volume
    /// @param function to integrate with parameters double r, args...
    /// @param rmin lower bound
    /// @param rmax upper bound
    /// @param exp_lambda exponent of lambda function(to be used in jacobian)
    /// @param mode integration mode
    /// @param r_step integration step (for rectangular integration). Defaults to 10 meters
    /// @return integral value as a function of args...
    template <class... Args>
    std::function<double(Args... args)> integrate_volume(const std::function<double(double, Args...)> &function, double rmin, double rmax, const std::function<double(double)> &exp_lambda, IntegrationMode mode, double r_step = 5.06E16)
    {
        const double pi = 3.14159265359;
        return [=](Args... args) -> double
        {
            switch (mode)
            {
            case IntegrationMode::kRectangular:
            {
                double result = 0.0;
                for (double r = rmin + r_step / 2; r < rmax - r_step / 2; r += r_step)
                {
                    result += function(r, args...) * 4 * pi * r * r * r_step * exp_lambda(r);
                }
                return result;
            }
            case IntegrationMode::kGaussLegendre_6p:
            {
                std::vector<double> weights = {0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910, 0.3607615730481386, 0.1713244923791704},
                                    points = {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969, 0.2386191860831969, 0.6612093864662645, 0.9324695142031521};
                double result = 0.0;
                auto integrand = [=](double r) -> double
                {
                    return function(r, args...) * 4 * pi * r * r * exp_lambda(r);
                };
                for (size_t i = 0; i < weights.size(); i++)
                {
                    result += (rmax - rmin) / 2 * weights[i] * integrand((rmax - rmin) / 2 * points[i] + (rmax + rmin) / 2);
                }
                return result;
            }
            case IntegrationMode::kGaussLegendre_12p:
            {
                std::vector<double> weights = {0.04717533638651177, 0.10693932599531843, 0.16007832854334622, 0.20316742672306592, 0.23349253653835481, 0.24914704581340277, 0.24914704581340277, 0.23349253653835481, 0.20316742672306592, 0.16007832854334622, 0.10693932599531843, 0.04717533638651177},
                                    points = {-0.9815606342467192, -0.9041172563704749, -0.7699026741943047, -0.5873179542866175, -0.3678314989981802, -0.1252334085114699, 0.1252334085114699, 0.3678314989981802, 0.5873179542866175, 0.7699026741943047, 0.9041172563704749, 0.9815606342467192};
                double result = 0.0;
                auto integrand = [=](double r) -> double
                {
                    return function(r, args...) * 4 * pi * r * r * exp_lambda(r);
                };
                for (size_t i = 0; i < weights.size(); i++)
                {
                    result += (rmax - rmin) / 2 * weights[i] * integrand((rmax - rmin) / 2 * points[i] + (rmax + rmin) / 2);
                }
                return result;
            }
            default:
                throw std::runtime_error("Unimplemented integration mode; Encountered in auxiliaries::integrate_volume");
            }
        };
    }

    /// @brief Class that handles species-specific namings
    class Species
    {
    public:
        enum class ParticleType;
        enum class ParticleClassification;

    private:
        /// @brief namings for species
        ParticleType m_type;
        /// @brief species classification
        ParticleClassification m_classification;

    public:
        /// @brief Expected particle types
        enum class ParticleType
        {
            kElectron,
            kMuon,
            kTau,
            kProton,
            kNeutron,
            kPion,
            kKaon,
            kUquark,
            kDquark,
            kSquark
        };
        /// @brief Possible particle classifications
        enum class ParticleClassification
        {
            kLepton,
            kBaryon,
            kMeson,
            kQuark
        };
        /// @brief Constructor
        /// @param naming particle type
        /// @param classification particle classification
        Species(ParticleType naming, ParticleClassification classification);
        /// @brief Comparison operator
        /// @param other other species
        /// @return true if species are the same
        bool operator==(const Species &other) const;
        /// @brief Comparison operator (for lookup in std::map)
        /// @param other other species
        /// @return true if particle types are ordered as provided in ParticleType enum
        bool operator<(const Species &other) const;
        /// @brief Get particle type
        ParticleType identify() const
        { return m_type; }
        /// @brief Get particle classification
        ParticleClassification classify() const
        { return m_classification; }
    };
}
#endif