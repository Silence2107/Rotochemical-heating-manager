#ifndef AUXILIARIES_H
#define AUXILIARIES_H

#include <string>
#include <vector>
#include <functional>
#include <map>
#include <iostream>

/// @brief various auxiliary functionality
namespace auxiliaries
{
    /// @brief IO auxiliary functionality.
    /// @brief Contains means to read tabulated files and clear lines from adverse symbols
    namespace io
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

        /// @brief Logger
        class Logger
        {
        private:
            std::string header;

        public:
            enum class LogLevel
            {
                kTrace,
                kDebug,
                kInfo,
                kError
            };
            /// @brief public constructor of Logger
            Logger(const std::string &header = "") : header{header} {}
            /// @brief globally accessible log level (set by the instantiator)
            static LogLevel g_log_level;
            /// @brief globally accessible log stream pointer (set by the instantiator)
            static std::ostream *g_stream_ptr;
            /// @brief logger functionality
            /// @param lazy_condition logging lambda condition, which must evaluate to true to log
            /// @param level log level
            /// @param lazy_message lambda, that must evaluate to the log message
            /// @param appendix header appendix for flexibility
            void log(std::function<bool()> &&lazy_condition, LogLevel level, std::function<std::string()> &&lazy_message, std::string appendix = "") const;
        };

        /// @brief Syntaxic sugar to log an error and tell the compiler we are exiting.
        /// @brief Note it is not actually enclosed by the namespace!
        /// @param message log message
#define RHM_ERROR(message)                                                                                   \
    do                                                                                                       \
    {                                                                                                        \
        auxiliaries::io::Logger logger(std::string(__FILE__) +                                               \
                                       "//" + std::string(__func__) + ", line " + std::to_string(__LINE__)); \
        logger.log([]()                                                                                      \
                   { return true; }, auxiliaries::io::Logger::LogLevel::kError, [=]()                        \
                   { return message; });                                                                     \
        *auxiliaries::io::Logger::g_stream_ptr << "Terminating on error.\n";                                 \
        std::exit(EXIT_FAILURE);                                                                             \
    } while (false)

    }

    /// @brief Math auxiliary functionality.
    /// @brief Includes interpolation, integration, caching functionality and matrix operations
    namespace math
    {
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
            /// @brief Monotonic cubic interpolation; requires at least two points to interpolate
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

        /// @brief Prototype for cached interpolator
        using CachedInterpolatorWrap = CachedFunc<std::function<double(double)>, double, const std::vector<double> &, const std::vector<double> &, InterpolationMode, double, bool, bool>;

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
        /// @cite Quadrature coefficients - https://pomax.github.io/bezierinfo/legendre-gauss.html
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
                    RHM_ERROR("Unimplemented integration mode.");
                }
            };
        }

        /// @brief Matrix number class
        class MatrixD
        {
        private:
            /// @brief matrix elements
            std::vector<std::vector<double>> m_matrix;

        public:
            /// @brief size + allocation constructor
            /// @param rows number of rows
            /// @param columns number of columns
            /// @param aloc_value value to populate matrix with
            MatrixD(size_t rows, size_t cols, double aloc_value = 0.0);
            /// @brief direct constructor
            /// @param matrix matrix elements
            MatrixD(const std::vector<std::vector<double>> &matrix);
            /// @brief vector-column constructor
            /// @param column vector-column
            MatrixD(const std::vector<double> &column);

            /// @brief matrix element access
            /// @param i row index
            /// @param j column index
            /// @return matrix element
            double &at(size_t i, size_t j);

            /// @brief matrix element access
            /// @param i row index
            /// @param j column index
            /// @return matrix element
            const double &at(size_t i, size_t j) const;

            /// @brief matrix row access
            /// @param i row index
            /// @return matrix row
            const std::vector<double> &row(size_t i) const;

            /// @brief matrix column access
            /// @param j column index
            /// @return matrix column
            std::vector<double> column(size_t j) const;

            /// @brief matrix row setter
            /// @param i row index
            /// @param row new row elements
            void set_row(size_t i, const std::vector<double> &row);

            /// @brief matrix column setter
            /// @param j column index
            /// @param column new column elements
            void set_column(size_t j, const std::vector<double> &column);

            /// @brief matrix rows number
            /// @return matrix rows number
            size_t rows() const;

            /// @brief matrix columns number
            /// @return matrix columns number
            size_t columns() const;

            /// @brief matrix determinant
            /// @return matrix determinant
            double det() const;

            /// @brief matrix transpose
            /// @return matrix transpose
            MatrixD transpose() const;

            /// @brief matrix inverse
            /// @return matrix inverse
            MatrixD inverse() const;

            /// @brief solve matrix equation M*X = RHS
            /// @param rhs right hand side (RHS)
            /// @return solution (X)
            std::vector<double> solve(const std::vector<double> &rhs) const;

            /// @brief tridiagonal matrix inverse (i.e. diagonal and 1-off-diagonal elements only, check is not performed).
            /// @brief Only applicable for tridiagonal matrices with all diagonal elements non-zero.
            /// @return tridiagonal matrix inverse
            MatrixD tridiagonal_inverse() const;

            /// @brief solve tridiagonal matrix equation M*X = RHS (i.e. diagonal and 1-off-diagonal elements only, check is not performed).
            /// @brief Only applicable for tridiagonal matrices with all diagonal elements non-zero.
            /// @param rhs right hand side (RHS)
            /// @return solution (X)
            std::vector<double> tridiagonal_solve(const std::vector<double> &rhs) const;

            /// @brief matrix addition
            /// @param other matrix to add
            /// @return matrix sum
            MatrixD operator+(const MatrixD &other) const;

            /// @brief matrix subtraction
            /// @param other matrix to subtract
            /// @return matrix difference
            MatrixD operator-(const MatrixD &other) const;

            /// @brief matrix multiplication
            /// @param other matrix to multiply with
            /// @return matrix product
            MatrixD operator*(const MatrixD &other) const;
        };

        /// @brief matrix multiplication with scalar
        /// @param scalar scalar to multiply with
        /// @param matrix matrix to multiply with
        /// @return matrix product
        MatrixD operator*(double scalar, const MatrixD &matrix);

        /// @brief matrix multiplication with scalar
        /// @param matrix matrix to multiply with
        /// @param scalar scalar to multiply with
        /// @return matrix product
        MatrixD operator*(const MatrixD &matrix, double scalar);

        /// @brief print matrix to stdout
        std::ostream &operator<<(std::ostream &os, const MatrixD &matrix);
    }

    /// @brief Physics related functionality.
    /// @brief Includes a handler for particle species, critical phenomena related functions, as well as Te(Tb) relation and thermodynamic quantities
    namespace phys
    {
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
            /// @brief string namings for species
            std::string m_name;
            /// @brief particle mass
            double m_mass;
            /// @brief particle charge
            double m_qcharge;
            /// @brief particle baryon number
            double m_bcharge;

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
            Species(ParticleType naming, ParticleClassification classification, const std::string &name, double mass, double qcharge, double bcharge);
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
            {
                return m_type;
            }
            /// @brief Get particle classification
            ParticleClassification classify() const
            {
                return m_classification;
            }
            /// @brief Get particle name
            std::string name() const
            {
                return m_name;
            }
            /// @brief Get particle mass
            double mass() const
            {
                return m_mass;
            }
            /// @brief Get particle charge
            double qcharge() const
            {
                return m_qcharge;
            }
            /// @brief Get particle baryon number
            double bcharge() const
            {
                return m_bcharge;
            }
        };

        /// @brief Specific heat of substance, based on Fermi gas model
        /// @param k_fermi_of_nbar fermi momentum [GeV] of species as a function of baryon density [GeV^3]
        /// @param m_stars_of_nbar mass of stars [GeV] of species as a function of baryon density [GeV^3]
        /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
        /// @param nbar_sf_shift lowest baryon density [GeV^3] with triplet pairing
        /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
        /// @param superfluid_p_temp proton SF critical temperature [GeV] as a function of baryon density [GeV^3]
        /// @param superfluid_n_temp neutron SF critical temperature [GeV] as a function of baryon density [GeV^3]
        /// @param superconduct_q_gap quark superconductivity gap [GeV] as a function of baryon density [GeV^3]
        /// @return specific heat as a function of radius, time and T^inf [natural units]
        /// @cite Base density - Yanagi, 2020; superfluidity - Yakovlev, Levenfish, 1999; superconductivity - Blaschke, Grigorian 2001
        std::function<double(double, double, double)> fermi_specific_heat_density(
            const std::map<auxiliaries::phys::Species, std::function<double(double)>> &k_fermi_of_nbar,
            const std::map<auxiliaries::phys::Species, std::function<double(double)>> &m_stars_of_nbar, const std::function<double(double)> &nbar_of_r,
            double nbar_sf_shift, const std::function<double(double)> &exp_phi,
            const std::function<double(double)> &superfluid_p_temp, const std::function<double(double)> &superfluid_n_temp, const std::function<double(double)> &superconduct_q_gap);

        /// @brief Thermal conductivity of substance
        /// @param rho energy density of the substance [GeV^4] as a function of baryon density [GeV^3]
        /// @param nbar_of_r baryon density [GeV^3] as a function of radius [GeV^{-1}]
        /// @param exp_phi e^phi metric function of radius [GeV^{-1}]
        /// @return thermal conductivity as a function of radius, time and T^inf [natural units]
        /// @cite Base value - Flowers, Itoh, 1982
        std::function<double(double, double, double)> thermal_conductivity_FI(const std::function<double(double)> &rho, const std::function<double(double)> &nbar_of_r, const std::function<double(double)> &exp_phi);

        /// @brief Te-Tb relation
        /// @param Tb internal temperature [GeV], measured by distant observer (inf)
        /// @param R NS radius [GeV^{-1}]
        /// @param M NS mass [GeV]
        /// @param eta NS light element share
        /// @return corresponding surface temperature [GeV]
        /// @cite Te(Tb) - Potekhin, Chabrier, 1997
        double te_tb_relation(double Tb, double R, double M, double eta);

        /// @brief Nucleon critical temperature parametrization
        /// @param k_fermi fermi momentum [GeV], correspoding to the nucleon
        /// @param temp_ampl amplitude parameter [GeV]
        /// @param k_offs mean-like parameter [GeV]
        /// @param k_width variance-like parameter [GeV]
        /// @param quad_skew high order correction to the gaussian [dimensionless]
        /// @return T0 * exp[- ((k_fermi - k_offs) / (k_width))^2 - quad_skew * ((k_fermi - k_offs) / (k_width))^4]
        /// @cite Parametrization - Yanagi, Nagata, 2019
        double critical_temperature_smeared_guassian(double k_fermi, double temp_ampl, double k_offs, double k_width, double quad_skew);

        /// @brief Nucleon critical temperature parametrization
        /// @param k_fermi fermi momentum [GeV], correspoding to the nucleon
        /// @param t0 amplitude parameter [GeV]
        /// @param k0 fit parameter [GeV]
        /// @param k1 fit parameter [GeV]
        /// @param k2 fit parameter [GeV]
        /// @param k3 fit parameter [GeV]
        /// @return T0 * (k_fermi - k0)^2/((k_fermi - k0)^2 + k1^2) * (k_fermi - k2)^2/((k_fermi - k2)^2 + k3^2) * I(k0 < k_fermi < k2)
        /// @cite Parametrization - Kaminker, Haensel, Yakovlev 2001
        double critical_temperature_double_lorenzian(double k_fermi, double t0, double k0, double k1, double k2, double k3);

        /// @brief Enumerate choices for the nucleon critical temperature parametrization.
        /// @brief Naming convention: k{name}_{baryon}{state}
        enum class CriticalTemperatureModel
        {
            kGIPSF_NS,
            kMSH_NS,
            kAWP2_NS,
            kSFB_NS,
            kAO_NT,
            kTTOA_NT,
            kBEEHS_NT,
            kTTAV_NT,
            kA_NT,
            kB_NT,
            kC_NT,
            kCCDK_PS,
            kAO_PS,
            kBS_PS,
            kBCLL_PS
        };

        /// @brief Critical temperature in given model
        /// @param k_fermi fermi momentum [GeV], correspoding to the nucleon
        /// @param model choice of the model
        /// @cite Fit parameters - Ho, Elshamouty, 2015
        double critical_temperature(double k_fermi, CriticalTemperatureModel model);

        /// @brief Superfluid gap in 1S0 state (A)
        /// @param tau T/Tc
        /// @cite Parametrization - Yakovlev, Levenfish, 1999
        double superfluid_gap_1s0(double tau);

        /// @brief Superfluid gap in 3P2 state (B)
        /// @param tau T/Tc
        /// @cite Parametrization - Yakovlev, Levenfish, 1999
        double superfluid_gap_3p2(double tau);
    }
}
#endif