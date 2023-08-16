#include "../../include/auxiliaries.h"
#include "../../include/cooling.h"
#include "../../include/constants.h"
#include "../../include/tov_solver.h"
#include "../../include/inputfile.hpp"

#include <vector>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <fstream>

int main()
{
    using namespace inputfile;

    // RUN --------------------------------------------------------------------------

    // EoS definition

    // returns {r, m} pair at given center density. (Hopefully) cleans up all global cache that may spoil further calls
    auto get_m_r_at_density = [&](double edensity)
    {
        auto eos_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, double, double>(
            [&](std::vector<std::vector<double>> &cache, double rho)
            {
                if (rho < 0 || rho > edensity_upp)
                    THROW(std::runtime_error, "Data request out of range.");
                if (rho <= edensity_low)
                    return 0.0;
                if (cache.empty() || cache[0].size() != discr_size_EoS)
                {                                                                                        // then fill/refill cache
                    cache = std::vector<std::vector<double>>(2, std::vector<double>(discr_size_EoS, 0)); // initialize 2xdiscr_size_EoS matrix
                    std::vector<double> x(discr_size_EoS, 0);
                    for (int i = 1; i < discr_size_EoS - 1; ++i)
                    { // cache EoS for further efficiency
                        x[i] = i * (nbar_upp - nbar_low) / discr_size_EoS + nbar_low;
                        cache[0][i] = energy_density_of_nbar(x[i]);
                        cache[1][i] = pressure_of_nbar(x[i]);
                    }
                    x[0] = nbar_low;
                    x[x.size() - 1] = nbar_upp;
                    cache[0][0] = edensity_low;
                    cache[0][cache[0].size() - 1] = edensity_upp;
                    cache[1][0] = pressure_low;
                    cache[1][cache[1].size() - 1] = pressure_upp;
                    eos_interpolator_cached.erase(); // clean up cached interpolator
                }
                return eos_interpolator(cache[0], cache[1], rho);
            });

        // TOV solver

        auto tov_cached = auxiliaries::math::CachedFunc<std::vector<std::vector<double>>, std::vector<double>,
                                                  const std::function<double(double)> &, double, double, double, double>(tov_solver::tov_solution);
        auto tov = [&tov_cached, &eos_cached, edensity](double r)
        {
            // TOV solution cached
            return tov_cached(eos_cached, r, edensity, radius_step, density_step);
        };

        double r_ns = tov(0.0)[4];
        double m_ns = tov(r_ns)[0];
        // clean all global cache
        // None to clean actually

        return std::vector<double>({r_ns, m_ns});
    };

    size_t n = 500;
    size_t offset = 15;
    std::vector<double> x, y;
    // assemble data for different center densities
    for(size_t count = offset; count < n - offset + 1; ++count)
    {
        using namespace constants::conversion;
        double edensity = (count * edensity_upp) / n;
        auto point = get_m_r_at_density(edensity);
        x.push_back(point[0] / km_gev);
        y.push_back(point[1] * gev_over_msol);
        std::cout << "At " << count << " out of " << n << " cycles. M = " << y.back() << " Ms, R = " << x.back() << " km\n";
    }
}