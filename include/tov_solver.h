#ifndef TOV_SOLVER_H
#define TOV_SOLVER_H

#include <functional>
#include <vector>

/// @brief contains: function tov_solution
namespace tov_solver
{
	/// @brief Solves TOV system of equations based on eos and distance from centre of NS r (GeV^-1);
	/// TOV system reduces to m'(r) = 4pir^2 rho(r); rho'(r) = -G/(r^2 * P'_rho(rho(r)) ) * (rho(r)+P(rho(r))) * (m(r)+4pir^3 P(rho(r)))/(1-2Gm(r)/r),
	/// Phi'(r) = P'(r)/(rho(r)+P(rho(r))), P(r) = P(rho(r)),
	/// what gets solved by vector RK4 method with initial values m(0)=0, rho(0)=center_density, exp(2phi(R))=(1-2GM/R) until such r that rho(r) &lt; density_step
	/// ; !!WRAP IT WITH auxiliaries::CachedFunc!!
	/// @param cache structure: cache[0] - tabulated radius, cache[1] - tabulated mass, cache[2] - tabulated density, cache[3] - tabulated phi, cache[4][0] - center_density
	/// ; in case center_density differs from cached value for more than density_step, cache gets reevaluated
	/// ; in case you need to change steps for eos or radius, or eos itself, make sure to clean cache by yourself (call CachedFunc::erase)
	/// @param eos energy density (GeV^4) -> pressure (GeV^4)
	/// @param r distance from centre of NS to a considered point
	/// @param center_density initial value for energy density at r = 0
	/// @param radius_step is used for RK4 algorithm (typical scale is way less than 1E19)
	/// @param density_step (typical scale is way less than 1E-4 - 1E-2) is used for eos derivative, for caching purposes (if |cached_density-center_density|&lt;density_step then no need for recaching) and for finding the radius of NS
	/// @return [0] -> mass in given point r, [1] -> density -//-, [2] -> phi metric function -//-, [3] -> pressure -//-, [4] -> additionals(radius of NS)
	std::vector<double> tov_solution(std::vector<std::vector<double>> &cache, const std::function<double(double)> &eos, double r, double center_density, double radius_step, double density_step);
}

#endif