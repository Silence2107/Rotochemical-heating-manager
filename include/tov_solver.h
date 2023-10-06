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
	/// what gets solved by vector RK4 method with initial values m(0)=0, rho(0)=center_density, exp(2phi(R))=(1-2GM/R) until pressure at reached point is less than surface_pressure
	/// @param cache cache support via auxiliaries::math::CachedFunc
	/// @param eos energy density (GeV^4) -> pressure (GeV^4)
	/// @param r distance from centre of NS to a considered point
	/// @param center_density initial value for energy density at r = 0
	/// @param radius_step RK4 base radius step (may be adapted if the process desires so)
	/// @param surface_pressure desired pressure at the surface of NS
	/// @param adaption_limit maximum number of adaption cuts of radius_step. If is exceeded, the algorithm will conclude the simulation immediately
	/// @return [0] -> mass in given point r, [1] -> density -//-, [2] -> phi metric function -//-, [3] -> pressure -//-, [4] -> additionals(radius of NS)
	/// @cite TOV - Oppenheimer, Volkoff, 1939
	std::vector<double> tov_solution(std::vector<std::vector<double>> &cache, const std::function<double(double)> &eos, double r, double center_density, double radius_step, double surface_pressure, size_t adaption_limit);
}

#endif