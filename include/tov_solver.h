#ifndef TOV_SOLVER_H
#define TOV_SOLVER_H

#include "../include/auxiliaries.h"

#include <functional>
#include <vector>

/// @brief contains: TOV solvers
namespace tov_solver
{
	/// @brief Solves TOV system of equations based on eos and distance from centre of NS r (GeV^-1);
	/// TOV system reduces to m'(r) = 4pir^2 rho(P(r)); P'(r) = -G/r^2 * (rho(P(r))+P(r)) * (m(r)+4pir^3 P(r))/(1-2Gm(r)/r),
	/// Phi'_p(r) = -1/(rho(P(r))+P(r)), rho(r) = rho(P(r)),
	/// what gets solved by vector RK4 method with initial values m(0)=0, P(0)=center_pressure, exp(2phi(R))=(1-2GM/R) until pressure at reached point is less than surface_pressure
	/// @param cache cache support via auxiliaries::math::CachedFunc
	/// @param eos_inv pressure (GeV^4) -> energy density (GeV^4)
	/// @param r distance from centre of NS to a considered point
	/// @param center_pressure initial value for pressure at r = 0
	/// @param radius_step RK4 base radius step (may be adapted if the process desires so)
	/// @param surface_pressure desired pressure at the surface of NS
	/// @param adaption_limit maximum number of adaption cuts of radius_step. If is exceeded, the algorithm will conclude the simulation immediately
	/// @param mode interpolation mode for radial functions
	/// @return [0] -> mass in given point r, [1] -> density -//-, [2] -> phi metric function -//-, [3] -> pressure -//-, [4] -> additionals(radius of NS)
	/// @cite TOV - Oppenheimer, Volkoff, 1939
	std::vector<double> tov_solution(std::vector<std::function<double(double)>> &cache, const std::function<double(double)> &eos_inv, double r, double center_pressure, double radius_step, double surface_pressure, size_t adaption_limit, auxiliaries::math::InterpolationMode mode);
}

#endif