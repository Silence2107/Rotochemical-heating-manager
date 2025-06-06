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
	/// @param eos_inv pressure (GeV^4) -> energy density (GeV^4)
	/// @param center_pressure initial value for pressure at r = 0
	/// @param radius_step RK4 base radius step (may be adapted if the process desires so)
	/// @param surface_pressure desired pressure at the surface of NS
	/// @param lowest_pressure lowest pressure which EoS can provide
	/// @param adaption_limit maximum number of adaption cuts of radius_step. If is exceeded, the algorithm will conclude the simulation immediately
	/// @return vectors [0] -> radii, [1] -> mass in given point r, [2] -> pressure -//-, [3] -> phi metric function -//-
	/// @cite TOV - Oppenheimer, Volkoff, 1939
	std::vector<std::vector<double>> tov_solution(const std::function<double(double)> &eos_inv, double center_pressure, double radius_step, double surface_pressure, double lowest_pressure, size_t adaption_limit);

	/// @brief Solves for dimensionless auxiliary function for tidal deformability based on TOV solution and distance from centre of NS r (GeV^-1);
	/// dy(r)/dr = -y(r)^2/r -y(r)F(r)/r - Q(r)/r, where F(r) = [1 - 4piGr^2(rho(r) - p(r))] * [1 - 2Gm(r)/r]^-1;
	/// Q(r) = 4piGr^2[5rho(r) + 9p(r) + (rho(r) + p(r))*drho/dp - 6/(4piGr^2)] * [1 - 2Gm(r)/r]^-1
	///     − (4G^2m^2(r) / r^2) * [1 + (4pir^3 p(r)) / m(r)]^2 * [1 - 2Gm(r) / r]^-2
	/// @param edensity_of_r energy density (GeV^4) at given r (GeV^-1)
	/// @param pressure_of_r pressure (GeV^4) at given r (GeV^-1)
	/// @param mass_of_r enclosed mass (GeV) at given r (GeV^-1)
	/// @param radius_step RK4 radius step
	/// @param r_ns radius of the neutron star (GeV^-1)
	/// @cite Leung, Chu et al, 2022
	std::vector<std::vector<double>> tidal_solution(const std::function<double(double)> &edensity_of_r, const std::function<double(double)> &pressure_of_r, const std::function<double(double)> &mass_of_r, double radius_step, double r_ns);
}

#endif