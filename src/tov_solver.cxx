
#include "../include/tov_solver.h"
#include "../include/constants.h"
#include "../include/auxiliaries.h"

#include <vector>
#include <limits>
#include <cmath>

std::vector<double> tov_solver::tov_solution(std::vector<std::vector<double>> &cache, const std::function<double(double)> &eos, double r, double center_density, double radius_step, double surface_pressure, size_t adaption_limit)
{
	using constants::scientific::G;
	using constants::scientific::Pi;

	// for P'_rho I will use the following function
	auto eos_prime = [&eos, center_density](double dens)
	{
		double epsilon = std::numeric_limits<double>::epsilon();
		// hope no one puts negative density in here
		double dens_step = std::sqrt(epsilon) * (dens + std::sqrt(epsilon));
		// make sure to not go beyond center density requests, since we may easily get outside of EOS domain
		if (dens + dens_step > center_density)
			return (eos(dens) - eos(dens - dens_step)) / dens_step;
		// otherwise we can safely go forward
		return (eos(dens + dens_step) - eos(dens)) / dens_step;
	};

	// Here I apply RK4 to coupled equations {m'=f(rho,r), rho'=g(rho,m,r)}
	// For reference, see "Computational Quantum Mechanics" by Joshua Izaac, Jingbo Wang, section 5.9 and https://www.myphysicslab.com/explain/runge-kutta-en.html

	auto m_prime = [](double rho, double r)
	{
		if (rho < 0)
			RHM_THROW(std::runtime_error, "Negative density encountered.");
		return 4 * Pi * r * r * rho;
	};

	auto rho_prime = [&eos, &eos_prime, radius_step](double rho, double m, double r)
	{
		if (rho < 0)
			RHM_THROW(std::runtime_error, "Negative density encountered.");
		// zero radius approx -- explicitly get rid of singularity
		if (r < radius_step / 2)
			return -(4 * Pi * G * r) / eos_prime(rho) * (rho + eos(rho)) * (rho / 3 + eos(rho));
		// otherwise proceed
		return -G / (r * r * eos_prime(rho)) * (rho + eos(rho)) * (m + 4 * Pi * r * r * r * eos(rho)) / (1 - 2 * G * m / r);
	};

	auto phi_integrand = [&eos_prime, &eos](double rho)
	{
		return -eos_prime(rho) / (rho + eos(rho));
	}; // integrand of phi function wrt rho

	double adaptive_radius_step = radius_step; // this one radius step may get smaller if we encounter problems
	size_t adaption_count = 0; // if the algorithm fails to converge under adaption_limit, we conclude the simulation
							   // consider moving it to outer scope
	if (cache.empty())
	{
		cache = std::vector<std::vector<double>>(4, std::vector<double>());
		double m(0), rho(center_density), r(0), phi(0); // initial values; for phi we know phi(R) = 1/2 ln(1-2GM/R), so we will need to shift total function later
		cache[0].push_back(r);
		cache[1].push_back(m);
		cache[2].push_back(rho);
		cache[3].push_back(phi);

		while (adaption_count < adaption_limit)
		{ 
			// while pressure is above desided minima i.e. were not on surface (see below)
			// and we didn't exceed adaption limit
			// then proceed
			std::vector<double> m_rk(4), rho_rk(4); // RK4 variables

			// RK4 algorithm follows
			try
			{
				m_rk[0] = m_prime(rho, r);
				rho_rk[0] = rho_prime(rho, m, r);

				m_rk[1] = m_prime(rho + adaptive_radius_step / 2 * rho_rk[0], r + adaptive_radius_step / 2);
				rho_rk[1] = rho_prime(rho + adaptive_radius_step / 2 * rho_rk[0], m + adaptive_radius_step / 2 * m_rk[0], r + adaptive_radius_step / 2);

				m_rk[2] = m_prime(rho + adaptive_radius_step / 2 * rho_rk[1], r + adaptive_radius_step / 2);
				rho_rk[2] = rho_prime(rho + adaptive_radius_step / 2 * rho_rk[1], m + adaptive_radius_step / 2 * m_rk[1], r + adaptive_radius_step / 2);

				m_rk[3] = m_prime(rho + adaptive_radius_step * rho_rk[2], r + adaptive_radius_step);
				rho_rk[3] = rho_prime(rho + adaptive_radius_step * rho_rk[2], m + adaptive_radius_step * m_rk[2], r + adaptive_radius_step);
			}
			catch (std::runtime_error &e)
			{
				// if we encounter negative density, reevaluate with smaller step
				adaptive_radius_step /= 2;
				++adaption_count;
				continue;
			}

			rho += adaptive_radius_step / 6 * (rho_rk[0] + 2 * rho_rk[1] + 2 * rho_rk[2] + rho_rk[3]);
			if (rho < 0.0)
			{
				// if we encounter negative density, reevaluate with smaller step
				rho -= adaptive_radius_step / 6 * (rho_rk[0] + 2 * rho_rk[1] + 2 * rho_rk[2] + rho_rk[3]);
				adaptive_radius_step /= 2;
				++adaption_count;
				continue;
			}
			// if we are below surface pressure, we exit EARLY to avoid going beyond EoS domain
			if (eos(rho) < surface_pressure)
				break; 
			m += adaptive_radius_step / 6 * (m_rk[0] + 2 * m_rk[1] + 2 * m_rk[2] + m_rk[3]);
			r += adaptive_radius_step;
			phi += 0.5 * (rho - cache[2].back()) * (phi_integrand(rho) + phi_integrand(cache[2].back())); // trapezoid integral
			// memorize results
			cache[0].push_back(r);
			cache[1].push_back(m);
			cache[2].push_back(rho);
			cache[3].push_back(phi);
			adaptive_radius_step = radius_step; // reset adaptive radius step
			adaption_count = 0;					// reset adaption count
		}

		// shift phi function so that to satisfy phi(R) condition
		double phi_shift = 1.0 / 2 * log(1 - 2 * G * cache[1].back() / cache[0].back()) - cache[3].back();
		for (auto &elem : cache[3])
			elem += phi_shift;
	}

	// linear interpolation

	double mass = auxiliaries::math::interpolate(cache[0], cache[1], auxiliaries::math::InterpolationMode::kLinear, r);
	double density = auxiliaries::math::interpolate(cache[0], cache[2], auxiliaries::math::InterpolationMode::kLinear, r);
	double phi = auxiliaries::math::interpolate(cache[0], cache[3], auxiliaries::math::InterpolationMode::kLinear, r);
	return std::vector<double>({mass, density, phi, eos(density), cache[0].back()});
}