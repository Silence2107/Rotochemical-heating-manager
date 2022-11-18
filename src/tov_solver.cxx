
#include "../include/tov_solver.h"
#include "../include/constants.h"
#include "../include/auxiliaries.h"

#include <vector>
#include <iostream>
#include <cmath>

std::vector<double> tov_solver::tov_solution(std::vector<std::vector<double>> &cache, const std::function<double(double)> &eos, double r, double center_density, double radius_step, double density_step)
{

	using constants::scientific::G;
	using constants::scientific::Pi;

	// for P'_rho I will use the following function
	auto eos_prime = [&eos, density_step](double dens)
	{
		try
		{
			return (eos(dens + density_step / 2) - eos(dens - density_step / 2)) / (density_step);
		}
		catch (const std::exception &e)
		{
			throw std::runtime_error("Pressure derivative computation failed; Encountered in tov_solver::tov_solution; " + std::string(e.what()));
		}
		// maybe should implement leap-frog as well/instead
	};

	// Here I apply RK4 to coupled equations {m'=f(rho,r), rho'=g(rho,m,r)}
	// For reference, see "Computational Quantum Mechanics" by Joshua Izaac, Jingbo Wang, section 5.9 and https://www.myphysicslab.com/explain/runge-kutta-en.html

	auto m_prime = [](double rho, double r)
	{
		return 4 * Pi * r * r * rho;
	};

	auto rho_prime = [&eos, eos_prime, radius_step](double rho, double m, double r)
	{
		if (r < radius_step / 2)
			return 0.0;
		return -G / (r * r * eos_prime(rho)) * (rho + eos(rho)) * (m + 4 * Pi * r * r * r * eos(rho)) / (1 - 2 * G * m / r);
	};

	// check whether we need to recache
	bool emptyflag = false;
	if (cache.empty())
	{
		cache = std::vector<std::vector<double>>(5, std::vector<double>());
		cache[4].push_back(0); // make sure cache[4][0] exists so that further check does not fail
		emptyflag = true;
	}
	if (emptyflag || fabs(cache[4][0] - center_density) > density_step)
	{ // recache in case cache is empty or center_density differs from cached
		cache = std::vector<std::vector<double>>(5, std::vector<double>());
		double m(0), rho(center_density), r(0), phi(0); // initial values; for phi we know phi(R) = 1/2 ln(1-2GM/R), so we will need to shift total function later
		cache[0].push_back(r);
		cache[1].push_back(m);
		cache[2].push_back(rho);
		cache[3].push_back(phi);
		cache[4].push_back(center_density);

		auto phi_integrand = [eos_prime, &eos](double rho)
		{
			return -eos_prime(rho) / (rho + eos(rho));
		}; // integrand to find phi(r)

		while (rho > density_step)
		{ // while density is not close to zero i.e. were not on surface
			// then proceed
			std::vector<double> m_rk(4), rho_rk(4); // RK4 variables

			// std::cout << r << " : " << m << " " << rho << '\n';

			// RK4 algorithm follows
			m_rk[0] = m_prime(rho, r);
			rho_rk[0] = rho_prime(rho, m, r);

			if (rho + radius_step / 2 * rho_rk[0] < 0)
				break;
			m_rk[1] = m_prime(rho + radius_step / 2 * rho_rk[0], r + radius_step / 2);
			rho_rk[1] = rho_prime(rho + radius_step / 2 * rho_rk[0], m + radius_step / 2 * m_rk[0], r + radius_step / 2);

			if (rho + radius_step / 2 * rho_rk[1] < 0)
				break;
			m_rk[2] = m_prime(rho + radius_step / 2 * rho_rk[1], r + radius_step / 2);
			rho_rk[2] = rho_prime(rho + radius_step / 2 * rho_rk[1], m + radius_step / 2 * m_rk[1], r + radius_step / 2);

			if (rho + radius_step * rho_rk[2] < 0)
				break;
			m_rk[3] = m_prime(rho + radius_step * rho_rk[2], r + radius_step);
			rho_rk[3] = rho_prime(rho + radius_step * rho_rk[2], m + radius_step * m_rk[2], r + radius_step);

			rho += radius_step / 6 * (rho_rk[0] + 2 * rho_rk[1] + 2 * rho_rk[2] + rho_rk[3]);
			m += radius_step / 6 * (m_rk[0] + 2 * m_rk[1] + 2 * m_rk[2] + m_rk[3]);
			r += radius_step;
			phi += 0.5 * (rho - cache[2].back()) * (phi_integrand(rho) + phi_integrand(cache[2].back())); // trapezoid integral
			// memorize results
			cache[0].push_back(r);
			cache[1].push_back(m);
			cache[2].push_back(rho);
			cache[3].push_back(phi);
		}

		// shift phi function so that to satisfy phi(R) condition
		double phi_shift = 1.0 / 2 * log(1 - 2 * G * cache[1].back() / cache[0].back()) - cache[3].back();
		for (auto &elem : cache[3])
			elem += phi_shift;
	}

	/*if (cache[0].back() <= r)
	{ // output in case provided coordinate is out of range
		return std::vector<double>({cache[1].back(), cache[2].back(), cache[3].back(), eos(cache[2].back()), cache[0].back()});
	}*/

	// linear interpolation

	double mass = auxiliaries::interpolate(cache[0], cache[1], auxiliaries::InterpolationMode::kLinear, r);
	double density = auxiliaries::interpolate(cache[0], cache[2], auxiliaries::InterpolationMode::kLinear, r);
	double phi = auxiliaries::interpolate(cache[0], cache[3], auxiliaries::InterpolationMode::kLinear, r);
	return std::vector<double>({mass, density, phi, eos(density), cache[0].back()});
}