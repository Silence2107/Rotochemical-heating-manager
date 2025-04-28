
#include "../include/tov_solver.h"
#include "../include/constants.h"
#include "../include/auxiliaries.h"

#include <vector>
#include <limits>
#include <cmath>
#include <sstream>

std::vector<double> tov_solver::tov_solution(std::vector<auxiliaries::math::Interpolator> &cache, const std::function<double(double)> &eos_inv, double r, double center_pressure, double radius_step, double surface_pressure, double lowest_pressure, size_t adaption_limit, auxiliaries::math::Interpolator::InterpolationMode mode)
{
	using constants::scientific::G;
	using constants::scientific::Pi;
	auxiliaries::io::Logger logger(__func__);
	// Here I apply RK4 to coupled equations {m'=f(p,r), p'=g(p,m,r)}
	// For reference, see "Computational Quantum Mechanics" by Joshua Izaac, Jingbo Wang, section 5.9 and https://www.myphysicslab.com/explain/runge-kutta-en.html

	auto m_prime = [&eos_inv, lowest_pressure](double p, double r)
	{
		// p < lowest_pressure is handled by main loop
		double rho = eos_inv(p);
		return 4 * Pi * r * r * rho;
	};

	auto p_prime = [&eos_inv, radius_step, lowest_pressure](double p, double m, double r)
	{
		// p < lowest_pressure is handled by main loop
		double rho = eos_inv(p);
		// zero radius approx -- explicitly get rid of singularity
		if (r < radius_step / 2)
			return -(4 * Pi * G * r) * (rho + p) * (rho / 3 + p);
		// otherwise proceed
		return -G / (r * r) * (rho + p) * (m + 4 * Pi * r * r * r * p) / (1 - 2 * G * m / r);
	};

	auto phi_integrand = [&eos_inv](double p)
	{
		double rho = eos_inv(p);
		return -1.0 / (rho + p);
	}; // integrand of phi function wrt pressure

	double adaptive_radius_step = radius_step; // this one radius step may get smaller if we encounter problems
	size_t adaption_count = 0;				   // if the algorithm fails to converge under adaption_limit, we conclude the simulation
											   // consider moving it to outer scope
	if (cache.empty())
	{
		auto df = std::vector<std::vector<double>>(4, std::vector<double>());
		double m(0), p(center_pressure), r(0), phi(0); // initial values; for phi we know phi(R) = 1/2 ln(1-2GM/R), so we will need to shift total function later
		df[0].push_back(r);
		df[1].push_back(m);
		df[2].push_back(p);
		df[3].push_back(phi);

		while (adaption_count < adaption_limit)
		{
			// while pressure is above desided minima i.e. were not on surface (see below)
			// and we didn't exceed adaption limit
			// then proceed
			size_t rk_size = 4;
			std::vector<double> m_rk(rk_size), p_rk(rk_size); // RK4 variables

			bool adaption_flag = false;
			// RK4 algorithm follows
			auto weights = std::vector<double>({0.5, 0.5, 1});
			for (size_t rk_index = 0; rk_index < rk_size; ++rk_index)
			{
				double p_cur = p,
					   m_cur = m,
					   r_cur = r;
				if (rk_index > 0)
				{
					p_cur += weights[rk_index - 1] * adaptive_radius_step * p_rk[rk_index - 1];
					m_cur += weights[rk_index - 1] * adaptive_radius_step * m_rk[rk_index - 1];
					r_cur += weights[rk_index - 1] * adaptive_radius_step;
				}
				if (p_cur < lowest_pressure)
				{
					adaption_flag = true;
					logger.log([&]()
							   { return true; }, auxiliaries::io::Logger::LogLevel::kTrace,
							   [&]()
							   { 
									std::stringstream ss;
									ss << "Too rapid RK4. RK weight #" + std::to_string(rk_index + 1) + "/4 wandered below lowest pressure. ";
									ss << "r[km] = " << r / constants::conversion::km_gev << ", p[MeV/fm^3] = " << 
									p / constants::conversion::mev_over_fm3_gev4 << ", m[Ms] = " << m * constants::conversion::gev_over_msol;
									return ss.str(); }, "TOV loop");
					break;
				}
				m_rk[rk_index] = m_prime(p_cur, r_cur);
				p_rk[rk_index] = p_prime(p_cur, m_cur, r_cur);
			}

			if (adaption_flag)
			{
				// if we wander outside EoS, reevaluate with smaller step
				adaptive_radius_step /= 2;
				++adaption_count;
				logger.log([&]()
							{ return true; }, auxiliaries::io::Logger::LogLevel::kDebug,
							[&]()
							{ return "Too rapid RK4. Adaption #" + std::to_string(adaption_count) + "/" + std::to_string(adaption_limit) + " at r[km] = " + std::to_string(r / constants::conversion::km_gev); }, "TOV loop");
				continue;
			}
			p += adaptive_radius_step / 6 * (p_rk[0] + 2 * p_rk[1] + 2 * p_rk[2] + p_rk[3]);
			if (p < lowest_pressure)
			{
				// if we wander outside EoS, reevaluate with smaller step
				p -= adaptive_radius_step / 6 * (p_rk[0] + 2 * p_rk[1] + 2 * p_rk[2] + p_rk[3]);
				adaptive_radius_step /= 2;
				++adaption_count;
				logger.log([&]()
						   { return true; }, auxiliaries::io::Logger::LogLevel::kDebug,
						   [&]()
						   { return "Surface overshoot. Adaption #" + std::to_string(adaption_count) + "/" + std::to_string(adaption_limit) + " at r[km] = " + std::to_string(r / constants::conversion::km_gev); }, "TOV loop");
				continue;
			}
			double linear_coeff = 1;
			bool surfaced = false;
			// if we are below surface pressure, we exit EARLY to avoid going beyond EoS domain
			if (p < surface_pressure)
			{
				// find linear interpolation of pressure between last two points
				linear_coeff = (df[2].back() - surface_pressure) / (df[2].back() - p);
				p = surface_pressure;
				surfaced = true;
			}
			m += linear_coeff * adaptive_radius_step / 6 * (m_rk[0] + 2 * m_rk[1] + 2 * m_rk[2] + m_rk[3]);
			r += linear_coeff * adaptive_radius_step;
			phi += 0.5 * linear_coeff * (p - df[2].back()) * (phi_integrand(p) + phi_integrand(df[2].back())); // trapezoid integral
			// memorize results
			df[0].push_back(r);
			df[1].push_back(m);
			df[2].push_back(p);
			df[3].push_back(phi);

			adaptive_radius_step = radius_step; // reset adaptive radius step
			adaption_count = 0;					// reset adaption count
			// exit if reached surface
			if (surfaced)
				break;
		}

		logger.log([&]()
				   { return adaption_count >= adaption_limit; }, auxiliaries::io::Logger::LogLevel::kDebug,
				   [&]()
				   { return "TOV exits early (adaption limit exceeded). Perhaps raise TOVSolver.AdaptionLimit"; });

		// shift phi function so that to satisfy phi(R) condition
		double phi_shift = 1.0 / 2 * log(1 - 2 * G * df[1].back() / df[0].back()) - df[3].back();
		for (auto &elem : df[3])
			elem += phi_shift;

		cache = std::vector<auxiliaries::math::Interpolator>(4, auxiliaries::math::Interpolator());

		for (size_t i = 0; i < 3; ++i)
		{
			cache[i].instantiate(df[0], df[i + 1], mode);
		}

		// make it work for now, overhaul soon..
		cache[3].instantiate({df[0].front(), df[0].back()}, {df[0].back(), df[0].back()}, mode);
	}

	// Interpolation

	double mass = cache[0](r);
	double pressure = cache[1](r);
	double phi = cache[2](r);
	double r_ns = cache[3](0.0);
	return std::vector<double>({mass, eos_inv(pressure), phi, pressure, r_ns});
}