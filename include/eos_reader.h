#ifndef EOS_READER_H
#define EOS_READER_H

#include "../include/auxiliaries.h"

#include <fstream>
#include <vector>
#include <functional>

/// @brief Various functionality for reading EoS datafiles
namespace eos_reader
{
	/// @brief predefined EoS readers
	namespace predefined
	{
		/// @brief APR4 EoS reader
		/// @param input input[0] - barionic density in fm-3
		/// @param fstr File to read from
		/// @return (barionic density &gt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, electron fraction, muon -//-, neutron -//-, proton -//-)
		///  		(barionic density &lt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, Acell, Aion, Z)
		/// @deprecated limits usage of interpolation to kLinear; can also be slow. Use cached version instead
		std::vector<double> apr4(const std::vector<double> &input, std::ifstream &fstr);

		/// @brief APR4 EoS reader
		/// @param cache Cache support. Wrap this function with auxiliaries::CachedFunc for usage
		/// @param input input[0] - barionic density in fm-3
		/// @param fstr File to read from
		/// @param interpolator Interpolation function of type double(const std::vector<double> &input, const std::vector<double> &output, double val). Defaults to linear interpolation.
		/// Cached interpolator is highly discouraged(results in wrong values)
		/// @return (barionic density &gt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, electron fraction, muon -//-, neutron -//-, proton -//-)
		///			(barionic density &lt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, Acell, Aion, Z)
		std::vector<double> apr4_cached(
			std::vector<std::vector<double>> &cache, const std::vector<double> &input, std::ifstream &fstr,
			const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator = [](const std::vector<double> &x, const std::vector<double> &y, double val)
			{ return auxiliaries::interpolate(x, y, auxiliaries::InterpolationMode::kLinear, val); });

		/// @brief IST for NS EoS reader
		/// @param input input[0] - barionic density in fm-3
		/// @param fstr File to read from
		/// @return (barionic density &lt; 0.0999 fm-3) -> (pressure MeV/fm3, energy density MeV/fm3, barionic chem. potential MeV, barionic density fm-3)
		/// 		(barionic density &gt; 0.0999 fm-3) -> (pressure MeV/fm3, energy density MeV/fm3, barionic chem. potential MeV, barionic density fm-3, electric chem. potential MeV, electron density fm-3, neutron density fm-3, proton density fm-3)
		/// @deprecated limits usage of interpolation to kLinear; can also be slow. Use cached version instead
		std::vector<double> ist_for_ns(const std::vector<double> &input, std::ifstream &fstr);

		/// @brief IST for NS EoS reader
		/// @param cache Cache support. Wrap this function with auxiliaries::CachedFunc for usage
		/// @param input input[0] - barionic density in fm-3
		/// @param fstr File to read from
		/// @param interpolator Interpolation function of type double(const std::vector<double> &input, const std::vector<double> &output, double val). Defaults to linear interpolation.
		/// Cached interpolator is highly discouraged(results in wrong values)
		/// @return (barionic density &lt; 0.0999 fm-3) -> (pressure MeV/fm3, energy density MeV/fm3, barionic chem. potential MeV, barionic density fm-3)
		/// 		(barionic density &gt; 0.0999 fm-3) -> (pressure MeV/fm3, energy density MeV/fm3, barionic chem. potential MeV, barionic density fm-3, electric chem. potential MeV, electron density fm-3, neutron density fm-3, proton density fm-3)
		std::vector<double> ist_for_ns_cached(
			std::vector<std::vector<double>> &cache, const std::vector<double> &input, std::ifstream &fstr,
			const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator = [](const std::vector<double> &x, const std::vector<double> &y, double val)
			{ return auxiliaries::interpolate(x, y, auxiliaries::InterpolationMode::kLinear, val); });
	}

	/// @brief General purpose EoS reader; bases on known EoS filestyle and input stream with such file allows to extract various data
	/// @param input Input that is required for EoS
	/// @param eos An EoS function that takes input and fstr and returns EoS output
	/// @param fstr EoS datafile
	/// @return EoS output
	std::vector<double> eos_data(const std::vector<double> &input, const std::function<std::vector<double>(const std::vector<double> &, std::ifstream &)> &eos, std::ifstream &fstr);
}

#endif