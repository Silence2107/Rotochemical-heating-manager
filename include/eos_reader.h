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
		/// @param cache Cache support. Wrap this function with auxiliaries::CachedFunc for usage
		/// @param input input[0] - barionic density in fm-3
		/// @param fstr File to read from
		/// @param interpolator Interpolation function of type double(const std::vector<double> &input, const std::vector<double> &output, double val). 
		/// Cached interpolator must not be used, as well as one should limit themselves to linear interpolation, since same column may carry different info
		/// @return (barionic density &gt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, electron fraction, muon -//-, neutron -//-, proton -//-, lambda -//-, sigma- -//-, sigma0 -//-, sigma+ -//-, m star proton -//-, m star neutron -//-, m star lambda -//-, m star sigma- -//-, m star sigma0 -//-, m star sigma+ -//-)
		///			(barionic density &lt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, Acell, Aion, Z, [empty])
		std::vector<double> apr4_cached(
			std::vector<std::vector<double>> &cache, const std::vector<double> &input, std::ifstream &fstr,
			const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator);

		/* DEPRECATED VERSION FOR OLD FILE FORMAT
		/// @brief IST for NS EoS reader
		/// @param cache Cache support. Wrap this function with auxiliaries::CachedFunc for usage
		/// @param input input[0] - barionic density in fm-3
		/// @param fstr File to read from
		/// @param interpolator Interpolation function of type double(const std::vector<double> &input, const std::vector<double> &output, double val).
		/// Cached interpolator must not be used, as well as one should limit themselves to linear interpolation, since same column may carry different info
		/// @return (barionic density &lt; 0.0999 fm-3) -> (pressure MeV/fm3, energy density MeV/fm3, barionic chem. potential MeV, barionic density fm-3, [empty])
		/// 		(barionic density &gt; 0.0999 fm-3) -> (pressure MeV/fm3, energy density MeV/fm3, barionic chem. potential MeV, barionic density fm-3, electric chem. potential MeV, electron density fm-3, neutron density fm-3, proton density fm-3)
		std::vector<double> ist_for_ns_cached(
			std::vector<std::vector<double>> &cache, const std::vector<double> &input, std::ifstream &fstr,
			const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator);
		*/
		/// @brief IST for NS EoS reader
		/// @param cache Cache support. Wrap this function with auxiliaries::CachedFunc for usage
		/// @param input input[0] - barionic density in fm-3
		/// @param fstr File to read from
		/// @param interpolator Interpolation function of type double(const std::vector<double> &input, const std::vector<double> &output, double val).
		/// Cached interpolator must not be used, as well as one should limit themselves to linear interpolation, since same column may carry different info
		/// @return (barionic density &gt; 0.094 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, electron fraction, muon -//-, neutron -//-, proton -//-, lambda -//-, sigma- -//-, sigma0 -//-, sigma+ -//-, m star proton -//-, m star neutron -//-, m star lambda -//-, m star sigma- -//-, m star sigma0 -//-, m star sigma+ -//-)
		///			(barionic density &lt; 0.0789 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, Acell, Aion, Z, [empty])
		std::vector<double> ist_for_ns_cached(
			std::vector<std::vector<double>> &cache, const std::vector<double> &input, std::ifstream &fstr,
			const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator);
	}

	/// @brief General purpose EoS reader; bases on known EoS filestyle and input stream with such file allows to extract various data
	/// @param input Input that is required for EoS
	/// @param eos An EoS function that takes input and fstr and returns EoS output
	/// @param fstr EoS datafile
	/// @return EoS output
	std::vector<double> eos_data(const std::vector<double> &input, const std::function<std::vector<double>(const std::vector<double> &, std::ifstream &)> &eos, std::ifstream &fstr);
}

#endif