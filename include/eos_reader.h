#ifndef EOS_READER_H
#define EOS_READER_H

#include <fstream>
#include <vector>
#include <functional>

/// @brief contains:
/// enum EOS with values EOS::APR4;
/// function eos_data;
/// function apr4
/// function ist_for_ns
namespace eos_reader
{
	/// @brief predefined EoS calculators and corresponding enum values
	namespace predefined
	{ 
		/// @brief used to distinguish EoS when passing datafile to eos_data
		/// @deprecated use explicit EoS functions instead
		enum class EOS
		{
			APR4,
			IST
		};

		/// @brief APR4 EoS reader; Gets invoked when EOS::APR4 is passed to eos_data
		/// @param input input[0] - barionic density in fm-3
		/// @param fstr File to read from
		/// @return (barionic density &gt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, electron fraction, muon -//-, neutron -//-, proton -//-)
		/// 		(barionic density &lt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, Acell, Aion, Z)
		std::vector<double> apr4(const std::vector<double> &input, std::ifstream &fstr);

		/// @brief APR4 EoS reader; Gets invoked when EOS::APR4 is passed to eos_data
		/// @param cache Cache support. Wrap this function with auxiliaries::CachedFunc for usage
		/// @param input input[0] - barionic density in fm-3
		/// @param fstr File to read from
		/// @return (barionic density &gt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, electron fraction, muon -//-, neutron -//-, proton -//-)
		/// 		(barionic density &lt; 0.055 fm-3) -> (energy density g/cm3, pressure dyne/cm2, barionic density fm-3, Acell, Aion, Z)
		std::vector<double> apr4(std::vector<std::vector<double>> &cache, const std::vector<double> &input, std::ifstream &fstr);

		/// @brief IST for NS EoS reader; Gets invoked when EOS::IST is passed to eos_data
		/// @param input input[0] - barionic density in fm-3
		/// @param fstr File to read from
		/// @return (barionic density &lt; 0.0999 fm-3) -> (pressure MeV/fm3, energy density MeV/fm3, barionic chem. potential MeV, barionic density fm-3)
		/// 		(barionic density &gt; 0.0999 fm-3) -> (pressure MeV/fm3, energy density MeV/fm3, barionic chem. potential MeV, barionic density fm-3, electric chem. potential MeV, electron density fm-3, neutron density fm-3, proton density fm-3)
		std::vector<double> ist_for_ns(const std::vector<double> &input, std::ifstream &fstr);
	}

	/// @brief General purpose EoS reader; bases on known EoS filestyle and input stream with such file allows to extract various data
	/// @param input Input that is required for EoS
	/// @param eos EoS from the eos_reader::predefined::EOS
	/// @param fstr EoS datafile
	/// @return EoS output
	/// @deprecated use explicit EoS functions instead
	std::vector<double> eos_data(const std::vector<double> &input, eos_reader::predefined::EOS eos, std::ifstream &fstr);

	/// @brief General purpose EoS reader; bases on known EoS filestyle and input stream with such file allows to extract various data
	/// @param input Input that is required for EoS
	/// @param eos An EoS function that takes input and fstr and returns EoS output
	/// @param fstr EoS datafile
	/// @return EoS output
	std::vector<double> eos_data(const std::vector<double> &input, const std::function<std::vector<double>(const std::vector<double> &, std::ifstream &)> &eos, std::ifstream &fstr);
}

#endif