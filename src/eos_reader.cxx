
#include "../include/eos_reader.h"
#include "../include/constants.h"
#include "../include/auxiliaries.h"

#include <exception>
#include <string>
#include <sstream>

void eos_reader::apr4(const std::vector<double> &input, std::vector<double> &output, std::ifstream &fstr)
{
	using namespace constants::apr4;
	double nbar = input[0];					// barionic density (input[0])
	if (nbar > nbar_upp || nbar < nbar_low) // we do not have data beyond these values
		throw std::runtime_error("Data request out of range; Encountered in eos_reader::apr4");
	else if (nbar > nbar_core_limit)
	{
		std::string nextline, prevline;
		for (int i = 0; std::getline(fstr, nextline) && i < 6; ++i)
			; // skip 6 lines
		while (std::getline(fstr >> std::ws, nextline))
		{ // read line by line, boundary whitespace free
			nextline = auxiliaries::retrieve_cleared_line(nextline); // clear line
			std::stringstream strstr(nextline);
			std::string word;
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' '); // skip two first values so that to compare input[0] with barionic density
			if (std::stod(word) < nbar)		 // while our input is less that value in the line, skip
				break;
			prevline = nextline; // save two lines whose barionic density is higher and lower than input[0]
		}
		std::vector<double> prevdata, nextdata;
		std::stringstream nextstrstr(nextline), prevstrstr(prevline); // read those two lines thorough
		std::string nextword, prevword;
		for (int i = 0; std::getline(nextstrstr, nextword, ' '); ++i)
		{
			// std::getline(nextstrstr, nextword, ' ');
			std::getline(prevstrstr, prevword, ' ');
			if (i < 7)
			{ // extract data on positions you are interested in
				prevdata.push_back(std::stod(prevword));
				nextdata.push_back(std::stod(nextword));
			}
		}
		output.reserve(prevdata.size());
		for (int i = 0; i < prevdata.size(); ++i)
			output[i] = auxiliaries::interpolate({prevdata[2], nextdata[2]}, {prevdata[i], nextdata[i]}
			, auxiliaries::InterpolationMode::kLinear, nbar, true);
	}
	else if (nbar < nbar_crust_limit)
	{
		std::string nextline, prevline; // same as for nbar > nbar_core_limit, but different data extracted
		for (int i = 0; std::getline(fstr, nextline) && i < 185; ++i)
			;
		while (std::getline(fstr >> std::ws, nextline))
		{
			nextline = auxiliaries::retrieve_cleared_line(nextline); // clear line
			std::stringstream strstr(nextline);
			std::string word;
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' ');
			if (std::stod(word) < nbar)
				break;
			prevline = nextline;
		}
		std::vector<double> prevdata, nextdata;
		std::stringstream nextstrstr(nextline), prevstrstr(prevline);
		std::string nextword, prevword;
		for (int i = 0; std::getline(nextstrstr, nextword, ' '); ++i)
		{
			// std::getline(nextstrstr, nextword, ' ');
			std::getline(prevstrstr, prevword, ' ');
			if (i < 6)
			{
				prevdata.push_back(std::stod(prevword));
				nextdata.push_back(std::stod(nextword));
			}
		}
		output.reserve(prevdata.size());
		for (int i = 0; i < prevdata.size(); ++i)
			output[i] = auxiliaries::interpolate({prevdata[2], nextdata[2]}, {prevdata[i], nextdata[i]}, 
			auxiliaries::InterpolationMode::kLinear, nbar, true);
	}
	else
	{
		// extract data between i.e. at phase transition;
		// for densities, pressure and baryonic density, we introduce slight slope for monotony; other entries get copypasted depending on nbar
		if (nbar > (nbar_core_limit + nbar_crust_limit) / 2.0)
		{
			output = std::vector<double>({1.5197E+14, 9.2819E+32, 9.0000E-02, 3.1606E-02, 0.0000E+00, 9.6839E-01, 3.1606E-02});
			// I choose slopes by hand : split 10%/80%/10%
			output[0] -= 2.0 * (nbar_core_limit - nbar) / (nbar_core_limit - nbar_crust_limit) * (1.5197E+14 - 3.493E+13) / 10.0;
			output[1] -= 2.0 * (nbar_core_limit - nbar) / (nbar_core_limit - nbar_crust_limit) * (9.2819E+32 - 7.311E+31) / 10.0;
			output[2] -= 2.0 * (nbar_core_limit - nbar) / (nbar_core_limit - nbar_crust_limit) * (9.0000E-02 - 2.096E-02) / 10.0;
		}
		else
		{
			output = std::vector<double>({3.493E+13, 7.311E+31, 2.096E-02, 1127., 124., 26.});
			// I choose slopes by hand : split 10%/80%/10%
			output[0] += 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (1.5197E+14 - 3.493E+13) / 10.0;
			output[1] += 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (9.2819E+32 - 7.311E+31) / 10.0;
			output[2] += 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (9.0000E-02 - 2.096E-02) / 10.0;
		}
	}
}

void eos_reader::ist_for_ns(const std::vector<double> &input, std::vector<double> &output, std::ifstream &fstr)
{
	using namespace constants::ist_ns;
	double nbar = input[0];					// barionic density (input[0])
	if (nbar > nbar_upp || nbar < nbar_low) // we do not have data beyond these values
		throw std::runtime_error("Data request out of range; Encountered in eos_reader::ist_for_ns");
	else if (nbar > nbar_core_limit)
	{
		std::string nextline, prevline;
		while (std::getline(fstr >> std::ws, nextline))
		{ // read line by line, boundary whitespace free
			nextline = auxiliaries::retrieve_cleared_line(nextline); // clear line
			std::stringstream strstr(nextline);
			std::string word;
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' '); // skip four first values so that to compare input[0] with barionic density
			if (std::stod(word) > nbar)		 // while our input is less that value in the line, skip
				break;
			prevline = nextline; // save two lines whose barionic density is higher and lower than input[0]
		}
		std::vector<double> prevdata, nextdata;
		std::stringstream nextstrstr(nextline), prevstrstr(prevline); // read those two lines thorough
		std::string nextword, prevword;
		for (int i = 0; std::getline(nextstrstr, nextword, ' '); ++i)
		{
			// std::getline(nextstrstr, nextword, ' ');
			std::getline(prevstrstr, prevword, ' ');
			if (i < 8)
			{ // extract data on positions you are interested in
				prevdata.push_back(std::stod(prevword));
				nextdata.push_back(std::stod(nextword));
			}
		}
		output.reserve(prevdata.size());
		for (int i = 0; i < prevdata.size(); ++i)
			output[i] = auxiliaries::interpolate({prevdata[3], nextdata[3]}, {prevdata[i], nextdata[i]}, 
			auxiliaries::InterpolationMode::kLinear, nbar, true);
	}
	else if (nbar < nbar_core_limit)
	{
		std::string nextline, prevline;
		while (std::getline(fstr >> std::ws, nextline))
		{ // read line by line, boundary whitespace free
			nextline = auxiliaries::retrieve_cleared_line(nextline); // clear line
			std::stringstream strstr(nextline);
			std::string word;
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' '); // skip four first values so that to compare input[0] with barionic density
			if (std::stod(word) > nbar)		 // while our input is less that value in the line, skip
				break;
			prevline = nextline; // save two lines whose barionic density is higher and lower than input[0]
		}
		std::vector<double> prevdata, nextdata;
		std::stringstream nextstrstr(nextline), prevstrstr(prevline); // read those two lines thorough
		std::string nextword, prevword;
		for (int i = 0; std::getline(nextstrstr, nextword, ' '); ++i)
		{
			// std::getline(nextstrstr, nextword, ' ');
			std::getline(prevstrstr, prevword, ' ');
			if (i < 4)
			{ // extract data on positions you are interested in
				prevdata.push_back(std::stod(prevword));
				nextdata.push_back(std::stod(nextword));
			}
		}
		output.reserve(prevdata.size());
		for (int i = 0; i < prevdata.size(); ++i)
			output[i] = auxiliaries::interpolate({prevdata[3], nextdata[3]}, {prevdata[i], nextdata[i]},
												 auxiliaries::InterpolationMode::kLinear, nbar, true);
	}
	else
	{
		// extract data between
		// for densities, pressure and baryonic density, we introduce slight slope for monotony; other entries get copypasted depending on nbar
		if (nbar > (nbar_core_limit + nbar_crust_limit) / 2.0)
		{
			output = std::vector<double>({0.32265273667935968, 94.214131003471735, 945.36865713473674, 9.9999913289570197E-2, 55.307255336247522, 7.1273906305570442E-4, 9.9287083488028366E-2, 7.1282980154182770E-4});
			// I choose slopes by hand : split 10%/80%/10%
			output[0] -= 2.0 * (nbar_core_limit - nbar) / (nbar_core_limit - nbar_crust_limit) * (0.32265273667935968 - 0.32178668643178193) / 10.0;
			output[1] -= 2.0 * (nbar_core_limit - nbar) / (nbar_core_limit - nbar_crust_limit) * (94.214131003471735 - 94.023277698669276) / 10.0;
			output[2] -= 2.0 * (nbar_core_limit - nbar) / (nbar_core_limit - nbar_crust_limit) * (945.36865713473674 - 945.35998787187043) / 10.0;
			output[3] -= 2.0 * (nbar_core_limit - nbar) / (nbar_core_limit - nbar_crust_limit) * (9.9999913289570197E-2 - 9.9798029952044190E-2) / 10.0;
		}
		else
		{
			output = std::vector<double>({0.32178668643178193, 94.023277698669276, 945.35998787187043, 9.9798029952044190E-2});
			// I choose slopes by hand : split 10%/80%/10%
			output[0] += 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (0.32265273667935968 - 0.32178668643178193) / 10.0;
			output[1] += 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (94.214131003471735 - 94.023277698669276) / 10.0;
			output[2] += 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (945.36865713473674 - 945.35998787187043) / 10.0;
			output[3] -= 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (9.9999913289570197E-2 - 9.9798029952044190E-2) / 10.0;
		}
	}
}

std::vector<double> eos_reader::eos_data(const std::vector<double> &input, eos_reader::EOS eos, std::ifstream &fstr)
{
	std::vector<double> result; // empty placeholder for output array

	switch (eos) // choose EoS
	{
	case EOS::APR4:
		apr4(input, result, fstr);
		break;
	case EOS::IST:
		ist_for_ns(input, result, fstr);
		break;
	default:
		throw std::runtime_error("Unknown eos; Encountered in eos_reader::eos_data");
	}

	fstr.clear();
	fstr.seekg(0, std::ios::beg); // reset fstream to initial state for further usage
	return result;				  // output
}