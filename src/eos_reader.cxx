
#include "../include/eos_reader.h"
#include "../include/constants.h"
#include "../include/auxiliaries.h"

#include <exception>
#include <string>
#include <sstream>
#include <functional>
#include <fstream>
#include <vector>

std::vector<double> eos_reader::predefined::apr4(const std::vector<double> &input, std::ifstream &fstr)
{
	using namespace constants::apr4;
	std::vector<double> output;
	double nbar = input[0];					// barionic density (input[0])
	if (nbar > nbar_upp || nbar < nbar_low) // we do not have data beyond these values
		throw std::runtime_error("Data request out of range; Encountered in eos_reader::apr4");
	else if (nbar > nbar_core_limit)
	{
		std::string nextline, prevline;
		for (int i = 0; std::getline(fstr, nextline) && i < 6; ++i)
			; // skip 6 lines
		while (std::getline(fstr >> std::ws, nextline))
		{															 // read line by line, boundary whitespace free
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
			output.push_back(auxiliaries::interpolate({prevdata[2], nextdata[2]}, {prevdata[i], nextdata[i]}, auxiliaries::InterpolationMode::kLinear, nbar, true));
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
			output.push_back(auxiliaries::interpolate({prevdata[2], nextdata[2]}, {prevdata[i], nextdata[i]}, auxiliaries::InterpolationMode::kLinear, nbar, true));
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
	return output;
}

std::vector<double> eos_reader::predefined::apr4_cached(std::vector<std::vector<double>> &cache, const std::vector<double> &input, std::ifstream &fstr, const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator)
{
	if (cache.empty())
	{
		std::string nextline;
		size_t line_number = 0;
		for (; std::getline(fstr, nextline) && line_number < 6; ++line_number)
			; // skip 6 lines
		while (std::getline(fstr >> std::ws, nextline) && line_number < 235)
		{
			nextline = auxiliaries::retrieve_cleared_line(nextline); // clear line
			std::stringstream strstr(nextline);
			std::string word;
			std::vector<double> data;
			for (int i = 0; std::getline(strstr, word, ' '); ++i)
				data.push_back(std::stod(word));
			cache.push_back(data);
			++line_number;
		} // read all data before line 236 into cache
	}
	using namespace constants::apr4;
	std::vector<double> output;
	double nbar = input[0];					// barionic density (input[0])
	if (nbar > nbar_upp || nbar < nbar_low) // we do not have data beyond these values
		throw std::runtime_error("Data request out of range; Encountered in eos_reader::apr4");
	output.reserve(cache.front().size()); // output size estimate
	std::vector<double> x_interp(cache.size()), y_interp(cache.size());
	for (int j = 0; j < cache.size(); ++j)
		x_interp[j] = cache[j][2]; // array of baryonic densities
	if (nbar > nbar_core_limit)
	{
		for (int i = 0; i < 7; ++i) // 7 variables to extract
		{
			for (int j = 0; j < cache.size(); ++j)
				y_interp[j] = cache[j][i]; // array of data
			output.push_back(interpolator(x_interp, y_interp, nbar));
		}
	}
	else if (nbar < nbar_crust_limit)
	{
		for (int i = 0; i < 6; ++i) // 6 variables to extract
		{
			for (int j = 0; j < cache.size(); ++j)
				y_interp[j] = cache[j][i]; // array of data
			output.push_back(interpolator(x_interp, y_interp, nbar));
		}
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
	return output;
}

std::vector<double> eos_reader::predefined::ist_for_ns(const std::vector<double> &input, std::ifstream &fstr)
{
	using namespace constants::ist_ns;
	std::vector<double> output;
	double nbar = input[0];					// barionic density (input[0])
	if (nbar > nbar_upp || nbar < nbar_low) // we do not have data beyond these values
		throw std::runtime_error("Data request out of range; Encountered in eos_reader::ist_for_ns");
	else if (nbar > nbar_core_limit)
	{
		std::string nextline, prevline;
		while (std::getline(fstr >> std::ws, nextline))
		{															 // read line by line, boundary whitespace free
			nextline = auxiliaries::retrieve_cleared_line(nextline); // clear line
			std::stringstream strstr(nextline);
			std::string word;
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' ');
			std::getline(strstr, word, ' '); // skip three first values so that to compare input[0] with barionic density
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
			output.push_back(auxiliaries::interpolate({prevdata[3], nextdata[3]}, {prevdata[i], nextdata[i]},
													  auxiliaries::InterpolationMode::kLinear, nbar, true));
	}
	else if (nbar < nbar_core_limit)
	{
		std::string nextline, prevline;
		while (std::getline(fstr >> std::ws, nextline))
		{															 // read line by line, boundary whitespace free
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
			output.push_back(auxiliaries::interpolate({prevdata[3], nextdata[3]}, {prevdata[i], nextdata[i]},
													  auxiliaries::InterpolationMode::kLinear, nbar, true));
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
			output[3] += 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (9.9999913289570197E-2 - 9.9798029952044190E-2) / 10.0;
		}
	}
	return output;
}

std::vector<double> eos_reader::predefined::ist_for_ns_cached(std::vector<std::vector<double>> &cache, const std::vector<double> &input, std::ifstream &fstr, const std::function<double(const std::vector<double> &, const std::vector<double> &, double)> &interpolator)
{
	if (cache.empty())
	{
		std::string nextline;
		while (std::getline(fstr >> std::ws, nextline))
		{
			nextline = auxiliaries::retrieve_cleared_line(nextline); // clear line
			std::stringstream strstr(nextline);
			std::string word;
			std::vector<double> data;
			for (int i = 0; std::getline(strstr, word, ' '); ++i)
				data.push_back(std::stod(word));
			cache.push_back(data);
		} // read all data into cache
	}
	using namespace constants::ist_ns;
	std::vector<double> output;
	double nbar = input[0];					// barionic density (input[0])
	if (nbar > nbar_upp || nbar < nbar_low) // we do not have data beyond these values
		throw std::runtime_error("Data request out of range; Encountered in eos_reader::apr4");
	output.reserve(cache.front().size()); // output size estimate
	std::vector<double> x_interp(cache.size()), y_interp(cache.size());
	for (int j = 0; j < cache.size(); ++j)
		x_interp[j] = cache[j][3]; // array of baryonic densities
	if (nbar > nbar_core_limit)
	{
		for (int i = 0; i < 8; ++i) // 8 variables to extract
		{
			for (int j = 0; j < cache.size(); ++j)
				y_interp[j] = cache[j][i]; // array of data
			output.push_back(interpolator(x_interp, y_interp, nbar));
		}
	}
	else if (nbar < nbar_crust_limit)
	{
		for (int i = 0; i < 4; ++i) // 4 variables to extract
		{
			for (int j = 0; j < cache.size(); ++j)
				y_interp[j] = cache[j][i]; // array of data
			output.push_back(interpolator(x_interp, y_interp, nbar));
		}
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
			output[3] += 2.0 * (nbar - nbar_crust_limit) / (nbar_core_limit - nbar_crust_limit) * (9.9999913289570197E-2 - 9.9798029952044190E-2) / 10.0;
		}
	}
	return output;
}

std::vector<double> eos_reader::eos_data(const std::vector<double> &input, eos_reader::predefined::EOS eos, std::ifstream &fstr)
{
	std::vector<double> result; // empty placeholder for output array

	using namespace eos_reader::predefined;
	switch (eos) // choose EoS
	{
	case EOS::APR4:
		result = apr4(input, fstr);
		break;
	case EOS::IST:
		result = ist_for_ns(input, fstr);
		break;
	default:
		throw std::runtime_error("Unknown eos; Encountered in eos_reader::eos_data");
	}

	fstr.clear();
	fstr.seekg(0, std::ios::beg); // reset fstream to initial state for further usage
	return result;				  // output
}

std::vector<double> eos_reader::eos_data(const std::vector<double> &input, const std::function<std::vector<double>(const std::vector<double> &, std::ifstream &)> &eos, std::ifstream &fstr)
{
	std::vector<double> result; // empty placeholder for output array

	result = eos(input, fstr);

	fstr.clear();
	fstr.seekg(0, std::ios::beg); // reset fstream to initial state for further usage
	return result;				  // output
}
