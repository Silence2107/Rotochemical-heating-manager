
#include "../../include/auxiliaries.h"

#include "../../3rd-party/argparse/argparse.hpp"
#include "../../3rd-party/json/single_include/nlohmann/json.hpp"

#include <vector>
#include <algorithm>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>

int main(int argc, char **argv)
{
    argparse::ArgumentParser parser("bridge_maxwell_construction", "Ensure recognition by RHM of a Maxwell-like construction at the H-QGP boundary (provides \"patched\" EoS)", "Argparse powered by SiLeader");

    parser.addArgument({"--inputfile"}, "json input file path (required)");
    parser.addArgument({"--pressure_fraction"}, "(Pmid-Ph) / (Pq-Ph) at the boundary, where Pmid is the transition pressure (optional, default: midpoint)");
    parser.addArgument({"--in_place"}, "whether to rewrite original EoS with modified EoS (optional, value-free, default: only console output)", argparse::ArgumentType::StoreTrue);
    parser.addArgument({"--n_digits"}, "how many digits to keep at printing (optional, value-free, default: 10)");
    auto args = parser.parseArgs(argc, argv);
    double pressure_frac = args.safeGet<double>("pressure_fraction", 0.5);
    bool in_place_rewrite = args.has("in_place");
    size_t n_digits = args.safeGet<size_t>("n_digits", 10);

    // this is non-standard RHM program, that tampers with EoS file. Instatiator protects against this and hence is not appropriate here
    // relevant information : eos datatable (path, rows, columns), pressure column, nbar column, uquark nbar column

    // START copypaste

    // eos table
    std::ifstream i(args.get<std::string>("inputfile"));
    if (!i.is_open())
    {
        RHM_THROW(std::runtime_error, "UI inputfile requested, but the path is invalid.");
    }
    using json = nlohmann::json;
    json j = json::parse(i);

    auto eos_datafile = j["EoSSetup"]["Datafile"]["Path"];
    if (!(eos_datafile.is_string()))
        RHM_THROW(std::runtime_error, "UI error: Datafile path must be provided as a string.");
    auto eos_datafile_rows = j["EoSSetup"]["Datafile"]["Rows"];
    if (!(eos_datafile_rows.size() == 2))
    {
        if (eos_datafile_rows.is_null())
            eos_datafile_rows = {0, 0};
        else
            RHM_THROW(std::runtime_error, "UI error: Datafile rows must be a pair-array.");
    }
    else if (!(eos_datafile_rows[0].is_number_integer() && eos_datafile_rows[1].is_number_integer()))
        RHM_THROW(std::runtime_error, "UI error: Datafile rows must be provided as integers.");

    auto eos_datafile_cols = j["EoSSetup"]["Datafile"]["Columns"];
    if (!(eos_datafile_cols.size() == 2))
    {
        if (eos_datafile_cols.is_null())
            eos_datafile_cols = {0, 0};
        else
            RHM_THROW(std::runtime_error, "UI error: Datafile cols must be a pair-array.");
    }
    else if (!(eos_datafile_cols[0].is_number_integer() && eos_datafile_cols[1].is_number_integer()))
        RHM_THROW(std::runtime_error, "UI error: Datafile cols must be provided as integers.");

    auto table = auxiliaries::io::read_tabulated_file(eos_datafile, eos_datafile_cols, eos_datafile_rows);

    // nbar column
    auto nbar_column = j["EoSSetup"]["Quantities"]["BarionicDensity"]["Column"];
    if (!(nbar_column.is_number_integer()))
        RHM_THROW(std::runtime_error, "UI error: Barionic density column number must be provided as an integer.");

    // pressure column
    auto pressure_column = j["EoSSetup"]["Quantities"]["Pressure"]["Column"];
    if (!(pressure_column.is_number_integer()))
        RHM_THROW(std::runtime_error, "UI error: Pressure column number must be provided as an integer.");
    
    // uquark number density column
    auto uquark_density_column_read = j["EoSSetup"]["Quantities"]["NumberDensities"]["Uquark"]["Column"];
    if (!(uquark_density_column_read.is_number_integer()))
        RHM_THROW(std::runtime_error, "UI error: Uquark density column number must be provided as an integer, which is required for PT determination.");

    // END copypaste

    // locate PT

    // is the EoS df descending/ascending?
    auto nbars = table[nbar_column],
         uquark_densities = table[uquark_density_column_read],
         pressures = table[pressure_column];
    if(nbars.size() < 4)
        RHM_THROW(std::runtime_error, "EoS datafile must have plenty of rows to perform meaningful Maxwell brinding.");
    bool ascending = (nbars[0] < nbars[1]);
    if (!ascending)
    {
        std::reverse(nbars.begin(), nbars.end());
        std::reverse(uquark_densities.begin(), uquark_densities.end());
        std::reverse(pressures.begin(), pressures.end());
    }

    size_t last_H = 0, first_Q = 0;
    for (size_t row = 0; row < uquark_densities.size(); ++row)
    {
        bool entered_q_phase = (uquark_densities[row] != 0.0);
        if (entered_q_phase && row < 2)
        {
            std::cout << "HP thin or nonexistent, Maxwell construction redundant\n";
            return 0;
        }
        if(entered_q_phase)
        {
            last_H = row - 1;
            break;
        }
        if (!entered_q_phase && row > uquark_densities.size() - 2)
        {
            std::cout << "QP thin or nonexistent, Maxwell construction redundant\n";
            return 0;
        }
    }
    first_Q = last_H + 1;

    if(pressures[first_Q] == pressures[last_H])
    {
        std::cout << "Maxwell-like construction is already at place\n";
        return 0;
    }

    // transition in (P, n)
    double p_transition = pressure_frac * (pressures[first_Q] - pressures[last_H]) + pressures[last_H];

    // linearly extrapolate hadronic and quark branches
    double k_h_pres = (pressures[last_H] - pressures[last_H - 1]) / (nbars[last_H] - nbars[last_H - 1]),
           k_q_pres = (pressures[first_Q + 1] - pressures[first_Q]) / (nbars[first_Q + 1] - nbars[first_Q]);

    double nbar_new_H_limit, nbar_new_Q_limit;
    if (k_h_pres * k_q_pres == 0)
    {
        // handle as large derivative discontinuity
        nbar_new_H_limit = nbars[last_H];
        nbar_new_Q_limit = nbars[first_Q];
    }
    else
    {
        nbar_new_H_limit = (p_transition - pressures[last_H]) / k_h_pres + nbars[last_H];
        nbar_new_Q_limit = (p_transition - pressures[first_Q]) / k_q_pres + nbars[first_Q];
    }

    //rewrite EoS with two additional lines at PT
    auto modified_table = table;
    for (size_t col = 0; col < table.size(); ++col) 
    {
        // linearly extrapolate every quantity
        auto quantities = modified_table[col];
        if (!ascending)
            std::reverse(quantities.begin(), quantities.end());
        double k_h = (quantities[last_H] - quantities[last_H - 1]) / (nbars[last_H] - nbars[last_H - 1]),
               k_q = (quantities[first_Q + 1] - quantities[first_Q]) / (nbars[first_Q + 1] - nbars[first_Q]);

        double quantity_new_H_limit = quantities[last_H] + k_h * (nbar_new_H_limit - nbars[last_H]),
                quantity_new_Q_limit = quantities[first_Q] + k_q * (nbar_new_Q_limit - nbars[first_Q]);
        // properly insert them into modified table
        if (ascending)
            modified_table[col].insert(modified_table[col].begin() + first_Q, {quantity_new_H_limit, quantity_new_Q_limit});
        else
            modified_table[col].insert(modified_table[col].end() - first_Q, {quantity_new_Q_limit, quantity_new_H_limit});
        }
    size_t indent = n_digits + 10;
    std::cout << std::scientific << std::setprecision(n_digits);
    for(size_t row = 0; row < modified_table[0].size(); ++row)
    {
        std::cout << std::left;
        for (size_t col = 0; col < modified_table.size(); ++col)
        {
            std::cout << std::setw(indent) << modified_table[col][row];
        }
        std::cout << '\n';
    }
    if(in_place_rewrite)
    {
        std::ofstream out(eos_datafile);
        out << std::scientific << std::setprecision(n_digits);
        for (size_t row = 0; row < modified_table[0].size(); ++row)
        {
            out << std::left;
            for (size_t col = 0; col < modified_table.size(); ++col)
            {
                out << std::setw(indent) << modified_table[col][row];
            }
            out << '\n';
        }
    }
}