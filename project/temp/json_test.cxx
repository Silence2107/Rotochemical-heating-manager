
#include "../../include/inputfile.hpp"

#include <iostream>
#include <fstream>
#include <typeinfo>

int main()
{
    inputfile::instantiate_system("../../data/json_test.json");
    std::cout << "pressure_upp: " << inputfile::pressure_upp << std::endl;
}