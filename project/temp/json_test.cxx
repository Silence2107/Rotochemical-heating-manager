
#include "../../include/instantiator.hpp"

#include <iostream>
#include <fstream>
#include <typeinfo>

int main()
{
    instantiator::instantiate_system("../../data/json_test.json");
    std::cout << "pressure_upp: " << instantiator::pressure_upp << std::endl;
}