#include "../../3rd-party/json/single_include/nlohmann/json.hpp"
#include <iostream>
#include <fstream>
#include <typeinfo>

int main()
{
    using json = nlohmann::json;
    std::ifstream i("../../data/json_test.json");
    if (!i.is_open()) {
        std::cout << "open file failed" << std::endl;
        return -1;
    }
    json j = json::parse(i);
    bool b = j["happy"];
    std::cout << j["pi"] << std::endl;
    std::cout << b << std::endl;
}