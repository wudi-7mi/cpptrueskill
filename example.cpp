#include <iostream>
#include "src/trueskill.hpp"

int main() {
    trueskill::Rating alice(80, 2.5);
    trueskill::Rating bob(75, 2.5);
    std::cout << "alice: " << alice << std::endl;
    std::cout << "bob: " << bob << std::endl;
    double quality = trueskill::quality_1vs1(alice, bob);
    std::cout << "quality: " << quality << std::endl;
    auto result = trueskill::rate_1vs1(alice, bob);
    for (const auto& team : result) {
        for (const auto& rating : team) {
            std::cout << rating << " ";
        }
        std::cout << std::endl;
    }
}
