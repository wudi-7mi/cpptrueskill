#include <iostream>
#include "src/trueskill.hpp"

int main() {

    // 1 vs 1
    std::cout << "[Example] 1 vs 1" << std::endl;

    trueskill::Rating p1(80, 2.5);
    trueskill::Rating p2(75, 2.5);
    double quality1v1 = trueskill::quality_1vs1(p1, p2);
    std::cout << "quality1v1: " << quality1v1 << std::endl;

    auto result1v1 = trueskill::rate_1vs1(p1, p2);

    for (const auto& team : result1v1) {
        for (const auto& rating : team) {
            std::cout << rating << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    std::cout << "If match is drawn" << std::endl;
    auto result1v1_drawn = trueskill::rate_1vs1(p1, p2, true);
    for (const auto& team : result1v1_drawn) {
        for (const auto& rating : team) {
            std::cout << rating << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;


    // 2 vs 2 vs 3
    std::cout << "[Example] 2 vs 2 vs 3 with weights" << std::endl;

    trueskill::TrueSkill ts(25.0, 8.333333, 4.2, 0.3, 0.1);
    std::cout << ts << std::endl;

    trueskill::Rating alice(34.23, 7.22);
    trueskill::Rating bob(24.23, 3.22);   
    trueskill::Rating charlie(12.23, 3.9);
    trueskill::Rating david(26.23, 6.11);
    trueskill::Rating eve(41.2, 4.22);
    trueskill::Rating frank(24.23, 5.08);
    trueskill::Rating gloria(4.23, 3.08);

    std::vector<std::vector<double>> weights = {{1.0, 1.0}, {0.8, 1.0}, {1.0, 1.0, 1.0}};
    std::vector<std::vector<trueskill::Rating>> teams = {{alice, bob}, {charlie, david}, {eve, frank, gloria}};

    double quality = ts.quality(teams, weights);
    std::cout << "quality: " << quality << std::endl;

    auto newTeams = ts.rate(teams, {2, 1, 3}, weights);
    for (const auto& team : newTeams) {
        for (const auto& rating : team) {
            std::cout << rating << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    return 0;
}
