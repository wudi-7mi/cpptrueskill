#include <iostream>
#include "src/trueskill.hpp"  // 假设 TrueSkill 实现已经封装在这个头文件中

int main() {
    // 假设 TrueSkill 的 Rating 类有一个构造函数接受两个参数 mu 和 sigma
    Rating alice(25.0, 8.0);
    Rating bob(30.0, 3.0);

    // 计算两者比赛的质量
    double quality = quality_1vs1(alice, bob);
    
    if (quality < 0.50) {
        std::cout << "This match seems to be not so fair" << std::endl;
    }

    // 假设 rate_1vs1 返回一个包含更新后 Rating 的 pair
    auto result = rate_1vs1(alice, bob);
    
    alice = result.first;  // 更新 Alice 的 Rating
    bob = result.second;   // 更新 Bob 的 Rating

    // 输出更新后的结果
    std::cout << "Alice's new rating: mu = " << alice.mu() << ", sigma = " << alice.sigma() << std::endl;
    std::cout << "Bob's new rating: mu = " << bob.mu() << ", sigma = " << bob.sigma() << std::endl;
    std::cout << "Match quality: " << quality << std::endl;

    return 0;
}
