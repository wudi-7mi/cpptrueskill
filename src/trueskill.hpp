#ifndef TRUESKILL_HPP
#define TRUESKILL_HPP

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <memory>
#include "factors.hpp"
#include "mathematics.hpp"

namespace trueskill {

class Rating : public trueskill::math::Gaussian {
public:
    Rating(double mu, double sigma) : trueskill::math::Gaussian(mu, sigma) {}
    
    operator int() const {
        return static_cast<int>(mu());
    }

    operator long() const {
        return static_cast<long>(mu());
    }

    operator float() const {
        return static_cast<float>(mu());
    }

    friend std::ostream& operator<<(std::ostream& os, const Rating& r) {
        os << "Rating(mu=" << r.mu() << ", sigma=" << r.sigma() << ")";
        return os;
    }
};

class TrueSkill {
public:
    TrueSkill(double mu, double sigma, double beta, double tau, double draw_probability = 0.6)
        : mu_(mu), sigma_(sigma), beta_(beta), tau_(tau), draw_probability_(draw_probability) {}

    std::vector<std::pair<double, double>> rate(const std::vector<std::vector<Player>>& teams, bool draw = false) {
        std::vector<std::pair<double, double>> new_ratings;
        
        FactorGraph graph;
        
        // 创建变量和因子
        for (const auto& team : teams) {
            for (const auto& player : team) {
                auto var = graph.createVariable(player.mu(), player.sigma());
                graph.addFactor(std::make_shared<PriorFactor>(var, player.mu(), player.sigma()));
            }
        }
        
        // 添加队伍性能因子和比较因子
        for (size_t i = 0; i < teams.size() - 1; ++i) {
            auto team1_perf = graph.createVariable();
            auto team2_perf = graph.createVariable();
            
            graph.addFactor(std::make_shared<PerformanceFactor>(team1_perf, teams[i]));
            graph.addFactor(std::make_shared<PerformanceFactor>(team2_perf, teams[i+1]));
            
            graph.addFactor(std::make_shared<ComparisonFactor>(team1_perf, team2_perf, draw));
        }

        // 执行消息传递算法
        graph.runInference();

        // 更新玩家评分
        for (const auto& var : graph.getVariables()) {
            const auto& msg = var->getValue();
            new_ratings.emplace_back(msg.mu(), msg.sigma());
        }
        
        return new_ratings;
    }

    // 其他辅助函数保持不变
    // ...

    double calculateWinProbability(const Rating& r1, const Rating& r2) {
        double deltaMu = r1.mu() - r2.mu();
        double sqrtSigma = std::sqrt(2 * std::pow(beta_, 2) + std::pow(r1.sigma(), 2) + std::pow(r2.sigma(), 2));
        return trueskill::math::cdf(deltaMu / sqrtSigma);
    }

private:
    double mu_;
    double sigma_;
    double beta_;
    double tau_;
    double draw_probability_;
};

} // namespace trueskill

#endif // TRUESKILL_HPP
