#ifndef TRUESKILL_HPP
#define TRUESKILL_HPP

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <memory>
#include <stdexcept>
#include <numeric>
#include <functional>
#include <iterator>

#include "mathematics.hpp"
#include "factors.hpp"

namespace trueskill {

const double DEFAULT_MU = 25.0;
const double DEFAULT_SIGMA = DEFAULT_MU / 3.0;
const double DEFAULT_BETA = DEFAULT_SIGMA / 2.0;
const double DEFAULT_TAU = DEFAULT_SIGMA / 100.0;
const double DEFAULT_DRAW_PROBABILITY = 0.1;
const double DEFAULT_DELTA = 1e-4;
const double DEFAULT_EPSILON = 1e-6;

double calc_draw_probability(double draw_margin, int size, double beta) {
    return 2 * trueskill::math::cdf(draw_margin / (std::sqrt(size) * beta)) - 1;
}

double calc_draw_margin(double draw_probability, int size, double beta) {
    return beta * std::sqrt(size) * trueskill::math::ppf((draw_probability + 1) / 2);
}

std::vector<int> _team_sizes(const std::vector<std::vector<Rating>>& rating_groups) {
    std::vector<int> team_sizes;
    team_sizes.push_back(0);
    for (const auto& group : rating_groups) {
        team_sizes.push_back(group.size() + team_sizes.back());
    }
    team_sizes.erase(team_sizes.begin());
    return team_sizes;
}

class Rating : public math::Gaussian {
public:
    Rating() : math::Gaussian(DEFAULT_MU, DEFAULT_SIGMA) {}
    Rating(double mu, double sigma) : math::Gaussian(mu, sigma) {}
    
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
    TrueSkill(
        double mu = DEFAULT_MU, 
        double sigma = DEFAULT_SIGMA, 
        double beta = DEFAULT_BETA, 
        double tau = DEFAULT_TAU, 
        double draw_probability = DEFAULT_DRAW_PROBABILITY
    ) : _mu(mu), _sigma(sigma), _beta(beta), _tau(tau), _draw_probability(draw_probability) {}

    Rating createRating(double mu = DEFAULT_MU, double sigma = DEFAULT_SIGMA) {
        if (mu == DEFAULT_MU) {
            mu = _mu;
        }
        if (sigma == DEFAULT_SIGMA) {
            sigma = _sigma;
        }
        return Rating(mu, sigma);
    }

    double v_win(double diff, double draw_margin) {
        double x = diff - draw_margin;
        double denom = trueskill::math::cdf(x);
        return (std::abs(denom) > DEFAULT_EPSILON) ? (trueskill::math::pdf(x) / denom) : -x;
    }

    double v_draw(double diff, double draw_margin) {
        double abs_diff = std::abs(diff);
        double a = draw_margin - abs_diff;
        double b = -draw_margin - abs_diff;
        double denom = trueskill::math::cdf(a) - trueskill::math::cdf(b);
        double numer = trueskill::math::pdf(b) - trueskill::math::pdf(a);
        double f1 = (std::abs(denom) > DEFAULT_EPSILON) ? (numer / denom) : a;
        double f2 = (diff < 0) ? -1 : +1;
        return f1 * f2;
    }

    double w_win(double diff, double draw_margin) {
        double x = diff - draw_margin;
        double v = v_win(diff, draw_margin);
        double w = v * (v + x);
        if (0 < w && w < 1) {
            return w;
        }
        throw std::runtime_error("w_win: w is out of range");
    }

    double w_draw(double diff, double draw_margin) {
        double abs_diff = std::abs(diff);
        double a = draw_margin - abs_diff;
        double b = -draw_margin - abs_diff;
        double denom = trueskill::math::cdf(a) - trueskill::math::cdf(b);
        if (std::abs(denom) < DEFAULT_EPSILON) {
            throw std::runtime_error("w_draw: denominator is too small");
        }
        double v = v_draw(abs_diff, draw_margin);
        return v * v + (a * trueskill::math::pdf(a) - b * trueskill::math::pdf(b)) / denom;
    }

    std::vector<std::vector<Rating>> rate(const std::vector<std::vector<Rating>>& rating_groups, 
                                          const std::vector<int>& ranks = std::vector<int>(),
                                          const std::vector<std::vector<double>>& weights = std::vector<std::vector<double>>(),
                                          double min_delta = DEFAULT_DELTA) {
        auto validated_rating_groups = validate_rating_groups(rating_groups);
        auto validated_weights = validate_weights(weights, validated_rating_groups);
        
        size_t group_size = validated_rating_groups.size();
        std::vector<int> default_ranks(group_size);
        std::iota(default_ranks.begin(), default_ranks.end(), 0);
        const auto& used_ranks = ranks.empty() ? default_ranks : ranks;

        if (used_ranks.size() != group_size) {
            throw std::invalid_argument("Wrong ranks");
        }

        // Sort rating groups
        std::vector<size_t> sorting_indices(group_size);
        std::iota(sorting_indices.begin(), sorting_indices.end(), 0);
        std::sort(sorting_indices.begin(), sorting_indices.end(),
                  [&used_ranks](size_t i1, size_t i2) { return used_ranks[i1] < used_ranks[i2]; });

        std::vector<std::vector<Rating>> sorted_rating_groups;
        std::vector<int> sorted_ranks;
        std::vector<std::vector<double>> sorted_weights;

        for (size_t i : sorting_indices) {
            sorted_rating_groups.push_back(validated_rating_groups[i]);
            sorted_ranks.push_back(used_ranks[i]);
            std::vector<double> adjusted_weights;
            for (double w : validated_weights[i]) {
                adjusted_weights.push_back(std::max(min_delta, w));
            }
            sorted_weights.push_back(adjusted_weights);
        }

        // Build factor graph
        auto graph = build_factor_graph(sorted_rating_groups, sorted_ranks, sorted_weights);
        run_schedule(graph, min_delta);

        // Process results
        std::vector<std::vector<Rating>> result(group_size);
        for (size_t i = 0; i < group_size; ++i) {
            const auto& team = sorted_rating_groups[i];
            std::vector<Rating> new_ratings;
            for (const auto& rating : team) {
                new_ratings.emplace_back(rating.mu(), rating.sigma());
            }
            result[sorting_indices[i]] = std::move(new_ratings);
        }

        return result;
    }

    std::tuple<std::function<void()>, std::function<void()>, 
               std::function<void()>, std::function<void()>, 
               std::function<void()>> 
    factor_graph_builders(const std::vector<std::vector<Rating>>& rating_groups, 
                          const std::vector<int>& ranks, 
                          const std::vector<std::vector<double>>& weights) {
        std::vector<double> flatten_ratings;
        std::vector<double> flatten_weights;
        
        for (const auto& group : rating_groups) {
            flatten_ratings.insert(flatten_ratings.end(), group.begin(), group.end());
        }

        for (const auto& group : weights) {
            flatten_weights.insert(flatten_weights.end(), group.begin(), group.end());
        }

        int size = flatten_ratings.size();
        int group_size = rating_groups.size();
        
        // 创建变量
        std::vector<Variable> rating_vars(size);
        std::vector<Variable> perf_vars(size);
        std::vector<Variable> team_perf_vars(group_size);
        std::vector<Variable> team_diff_vars(group_size - 1);
        std::vector<int> team_sizes = _team_sizes(rating_groups);

        auto build_rating_layer = [&]() {
            for (size_t i = 0; i < size; ++i) {
                yield PriorFactor(rating_vars[i], flatten_ratings[i], this->tau);
            }
        };

        auto build_perf_layer = [&]() {
            for (size_t i = 0; i < size; ++i) {
                yield LikelihoodFactor(rating_vars[i], perf_vars[i], this->beta * this->beta);
            }
        };

        auto build_team_perf_layer = [&]() {
            for (int team = 0; team < group_size; ++team) {
                int start = (team > 0) ? team_sizes[team - 1] : 0;
                int end = team_sizes[team];
                std::vector<Variable> child_perf_vars(perf_vars.begin() + start, perf_vars.begin() + end);
                std::vector<double> coeffs(flatten_weights.begin() + start, flatten_weights.begin() + end);
                yield SumFactor(team_perf_vars[team], child_perf_vars, coeffs);
            }
        };

        auto build_team_diff_layer = [&]() {
            for (int team = 0; team < team_diff_vars.size(); ++team) {
                yield SumFactor(team_diff_vars[team], 
                                std::vector<Variable>{team_perf_vars[team], team_perf_vars[team + 1]}, 
                                std::vector<double>{1.0, -1.0});
            }
        };

        auto build_trunc_layer = [&]() {
            for (size_t x = 0; x < team_diff_vars.size(); ++x) {
                double draw_probability_value = (draw_probability) ? draw_probability() : 0.0;
                int size = rating_groups[x].size() + rating_groups[x + 1].size();
                double draw_margin = calc_draw_margin(draw_probability_value, size, this);
                // 处理排名相等的情况
                std::function<double()> v_func = (ranks[x] == ranks[x + 1]) ? v_draw : v_win;
                std::function<double()> w_func = (ranks[x] == ranks[x + 1]) ? w_draw : w_win;
                yield TruncateFactor(team_diff_vars[x], v_func, w_func, draw_margin);
            }
        };

        // 返回生成器
        return std::make_tuple(build_rating_layer, build_perf_layer, build_team_perf_layer, build_team_diff_layer, build_trunc_layer);
    }

    std::vector<std::vector<Factor*>> run_schedule(std::function<void()> build_rating_layer, 
                                                    std::function<void()> build_perf_layer, 
                                                    std::function<void()> build_team_perf_layer, 
                                                    std::function<void()> build_team_diff_layer, 
                                                    std::function<void()> build_trunc_layer, 
                                                    double min_delta) {
        if (min_delta <= 0) {
            throw std::invalid_argument("min_delta must be greater than 0");
        }
        std::vector<std::vector<Factor*>> layers;

        auto build = [&](std::vector<std::function<void()>> builders) {
            std::vector<std::vector<Factor*>> layers_built;
            for (auto& builder : builders) {
                layers_built.push_back(builder());
            }
            layers.insert(layers.end(), layers_built.begin(), layers_built.end());
            return layers_built;
        };

        // 构建层
        auto layers_built = build({ build_rating_layer, build_perf_layer, build_team_perf_layer });
        auto rating_layer = layers_built[0];
        auto perf_layer = layers_built[1];
        auto team_perf_layer = layers_built[2];

        for (auto& f : rating_layer) {
            f->down();
        }
        // 处理后续箭头
        auto team_diff_layer = build_team_diff_layer();
        auto trunc_layer = build_trunc_layer();
        int team_diff_len = team_diff_layer.size();
        
        for (int x = 0; x < 10; ++x) {
            if (team_diff_len == 1) {
                team_diff_layer[0]->down();
                double delta = trunc_layer[0]->up();
            } else {
                double delta = 0;
                for (int x = 0; x < team_diff_len - 1; ++x) {
                    team_diff_layer[x]->down();
                    delta = std::max(delta, trunc_layer[x]->up());
                    team_diff_layer[x]->up(1);
                }
                for (int x = team_diff_len - 1; x > 0; --x) {
                    team_diff_layer[x]->down();
                    delta = std::max(delta, trunc_layer[x]->up());
                    team_diff_layer[x]->up(0);
                }
            }
            if (delta <= min_delta) {
                break;
            }
        }

        // 处理剩余的箭头
        team_diff_layer[0]->up(0);
        team_diff_layer[team_diff_len - 1]->up(1);
        
        for (auto& f : team_perf_layer) {
            for (size_t x = 0; x < f->vars.size() - 1; ++x) {
                f->up(x);
            }
        }
        for (auto& f : perf_layer) {
            f->up();
        }
        return layers;
    }

    double calculateWinProbability(const Rating& r1, const Rating& r2) {
        double deltaMu = r1.mu() - r2.mu();
        double sqrtSigma = std::sqrt(2 * std::pow(_beta, 2) + std::pow(r1.sigma(), 2) + std::pow(r2.sigma(), 2));
        return trueskill::math::cdf(deltaMu / sqrtSigma);
    }

private:
    double _mu;
    double _sigma;
    double _beta;
    double _tau;
    double _draw_probability;

    std::vector<std::vector<Rating>> validate_rating_groups(const std::vector<std::vector<Rating>>& rating_groups) {
        if (rating_groups.size() < 2) {
            throw std::invalid_argument("Need multiple rating groups");
        }
        if (std::any_of(rating_groups.begin(), rating_groups.end(), 
                        [](const auto& group) { return group.empty(); })) {
            throw std::invalid_argument("Each group must contain multiple ratings");
        }
        return rating_groups;
    }

    std::vector<std::vector<double>> validate_weights(
        const std::vector<std::vector<double>>& weights,
        const std::vector<std::vector<Rating>>& rating_groups) {
        if (weights.empty()) {
            std::vector<std::vector<double>> default_weights;
            for (const auto& group : rating_groups) {
                default_weights.push_back(std::vector<double>(group.size(), 1.0));
            }
            return default_weights;
        }
        return weights;
    }

    math::Matrix build_rotated_a_matrix(const std::vector<std::vector<Rating>>& rating_groups,
                                        const std::vector<double>& flatten_weights) {
        size_t height = rating_groups.size() - 1;
        size_t width = 0;
        for (const auto& group : rating_groups) {
            width += group.size();
        }

        math::Matrix matrix(height, width);
        size_t t = 0;
        for (size_t r = 0; r < height; ++r) {
            const auto& cur = rating_groups[r];
            const auto& next = rating_groups[r + 1];
            
            for (size_t x = t; x < t + cur.size(); ++x) {
                matrix(r, x) = flatten_weights[x];
            }
            t += cur.size();
            
            for (size_t x = t; x < t + next.size(); ++x) {
                matrix(r, x) = -flatten_weights[x];
            }
        }

        return matrix;
    }

    double expose(Rating rating) {
        double k = _mu / _sigma;
        return rating.mu() - k * rating.sigma();
    }

    // Add other necessary private helper functions here
};

} // namespace trueskill

#endif // TRUESKILL_HPP
