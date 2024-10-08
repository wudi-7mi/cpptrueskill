#ifndef TRUESKILL_HPP
#define TRUESKILL_HPP

#include <vector>
#include <cmath>
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
        std::vector<Rating> flatten_ratings;
        std::vector<double> flatten_weights;
        for (const auto& group : sorted_rating_groups) {
            flatten_ratings.insert(flatten_ratings.end(), group.begin(), group.end());
        }

        for (const auto& group : sorted_weights) {
            flatten_weights.insert(flatten_weights.end(), group.begin(), group.end());
        }

        size_t size = flatten_ratings.size();
        size_t rating_group_size = rating_groups.size();

        std::vector<Variable> rating_vars(size);
        std::vector<Variable> perf_vars(size);
        std::vector<Variable> team_perf_vars(rating_group_size);
        std::vector<Variable> team_diff_vars(rating_group_size - 1);
        std::vector<int> team_sizes = _team_sizes(rating_groups);

        std::vector<PriorFactor*> rating_layer;
        for (size_t i = 0; i < size; ++i) {
            auto* pf = new PriorFactor(&rating_vars[i], flatten_ratings[i], _tau);
            rating_layer.push_back(pf);
        }

        std::vector<LikelihoodFactor*> perf_layer;
        for (size_t i = 0; i < size; ++i) {
            auto* lf = new LikelihoodFactor(&rating_vars[i], &perf_vars[i], _beta * _beta);
            perf_layer.push_back(lf);
        }

        std::vector<SumFactor*> team_perf_layer;
        for (size_t team = 0; team < rating_group_size; ++team) {
            size_t start = (team > 0) ? team_sizes[team - 1] : 0;
            size_t end = team_sizes[team];
            std::vector<Variable*> child_perf_vars;
            for (size_t i = start; i < end; ++i) {
                child_perf_vars.push_back(&perf_vars[i]);
            }
            std::vector<double> coeffs;
            for (size_t i = start; i < end; ++i) {
                coeffs.push_back(flatten_weights[i]);
            }
            auto* sf = new SumFactor(&team_perf_vars[team], child_perf_vars, coeffs);
            team_perf_layer.push_back(sf);
        }

        std::vector<SumFactor*> team_diff_layer;
        for (size_t team = 0; team < rating_group_size - 1; ++team) {
            auto* sf = new SumFactor(&team_diff_vars[team], {&team_perf_vars[team], &team_perf_vars[team + 1]}, {1, -1});
            team_diff_layer.push_back(sf);
        }

        std::vector<TruncateFactor*> trunc_layer;
        for (size_t x = 0; x < team_diff_vars.size(); ++x) {
            int rg_size = rating_groups[x].size() + rating_groups[x + 1].size();
            double draw_margin = calc_draw_margin(_draw_probability, rg_size, _beta);
            auto v_func = (sorted_ranks[x] == sorted_ranks[x + 1]) ? v_draw : v_win;
            auto w_func = (sorted_ranks[x] == sorted_ranks[x + 1]) ? w_draw : w_win;
            auto* tf = new TruncateFactor(&team_diff_vars[x], v_func, w_func, draw_margin);
            trunc_layer.push_back(tf);
        }

        for (auto& f : rating_layer) {
            f->down();
        }
        for (auto& f : perf_layer) {
            f->down();
        }
        for (auto& f : team_perf_layer) {
            f->down();
        }

        size_t team_diff_len = team_diff_layer.size();
        for (int x = 0; x < 10; ++x) {
            double delta = 0;
            if (team_diff_len == 1) {
                team_diff_layer[0]->down();
                delta = trunc_layer[0]->up();
            } else {
                for (size_t y = 0; y < team_diff_len - 1; ++y) {
                    team_diff_layer[y]->down();
                    delta = std::max(delta, trunc_layer[y]->up());
                    team_diff_layer[y]->up(1);
                }
                for (size_t y = team_diff_len - 1; y > 0; --y) {
                    team_diff_layer[y]->down();
                    delta = std::max(delta, trunc_layer[y]->up());
                    team_diff_layer[y]->up(0);
                }
            }

            if (delta <= min_delta) {
                break;
            }
        }

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

        auto sorted_team_sizes = _team_sizes(sorted_rating_groups);

        // Create transformed groups
        std::vector<std::vector<Rating>> transformed_groups;
        for (size_t i = 0, start = 0; i < sorted_team_sizes.size(); ++i) {
            size_t end = sorted_team_sizes[i];
            std::vector<Rating> group;
            for (size_t j = start; j < end; ++j) {
                group.push_back(Rating(static_cast<double>(rating_layer[j]->vars[0]->mu()),
                                    static_cast<double>(rating_layer[j]->vars[0]->sigma())));
            }
            transformed_groups.push_back(group);
            start = end;
        }

        // Sort the groups based on the original sorting hint
        std::vector<std::tuple<int, std::vector<Rating>>> unsorted_groups;
        for (size_t i = 0; i < sorting_indices.size(); ++i) {
            unsorted_groups.push_back(std::make_tuple(sorting_indices[i], transformed_groups[i]));
        }

        // Sorting function based on the first element (original order)
        auto by_hint = [](const std::tuple<int, std::vector<Rating>>& a,
                        const std::tuple<int, std::vector<Rating>>& b) {
            return std::get<0>(a) < std::get<0>(b);
        };

        std::sort(unsorted_groups.begin(), unsorted_groups.end(), by_hint);

        std::vector<std::vector<Rating>> result;
        for (const auto& item : unsorted_groups) {
            result.push_back(std::get<1>(item));
        }
        return result;
    }

    double quality(const std::vector<std::vector<Rating>>& rating_groups, 
                  const std::vector<std::vector<double>>& weights = std::vector<std::vector<double>>()) {
        auto validated_rating_groups = validate_rating_groups(rating_groups);
        auto validated_weights = validate_weights(weights, validated_rating_groups);

        std::vector<Rating> flatten_ratings;
        for (const auto& group : validated_rating_groups) {
            flatten_ratings.insert(flatten_ratings.end(), group.begin(), group.end());
        }
        std::vector<double> flatten_weights;
        for (const auto& group : validated_weights) {
            flatten_weights.insert(flatten_weights.end(), group.begin(), group.end());
        }

        size_t length = flatten_ratings.size();

        auto mean_matrix = create_mean_matrix(flatten_ratings);
        auto variance_matrix = create_variance_matrix(flatten_ratings);
        auto rotated_a_matrix = create_rotated_a_matrix(validated_rating_groups, flatten_weights);
        auto a_matrix = rotated_a_matrix.transpose();

        math::Matrix _ata = (rotated_a_matrix * std::pow(_beta, 2)) * a_matrix;
        math::Matrix _atsa = rotated_a_matrix * variance_matrix * a_matrix;
        math::Matrix start = mean_matrix.transpose() * a_matrix;
        math::Matrix middle = _ata + _atsa;
        math::Matrix end = rotated_a_matrix * mean_matrix;

        auto m_earg = (start * middle.inverse() * end * -0.5);
        auto e_arg = m_earg.determinant();
        auto s_arg = _ata.determinant() / middle.determinant();

        return std::exp(e_arg) * std::sqrt(s_arg);
    }

    double calculateWinProbability(const Rating& r1, const Rating& r2) {
        double deltaMu = r1.mu() - r2.mu();
        double sqrtSigma = std::sqrt(2 * std::pow(_beta, 2) + std::pow(r1.sigma(), 2) + std::pow(r2.sigma(), 2));
        return trueskill::math::cdf(deltaMu / sqrtSigma);
    }

    double expose(Rating rating) {
        double k = _mu / _sigma;
        return rating.mu() - k * rating.sigma();
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

    math::Matrix create_mean_matrix(const std::vector<Rating>& ratings) {
            int length = ratings.size();
            math::Matrix mean_matrix(length, 1);
            for (int i = 0; i < length; ++i) {
                mean_matrix(i, 0) = ratings[i].mu();
            }
            return mean_matrix;
        }

    math::Matrix create_variance_matrix(const std::vector<Rating>& ratings) {
        int length = ratings.size();
        math::Matrix variance_matrix(length, length);
        for (int i = 0; i < length; ++i) {
            variance_matrix(i, i) = std::pow(ratings[i].sigma(), 2);
        }
        return variance_matrix;
    }

    math::Matrix create_rotated_a_matrix(const std::vector<std::vector<Rating>>& rating_groups,
                                        const std::vector<double>& flatten_weights) {
        int height = rating_groups.size() - 1;
        int width = 0;
        for (const auto& group : rating_groups) {
            width += group.size();
        }

        math::Matrix matrix(height, width);
        int t = 0;
        for (int r = 0; r < height; ++r) {
            const auto& cur = rating_groups[r];
            const auto& next = rating_groups[r + 1];
            
            for (int x = t; x < t + cur.size(); ++x) {
                matrix(r, x) = flatten_weights[x];
            }
            t += cur.size();

            for (int x = t; x < t + next.size(); ++x) {
                matrix(r, x) = -flatten_weights[x];
            }
        }

        return matrix;
    }

    // Add other necessary private helper functions here
};

double quality_1vs1(const Rating& r1, const Rating& r2) {
    TrueSkill ts;
    return ts.quality({{r1}, {r2}});
}

std::vector<std::vector<Rating>> rate_1vs1(const Rating& r1, const Rating& r2) {
    TrueSkill ts;
    return ts.rate({{r1}, {r2}});
}

} // namespace trueskill

#endif // TRUESKILL_HPP
