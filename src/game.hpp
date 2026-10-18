#ifndef GAME_HPP
#define GAME_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <numeric>
#include <stdexcept>
#include <vector>

#include "player.hpp"

namespace trueskill {

enum class Observation {
    Ordinal,
    Continuous,
    Discrete
};

class Game {
public:
    using Team = std::vector<Player>;
    using Teams = std::vector<Team>;
    using Distributions = std::vector<std::vector<math::Gaussian>>;
    using Weights = std::vector<std::vector<double>>;

    explicit Game(
        const Teams& teams,
        const std::vector<double>& result = {},
        double draw_probability = 0.0,
        const Weights& weights = {},
        Observation observation = Observation::Ordinal)
        : teams(teams),
          result(normalize_result(result, teams.size())),
          weights(normalize_weights(weights, teams)),
          draw_probability(draw_probability),
          observation(observation) {
        validate();
        order = sorted_order();
        team_performances = build_team_performances();
        differences = build_differences();
        build_outcomes();
        likelihoods = compute_likelihoods();
    }

    Distributions posteriors() const {
        Distributions result;
        result.reserve(teams.size());
        for (std::size_t team = 0; team < teams.size(); ++team) {
            result.emplace_back();
            result.back().reserve(teams[team].size());
            for (std::size_t player = 0; player < teams[team].size(); ++player) {
                result.back().push_back(teams[team][player].prior * likelihoods[team][player]);
            }
        }
        return result;
    }

    static math::Gaussian fixed_point_approx(
        double score,
        double mu,
        double sigma,
        std::size_t maximum_iterations = 16,
        double tolerance = 1e-6) {
        return fixed_point_approximation(
            score, mu, sigma, maximum_iterations, tolerance);
    }

    Teams teams;
    std::vector<double> result;
    Weights weights;
    double draw_probability;
    Observation observation;
    double evidence = 1.0;
    Distributions likelihoods;

private:
    struct TeamMessages {
        explicit TeamMessages(const math::Gaussian& prior) : prior(prior) {}

        math::Gaussian posterior_win() const {
            return prior * likelihood_lose;
        }

        math::Gaussian posterior_lose() const {
            return prior * likelihood_win;
        }

        math::Gaussian likelihood() const {
            return likelihood_win * likelihood_lose;
        }

        math::Gaussian prior;
        math::Gaussian likelihood_lose;
        math::Gaussian likelihood_win;
    };

    struct DifferenceMessages {
        explicit DifferenceMessages(const math::Gaussian& prior) : prior(prior) {}

        math::Gaussian prior;
        math::Gaussian likelihood;
    };

    static std::vector<double> normalize_result(
        const std::vector<double>& result,
        std::size_t team_count) {
        if (!result.empty()) {
            return result;
        }

        std::vector<double> default_result(team_count);
        for (std::size_t i = 0; i < team_count; ++i) {
            default_result[i] = static_cast<double>(team_count - i - 1);
        }
        return default_result;
    }

    static Weights normalize_weights(const Weights& weights, const Teams& teams) {
        if (!weights.empty()) {
            return weights;
        }

        Weights default_weights;
        default_weights.reserve(teams.size());
        for (const auto& team : teams) {
            default_weights.emplace_back(team.size(), 1.0);
        }
        return default_weights;
    }

    void validate() const {
        if (teams.size() < 2) {
            throw std::invalid_argument("a game requires at least two teams");
        }
        if (result.size() != teams.size()) {
            throw std::invalid_argument("result size must match team count");
        }
        if (weights.size() != teams.size()) {
            throw std::invalid_argument("weights size must match team count");
        }
        if (draw_probability < 0 || draw_probability > 1) {
            throw std::invalid_argument("draw probability must be in range [0, 1]");
        }
        for (std::size_t team = 0; team < teams.size(); ++team) {
            if (teams[team].empty()) {
                throw std::invalid_argument("teams must not be empty");
            }
            if (weights[team].size() != teams[team].size()) {
                throw std::invalid_argument("each weight group must match its team size");
            }
        }
    }

    std::vector<std::size_t> sorted_order() const {
        std::vector<std::size_t> indices(teams.size());
        std::iota(indices.begin(), indices.end(), 0);
        std::stable_sort(indices.begin(), indices.end(), [this](std::size_t lhs, std::size_t rhs) {
            return result[lhs] > result[rhs];
        });
        return indices;
    }

    Distributions individual_performances() const {
        Distributions performances;
        performances.reserve(teams.size());
        for (std::size_t team = 0; team < teams.size(); ++team) {
            performances.emplace_back();
            performances.back().reserve(teams[team].size());
            for (std::size_t player = 0; player < teams[team].size(); ++player) {
                performances.back().push_back(
                    teams[team][player].performance() * weights[team][player]);
            }
        }
        return performances;
    }

    std::vector<TeamMessages> build_team_performances() const {
        const auto individuals = individual_performances();
        std::vector<TeamMessages> performances;
        performances.reserve(teams.size());
        for (const auto& team : individuals) {
            math::Gaussian total(0, 0);
            for (const auto& performance : team) {
                total = total + performance;
            }
            performances.emplace_back(total);
        }
        return performances;
    }

    std::vector<DifferenceMessages> build_differences() const {
        std::vector<DifferenceMessages> result;
        result.reserve(teams.size() - 1);
        for (std::size_t i = 0; i + 1 < order.size(); ++i) {
            result.emplace_back(
                team_performances[order[i]].prior - team_performances[order[i + 1]].prior);
        }
        return result;
    }

    void build_outcomes() {
        ties.reserve(differences.size());
        margins.reserve(differences.size());
        for (std::size_t i = 0; i < differences.size(); ++i) {
            const std::size_t first = order[i];
            const std::size_t second = order[i + 1];
            ties.push_back(
                observation == Observation::Ordinal && result[first] == result[second]);

            double variance = 0;
            for (const auto& player : teams[first]) {
                variance += player.beta * player.beta;
            }
            for (const auto& player : teams[second]) {
                variance += player.beta * player.beta;
            }
            const double standard_deviation = std::sqrt(variance);
            margins.push_back(observation != Observation::Ordinal || draw_probability == 0
                ? 0
                : std::abs(math::ppf(
                    0.5 - draw_probability / 2.0, 0, standard_deviation)));
        }
    }

    static math::Gaussian approximate_truncation(
        const math::Gaussian& prior,
        double margin,
        bool tie) {
        const double mu = prior.mu();
        const double sigma = prior.sigma();
        double v = 0;
        double w = 0;

        if (!tie) {
            const double alpha = (margin - mu) / sigma;
            const double denominator = math::cdf(-alpha);
            v = math::pdf(-alpha) / denominator;
            w = v * (v - alpha);
        } else {
            const double alpha = (-margin - mu) / sigma;
            const double beta = (margin - mu) / sigma;
            const double denominator = math::cdf(beta) - math::cdf(alpha);
            v = (math::pdf(alpha) - math::pdf(beta)) / denominator;
            const double u = (alpha * math::pdf(alpha) - beta * math::pdf(beta))
                / denominator;
            w = -(u - v * v);
        }

        return math::Gaussian(
            mu + sigma * v,
            sigma * std::sqrt(std::max(0.0, 1.0 - w)));
    }

    void update_evidence(std::size_t index) {
        const auto& prior = differences[index].prior;
        const double observed_difference =
            result[order[index]] - result[order[index + 1]];
        switch (observation) {
        case Observation::Ordinal:
            if (ties[index]) {
                evidence *= math::cdf(margins[index], prior.mu(), prior.sigma())
                    - math::cdf(-margins[index], prior.mu(), prior.sigma());
            } else {
                evidence *= 1.0 - math::cdf(margins[index], prior.mu(), prior.sigma());
            }
            break;
        case Observation::Continuous:
            evidence *= math::pdf(observed_difference, prior.mu(), prior.sigma());
            break;
        case Observation::Discrete:
            evidence *= discrete_evidence(
                observed_difference, prior.mu(), prior.sigma());
            break;
        }
    }

    math::Gaussian difference_likelihood(std::size_t index) const {
        const double observed_difference =
            result[order[index]] - result[order[index + 1]];
        switch (observation) {
        case Observation::Ordinal:
            return approximate_truncation(
                differences[index].prior, margins[index], ties[index])
                / differences[index].prior;
        case Observation::Continuous:
            return math::Gaussian(observed_difference, 0);
        case Observation::Discrete:
            return fixed_point_approximation(
                observed_difference,
                differences[index].prior.mu(),
                differences[index].prior.sigma())
                / differences[index].prior;
        }
        throw std::logic_error("unknown observation model");
    }

    static math::Gaussian fixed_point_approximation(
        double score,
        double mu,
        double sigma,
        std::size_t maximum_iterations = 16,
        double tolerance = 1e-6) {
        const double variance = sigma * sigma;
        if (variance == 0) {
            return math::Gaussian(mu, 0);
        }
        double kappa = 1;
        for (std::size_t iteration = 0; iteration < maximum_iterations; ++iteration) {
            const double term = kappa - mu - score * variance - 1;
            const double root = std::sqrt(term * term + 2 * variance);
            const double numerator = mu + score * variance - 1 - kappa + root;
            const double next = std::log(numerator / (2 * variance));
            if (std::abs(next - kappa) <= tolerance) {
                kappa = next;
                break;
            }
            kappa = next;
        }

        const double rate = std::exp(kappa);
        return math::Gaussian(
            mu + variance * (score - rate),
            std::sqrt(variance / (1 + variance * rate)));
    }

    static double discrete_evidence(double score, double mu, double sigma) {
        if (score < 0 || std::floor(score) != score) {
            return 1e-12;
        }

        constexpr int intervals = 256;
        constexpr double lower = -8;
        constexpr double upper = 8;
        const double step = (upper - lower) / intervals;
        auto integrand = [score, mu, sigma](double z) {
            const double log_rate = mu + sigma * z;
            const double rate = std::exp(std::min(log_rate, 700.0));
            const double log_poisson = score * log_rate - rate - std::lgamma(score + 1);
            return std::exp(log_poisson) * math::pdf(z);
        };

        double sum = integrand(lower) + integrand(upper);
        for (int i = 1; i < intervals; ++i) {
            sum += (i % 2 == 0 ? 2 : 4) * integrand(lower + i * step);
        }
        return std::max(1e-12, sum * step / 3.0);
    }

    std::vector<math::Gaussian> converge_team_likelihoods() {
        for (int iteration = 0; iteration < 5; ++iteration) {
            for (std::size_t edge = 0; edge + 1 < differences.size(); ++edge) {
                differences[edge].prior = team_performances[order[edge]].posterior_win()
                    - team_performances[order[edge + 1]].posterior_lose();
                if (iteration == 0) {
                    update_evidence(edge);
                }
                differences[edge].likelihood = difference_likelihood(edge);
                team_performances[order[edge + 1]].likelihood_lose =
                    team_performances[order[edge]].posterior_win()
                    - differences[edge].likelihood;
            }

            for (std::size_t edge = differences.size() - 1; edge > 0; --edge) {
                differences[edge].prior = team_performances[order[edge]].posterior_win()
                    - team_performances[order[edge + 1]].posterior_lose();
                if (iteration == 0 && edge == differences.size() - 1) {
                    update_evidence(edge);
                }
                differences[edge].likelihood = difference_likelihood(edge);
                team_performances[order[edge]].likelihood_win =
                    team_performances[order[edge + 1]].posterior_lose()
                    + differences[edge].likelihood;
            }
        }

        if (differences.size() == 1) {
            update_evidence(0);
            differences[0].prior = team_performances[order[0]].posterior_win()
                - team_performances[order[1]].posterior_lose();
            differences[0].likelihood = difference_likelihood(0);
        }

        team_performances[order.front()].likelihood_win =
            team_performances[order[1]].posterior_lose() + differences.front().likelihood;
        team_performances[order.back()].likelihood_lose =
            team_performances[order[order.size() - 2]].posterior_win()
            - differences.back().likelihood;

        std::vector<math::Gaussian> result;
        result.reserve(team_performances.size());
        for (const auto& performance : team_performances) {
            result.push_back(performance.likelihood());
        }
        return result;
    }

    Distributions performance_likelihoods() {
        const auto individuals = individual_performances();
        const auto team_likelihoods = converge_team_likelihoods();
        Distributions result;
        result.reserve(teams.size());

        for (std::size_t team = 0; team < teams.size(); ++team) {
            result.emplace_back();
            result.back().reserve(teams[team].size());
            for (std::size_t player = 0; player < teams[team].size(); ++player) {
                if (weights[team][player] == 0) {
                    result.back().emplace_back();
                    continue;
                }

                const double remaining_variance = std::max(
                    0.0,
                    team_performances[team].prior.variance()
                        - individuals[team][player].variance());
                const math::Gaussian team_without_player(
                    team_performances[team].prior.mu() - individuals[team][player].mu(),
                    std::sqrt(remaining_variance));
                result.back().push_back(
                    (team_likelihoods[team] - team_without_player)
                    * (1.0 / weights[team][player]));
            }
        }
        return result;
    }

    Distributions compute_likelihoods() {
        auto result = performance_likelihoods();
        for (std::size_t team = 0; team < teams.size(); ++team) {
            for (std::size_t player = 0; player < teams[team].size(); ++player) {
                result[team][player] = result[team][player]
                    - math::Gaussian(0, teams[team][player].beta);
            }
        }
        return result;
    }

    std::vector<std::size_t> order;
    std::vector<TeamMessages> team_performances;
    std::vector<DifferenceMessages> differences;
    std::vector<bool> ties;
    std::vector<double> margins;
};

} // namespace trueskill

#endif // GAME_HPP
