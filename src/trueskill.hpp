#ifndef TRUESKILL_HPP
#define TRUESKILL_HPP

#include <vector>
#include <cmath>
#include <random>
#include <algorithm>
#include <memory>

#include "mathematics.hpp"

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

std::vector<int> _team_sizes(const std::vector<std::vector<int>>& rating_groups) {
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
};

} // namespace trueskill

#endif // TRUESKILL_HPP
