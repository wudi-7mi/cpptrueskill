#ifndef PLAYER_HPP
#define PLAYER_HPP

#include <ostream>
#include <stdexcept>

#include "mathematics.hpp"

namespace trueskill {

inline constexpr double TTT_DEFAULT_MU = 0.0;
inline constexpr double TTT_DEFAULT_SIGMA = 6.0;
inline constexpr double TTT_DEFAULT_BETA = 1.0;
inline constexpr double TTT_DEFAULT_GAMMA = 0.03;

class Player {
public:
    explicit Player(
        const math::Gaussian& prior = math::Gaussian(TTT_DEFAULT_MU, TTT_DEFAULT_SIGMA),
        double beta = TTT_DEFAULT_BETA,
        double gamma = TTT_DEFAULT_GAMMA)
        : prior(prior), beta(beta), gamma(gamma) {
        if (beta < 0) {
            throw std::invalid_argument("beta should be greater than or equal to 0");
        }
        if (gamma < 0) {
            throw std::invalid_argument("gamma should be greater than or equal to 0");
        }
    }

    math::Gaussian performance() const {
        return prior + math::Gaussian(0, beta);
    }

    friend std::ostream& operator<<(std::ostream& os, const Player& player) {
        os << "Player(prior=" << player.prior
           << ", beta=" << player.beta
           << ", gamma=" << player.gamma << ")";
        return os;
    }

    math::Gaussian prior;
    double beta;
    double gamma;
};

} // namespace trueskill

#endif // PLAYER_HPP
