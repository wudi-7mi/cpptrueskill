#ifndef SKILL_HPP
#define SKILL_HPP

#include <algorithm>
#include <cstddef>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "mathematics.hpp"

namespace trueskill {

class Skill {
public:
    Skill() = default;

    explicit Skill(std::size_t event)
        : events{event}, likelihoods{math::Gaussian()} {}

    math::Gaussian posterior() const {
        return forward * backward * likelihood_product();
    }

    math::Gaussian forward_posterior() const {
        return forward * likelihood_product();
    }

    math::Gaussian backward_posterior() const {
        return backward * likelihood_product();
    }

    const math::Gaussian& likelihood(std::size_t event) const {
        return likelihoods[event_position(event)];
    }

    std::pair<double, double> update_likelihood(
        std::size_t event,
        const math::Gaussian& likelihood) {
        const std::size_t position = event_position(event);
        const auto step = likelihood.delta(likelihoods[position]);
        likelihoods[position] = likelihood;
        return step;
    }

    void add_event(std::size_t event) {
        events.push_back(event);
        likelihoods.emplace_back();
    }

    std::vector<std::size_t> events;
    math::Gaussian forward;
    math::Gaussian backward;
    std::vector<math::Gaussian> likelihoods;
    std::optional<math::Gaussian> online;

private:
    std::size_t event_position(std::size_t event) const {
        const auto found = std::find(events.begin(), events.end(), event);
        if (found == events.end()) {
            throw std::out_of_range("event does not belong to this skill");
        }
        return static_cast<std::size_t>(std::distance(events.begin(), found));
    }

    math::Gaussian likelihood_product() const {
        math::Gaussian result;
        for (const auto& likelihood : likelihoods) {
            result = result * likelihood;
        }
        return result;
    }
};

} // namespace trueskill

#endif // SKILL_HPP
