#ifndef HISTORY_HPP
#define HISTORY_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <istream>
#include <iterator>
#include <limits>
#include <map>
#include <numeric>
#include <optional>
#include <ostream>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "game.hpp"
#include "skill.hpp"

namespace trueskill {

class History {
public:
    using PlayerId = std::string;
    using TeamComposition = std::vector<PlayerId>;
    using GameComposition = std::vector<TeamComposition>;
    using Composition = std::vector<GameComposition>;
    using Results = std::vector<std::vector<double>>;
    using Times = std::vector<double>;
    using Priors = std::map<PlayerId, Player>;
    using HistoryWeights = std::vector<Game::Weights>;
    using Observations = std::vector<Observation>;
    using LearningCurve = std::vector<std::pair<double, math::Gaussian>>;
    using LearningCurves = std::map<PlayerId, LearningCurve>;
    using Step = std::pair<double, double>;

    History(
        const Composition& composition,
        const Results& results = {},
        const Times& times = {},
        const Priors& custom_priors = {},
        double mu = 0.0,
        double sigma = 3.0,
        double beta = 1.0,
        double gamma = 0.15,
        double draw_probability = 0.0,
        bool online = false,
        const HistoryWeights& weights = {},
        const Observations& observations = {})
        : mu(mu),
          sigma(sigma),
          beta(beta),
          gamma(gamma),
          draw_probability(draw_probability),
          online_mode(online),
          priors(custom_priors) {
        validate_parameters();
        initialize(composition, results, times, weights, observations);
        revealed_batches = online_mode ? 0 : batches.size();
    }

    Step forward_propagation() {
        std::map<PlayerId, math::Gaussian> last_messages;
        std::map<PlayerId, std::optional<double>> last_times;
        Step maximum{0, 0};

        for (std::size_t batch = 0; batch < revealed_batches; ++batch) {
            receive_messages(batch, true, last_messages, last_times);
            for (std::size_t event = 0; event < batches[batch].size(); ++event) {
                maximum = max_step(maximum, update_event(batch, event));
            }
            send_messages(batch, true, last_messages, last_times);
        }
        return maximum;
    }

    Step backward_propagation() {
        if (revealed_batches == 0) {
            return {0, 0};
        }

        std::map<PlayerId, math::Gaussian> last_messages;
        std::map<PlayerId, std::optional<double>> last_times;
        Step maximum{0, 0};
        send_messages(revealed_batches - 1, false, last_messages, last_times);

        for (std::size_t batch = revealed_batches - 1; batch-- > 0;) {
            receive_messages(batch, false, last_messages, last_times);
            for (std::size_t event = 0; event < batches[batch].size(); ++event) {
                maximum = max_step(maximum, update_event(batch, event));
            }
            send_messages(batch, false, last_messages, last_times);
        }
        return maximum;
    }

    Step iteration() {
        return max_step(forward_propagation(), backward_propagation());
    }

    void unveil_batch() {
        if (revealed_batches < batches.size()) {
            ++revealed_batches;
        }
    }

    std::pair<Step, std::size_t> convergence(
        std::size_t iterations = 8,
        double epsilon = 1e-5) {
        Step final_step{0, 0};
        std::size_t final_iterations = 0;
        bool processed = false;

        do {
            if (online_mode && revealed_batches < batches.size()) {
                ++revealed_batches;
            } else if (!online_mode) {
                revealed_batches = batches.size();
            }

            double delta = std::numeric_limits<double>::infinity();
            final_iterations = 0;
            while (final_iterations < iterations && delta > epsilon) {
                final_step = iteration();
                delta = std::max(final_step.first, final_step.second);
                ++final_iterations;
            }
            processed = true;
        } while (online_mode && revealed_batches < batches.size());

        if (!processed) {
            return {{0, 0}, 0};
        }
        return {final_step, final_iterations};
    }

    void add_history(
        const Composition& composition,
        const Results& results = {},
        const Times& times = {},
        const Priors& custom_priors = {},
        const HistoryWeights& weights = {},
        const Observations& observations = {}) {
        validate_data(composition, results, times, weights, observations);
        for (const auto& prior : custom_priors) {
            priors.emplace(prior.first, prior.second);
        }

        Times normalized_times = times;
        if (normalized_times.empty()) {
            const double start = batch_times.empty() ? 0 : batch_times.back() + 1;
            normalized_times.resize(composition.size());
            for (std::size_t i = 0; i < normalized_times.size(); ++i) {
                normalized_times[i] = start + static_cast<double>(i);
            }
        }

        if (online_mode && revealed_batches > 0 && !normalized_times.empty()
            && *std::min_element(normalized_times.begin(), normalized_times.end())
                <= batch_times[revealed_batches - 1]) {
            throw std::invalid_argument(
                "online history can only append times after the last revealed batch");
        }

        for (std::size_t event = 0; event < composition.size(); ++event) {
            add_event(
                composition[event],
                results.empty() ? std::vector<double>{} : results[event],
                normalized_times[event],
                weights.empty() ? Game::Weights{} : weights[event],
                observations.empty()
                    ? Observation::Ordinal
                    : observations[observations.size() == 1 ? 0 : event]);
        }
        if (!online_mode) {
            revealed_batches = batches.size();
        }
    }

    LearningCurves learning_curves(
        const std::vector<PlayerId>& who = {},
        bool online = false) const {
        const std::set<PlayerId> selected(who.begin(), who.end());
        LearningCurves result;
        for (std::size_t batch = 0; batch < batch_skills.size(); ++batch) {
            for (const auto& entry : batch_skills[batch]) {
                if (!selected.empty() && selected.count(entry.first) == 0) {
                    continue;
                }
                const math::Gaussian estimate = online_mode && online && entry.second.online
                    ? *entry.second.online
                    : entry.second.posterior();
                result[entry.first].emplace_back(batch_times[batch], estimate);
            }
        }
        return result;
    }

    Game::Teams within_priors(
        std::size_t batch,
        std::size_t event,
        bool online = false) const {
        if (batch >= batches.size() || event >= batches[batch].size()) {
            throw std::out_of_range("history batch or event index is out of range");
        }

        Game::Teams teams;
        teams.reserve(batches[batch][event].size());
        for (const auto& composition_team : batches[batch][event]) {
            teams.emplace_back();
            teams.back().reserve(composition_team.size());
            for (const auto& name : composition_team) {
                const Skill& skill = batch_skills[batch].at(name);
                const math::Gaussian prior = online && skill.online
                    ? *skill.online
                    : skill.posterior() / skill.likelihood(event);
                const Player& attributes = priors.at(name);
                teams.back().emplace_back(prior, attributes.beta, attributes.gamma);
            }
        }
        return teams;
    }

    double log_evidence() const {
        double result = 0;
        for (const auto& batch : batch_evidence) {
            for (const auto& evidence : batch) {
                if (!evidence) {
                    continue;
                }
                if (*evidence == 0) {
                    return -std::numeric_limits<double>::infinity();
                }
                result += std::log(*evidence);
            }
        }
        return result;
    }

    double geometric_mean() const {
        std::size_t count = 0;
        for (const auto& batch : batch_evidence) {
            for (const auto& evidence : batch) {
                count += evidence.has_value();
            }
        }
        if (count == 0) {
            return std::numeric_limits<double>::quiet_NaN();
        }
        return std::exp(log_evidence() / static_cast<double>(count));
    }

    void save(std::ostream& output) const {
        static constexpr char magic[] = "CPPTTT1";
        output.write(magic, sizeof(magic));
        write_value(output, mu);
        write_value(output, sigma);
        write_value(output, beta);
        write_value(output, gamma);
        write_value(output, draw_probability);
        write_value(output, online_mode);

        write_size(output, priors.size());
        for (const auto& entry : priors) {
            write_string(output, entry.first);
            write_gaussian(output, entry.second.prior);
            write_value(output, entry.second.beta);
            write_value(output, entry.second.gamma);
        }

        write_size(output, event_count);
        write_size(output, revealed_batches);
        write_size(output, batches.size());
        for (std::size_t batch = 0; batch < batches.size(); ++batch) {
            write_value(output, batch_times[batch]);
            write_size(output, batches[batch].size());
            for (std::size_t event = 0; event < batches[batch].size(); ++event) {
                write_composition(output, batches[batch][event]);
                write_doubles(output, batch_results[batch][event]);
                write_weights(output, batch_weights[batch][event]);
                write_value(output, static_cast<std::uint8_t>(
                    batch_observations[batch][event]));
                write_value(output, batch_evidence[batch][event].has_value());
                if (batch_evidence[batch][event]) {
                    write_value(output, *batch_evidence[batch][event]);
                }
            }

            write_size(output, batch_skills[batch].size());
            for (const auto& entry : batch_skills[batch]) {
                write_string(output, entry.first);
                const Skill& skill = entry.second;
                write_size(output, skill.events.size());
                for (const auto event : skill.events) {
                    write_size(output, event);
                }
                write_gaussian(output, skill.forward);
                write_gaussian(output, skill.backward);
                write_size(output, skill.likelihoods.size());
                for (const auto& likelihood : skill.likelihoods) {
                    write_gaussian(output, likelihood);
                }
                write_value(output, skill.online.has_value());
                if (skill.online) {
                    write_gaussian(output, *skill.online);
                }
            }
        }
        if (!output) {
            throw std::runtime_error("failed to serialize history");
        }
    }

    static History load(std::istream& input) {
        static constexpr char expected_magic[] = "CPPTTT1";
        char magic[sizeof(expected_magic)]{};
        input.read(magic, sizeof(magic));
        if (!input || !std::equal(
                std::begin(magic), std::end(magic), std::begin(expected_magic))) {
            throw std::runtime_error("invalid history serialization header");
        }

        const double mu = read_value<double>(input);
        const double sigma = read_value<double>(input);
        const double beta = read_value<double>(input);
        const double gamma = read_value<double>(input);
        const double draw_probability = read_value<double>(input);
        const bool online = read_value<bool>(input);
        Priors priors;
        const std::size_t prior_count = read_size(input);
        for (std::size_t i = 0; i < prior_count; ++i) {
            const auto name = read_string(input);
            const auto prior = read_gaussian(input);
            const double player_beta = read_value<double>(input);
            const double player_gamma = read_value<double>(input);
            priors.emplace(name, Player(prior, player_beta, player_gamma));
        }

        History history(
            {}, {}, {}, priors, mu, sigma, beta, gamma,
            draw_probability, online);
        history.event_count = read_size(input);
        history.revealed_batches = read_size(input);
        const std::size_t batch_count = read_size(input);
        history.batches.resize(batch_count);
        history.batch_results.resize(batch_count);
        history.batch_times.resize(batch_count);
        history.batch_skills.resize(batch_count);
        history.batch_weights.resize(batch_count);
        history.batch_observations.resize(batch_count);
        history.batch_evidence.resize(batch_count);

        for (std::size_t batch = 0; batch < batch_count; ++batch) {
            history.batch_times[batch] = read_value<double>(input);
            const std::size_t event_count = read_size(input);
            history.batches[batch].reserve(event_count);
            history.batch_results[batch].reserve(event_count);
            history.batch_weights[batch].reserve(event_count);
            history.batch_observations[batch].reserve(event_count);
            history.batch_evidence[batch].reserve(event_count);
            for (std::size_t event = 0; event < event_count; ++event) {
                history.batches[batch].push_back(read_composition(input));
                history.batch_results[batch].push_back(read_doubles(input));
                history.batch_weights[batch].push_back(read_weights(input));
                history.batch_observations[batch].push_back(
                    static_cast<Observation>(read_value<std::uint8_t>(input)));
                const bool has_evidence = read_value<bool>(input);
                history.batch_evidence[batch].push_back(
                    has_evidence
                        ? std::optional<double>(read_value<double>(input))
                        : std::nullopt);
            }

            const std::size_t skill_count = read_size(input);
            for (std::size_t i = 0; i < skill_count; ++i) {
                const auto name = read_string(input);
                Skill skill;
                const std::size_t skill_event_count = read_size(input);
                skill.events.reserve(skill_event_count);
                for (std::size_t j = 0; j < skill_event_count; ++j) {
                    skill.events.push_back(read_size(input));
                }
                skill.forward = read_gaussian(input);
                skill.backward = read_gaussian(input);
                const std::size_t likelihood_count = read_size(input);
                skill.likelihoods.reserve(likelihood_count);
                for (std::size_t j = 0; j < likelihood_count; ++j) {
                    skill.likelihoods.push_back(read_gaussian(input));
                }
                if (read_value<bool>(input)) {
                    skill.online = read_gaussian(input);
                }
                history.batch_skills[batch].emplace(name, std::move(skill));
            }
        }
        if (!input) {
            throw std::runtime_error("failed to deserialize history");
        }
        return history;
    }

    std::size_t size() const {
        return event_count;
    }

    std::size_t batch_count() const {
        return batches.size();
    }

    std::size_t revealed_batch_count() const {
        return revealed_batches;
    }

    double mu;
    double sigma;
    double beta;
    double gamma;
    double draw_probability;
    bool online_mode;
    Priors priors;
    std::vector<std::vector<GameComposition>> batches;
    std::vector<Results> batch_results;
    Times batch_times;
    std::vector<std::map<PlayerId, Skill>> batch_skills;
    std::vector<std::vector<Game::Weights>> batch_weights;
    std::vector<Observations> batch_observations;
    std::vector<std::vector<std::optional<double>>> batch_evidence;

private:
    void validate_parameters() const {
        if (sigma < 0 || beta < 0 || gamma < 0) {
            throw std::invalid_argument("sigma, beta and gamma must be non-negative");
        }
        if (draw_probability < 0 || draw_probability > 1) {
            throw std::invalid_argument("draw probability must be in range [0, 1]");
        }
    }

    static void validate_data(
        const Composition& composition,
        const Results& results,
        const Times& times,
        const HistoryWeights& weights,
        const Observations& observations) {
        if (!results.empty() && results.size() != composition.size()) {
            throw std::invalid_argument("results size must match composition size");
        }
        if (!times.empty() && times.size() != composition.size()) {
            throw std::invalid_argument("times size must match composition size");
        }
        if (!weights.empty() && weights.size() != composition.size()) {
            throw std::invalid_argument("weights size must match composition size");
        }
        if (!observations.empty() && observations.size() != 1
            && observations.size() != composition.size()) {
            throw std::invalid_argument(
                "observations must contain one model or one model per game");
        }
        for (std::size_t event = 0; event < composition.size(); ++event) {
            if (composition[event].size() < 2) {
                throw std::invalid_argument("each game requires at least two teams");
            }
            for (const auto& team : composition[event]) {
                if (team.empty()) {
                    throw std::invalid_argument("teams must not be empty");
                }
            }
            if (!results.empty() && !results[event].empty()
                && results[event].size() != composition[event].size()) {
                throw std::invalid_argument("each result must match its game team count");
            }
            if (!weights.empty() && !weights[event].empty()) {
                if (weights[event].size() != composition[event].size()) {
                    throw std::invalid_argument("each weight set must match its game team count");
                }
                for (std::size_t team = 0; team < composition[event].size(); ++team) {
                    if (weights[event][team].size() != composition[event][team].size()) {
                        throw std::invalid_argument(
                            "each weight group must match its team player count");
                    }
                }
            }
        }
    }

    void initialize(
        const Composition& composition,
        const Results& results,
        const Times& times,
        const HistoryWeights& weights,
        const Observations& observations) {
        validate_data(composition, results, times, weights, observations);
        Times normalized_times = times;
        if (normalized_times.empty()) {
            normalized_times.resize(composition.size());
            std::iota(normalized_times.begin(), normalized_times.end(), 0.0);
        }

        std::vector<std::size_t> order(composition.size());
        std::iota(order.begin(), order.end(), 0);
        std::stable_sort(order.begin(), order.end(), [&normalized_times](std::size_t lhs, std::size_t rhs) {
            return normalized_times[lhs] < normalized_times[rhs];
        });
        for (const std::size_t event : order) {
            add_event(
                composition[event],
                results.empty() ? std::vector<double>{} : results[event],
                normalized_times[event],
                weights.empty() ? Game::Weights{} : weights[event],
                observations.empty()
                    ? Observation::Ordinal
                    : observations[observations.size() == 1 ? 0 : event]);
        }
    }

    void add_event(
        const GameComposition& composition,
        const std::vector<double>& result,
        double time,
        const Game::Weights& weights,
        Observation observation) {
        auto position = std::lower_bound(batch_times.begin(), batch_times.end(), time);
        std::size_t batch = static_cast<std::size_t>(std::distance(batch_times.begin(), position));
        if (position == batch_times.end() || *position != time) {
            batch_times.insert(position, time);
            batches.insert(batches.begin() + batch, decltype(batches)::value_type{});
            batch_results.insert(
                batch_results.begin() + batch, decltype(batch_results)::value_type{});
            batch_skills.insert(
                batch_skills.begin() + batch, decltype(batch_skills)::value_type{});
            batch_weights.insert(
                batch_weights.begin() + batch, decltype(batch_weights)::value_type{});
            batch_observations.insert(
                batch_observations.begin() + batch,
                decltype(batch_observations)::value_type{});
            batch_evidence.insert(
                batch_evidence.begin() + batch, decltype(batch_evidence)::value_type{});
        }

        const std::size_t event = batches[batch].size();
        batches[batch].push_back(composition);
        batch_results[batch].push_back(result);
        batch_weights[batch].push_back(weights);
        batch_observations[batch].push_back(observation);
        batch_evidence[batch].push_back(std::nullopt);
        ++event_count;

        for (const auto& team : composition) {
            for (const auto& name : team) {
                ensure_prior(name);
                auto found = batch_skills[batch].find(name);
                if (found == batch_skills[batch].end()) {
                    batch_skills[batch].emplace(name, Skill(event));
                } else {
                    found->second.add_event(event);
                }
            }
        }
    }

    Player& ensure_prior(const PlayerId& name) {
        const auto found = priors.find(name);
        if (found != priors.end()) {
            return found->second;
        }
        return priors.emplace(
            name,
            Player(math::Gaussian(mu, sigma), beta, gamma)).first->second;
    }

    math::Gaussian initial_message(const PlayerId& name) const {
        const auto found = priors.find(name);
        return found == priors.end() ? math::Gaussian(mu, sigma) : found->second.prior;
    }

    void receive_messages(
        std::size_t batch,
        bool forward,
        std::map<PlayerId, math::Gaussian>& last_messages,
        std::map<PlayerId, std::optional<double>>& last_times) {
        for (auto& entry : batch_skills[batch]) {
            const PlayerId& name = entry.first;
            const auto time_found = last_times.find(name);
            const double elapsed = time_found == last_times.end() || !time_found->second
                ? 0
                : std::abs(batch_times[batch] - *time_found->second);
            const double player_gamma = priors.at(name).gamma;
            const double drift = std::min(
                std::sqrt(elapsed * player_gamma * player_gamma),
                1.67 * sigma);
            const auto message_found = last_messages.find(name);
            const math::Gaussian previous = message_found == last_messages.end()
                ? (forward ? initial_message(name) : math::Gaussian())
                : message_found->second;
            const math::Gaussian received = previous + math::Gaussian(0, drift);

            if (forward) {
                entry.second.forward = received;
                if (online_mode && !entry.second.online) {
                    entry.second.online = received;
                }
            } else {
                entry.second.backward = received;
            }
        }
    }

    void send_messages(
        std::size_t batch,
        bool forward,
        std::map<PlayerId, math::Gaussian>& last_messages,
        std::map<PlayerId, std::optional<double>>& last_times) const {
        for (const auto& entry : batch_skills[batch]) {
            last_times[entry.first] = batch_times[batch];
            last_messages[entry.first] = forward
                ? entry.second.forward_posterior()
                : entry.second.backward_posterior();
        }
    }

    Step update_event(std::size_t batch, std::size_t event) {
        Game game(
            within_priors(batch, event),
            batch_results[batch][event],
            draw_probability,
            batch_weights[batch][event],
            batch_observations[batch][event]);
        if (!online_mode || !batch_evidence[batch][event]) {
            batch_evidence[batch][event] = game.evidence;
        }

        Step maximum{0, 0};
        for (std::size_t team = 0; team < batches[batch][event].size(); ++team) {
            for (std::size_t player = 0; player < batches[batch][event][team].size(); ++player) {
                const auto& name = batches[batch][event][team][player];
                maximum = max_step(
                    maximum,
                    batch_skills[batch].at(name).update_likelihood(
                        event, game.likelihoods[team][player]));
            }
        }
        return maximum;
    }

    static Step max_step(const Step& lhs, const Step& rhs) {
        return {std::max(lhs.first, rhs.first), std::max(lhs.second, rhs.second)};
    }

    template <typename T>
    static void write_value(std::ostream& output, const T& value) {
        output.write(reinterpret_cast<const char*>(&value), sizeof(T));
    }

    template <typename T>
    static T read_value(std::istream& input) {
        T value{};
        input.read(reinterpret_cast<char*>(&value), sizeof(T));
        return value;
    }

    static void write_size(std::ostream& output, std::size_t size) {
        write_value(output, static_cast<std::uint64_t>(size));
    }

    static std::size_t read_size(std::istream& input) {
        return static_cast<std::size_t>(read_value<std::uint64_t>(input));
    }

    static void write_string(std::ostream& output, const std::string& value) {
        write_size(output, value.size());
        output.write(value.data(), static_cast<std::streamsize>(value.size()));
    }

    static std::string read_string(std::istream& input) {
        std::string result(read_size(input), '\0');
        input.read(result.data(), static_cast<std::streamsize>(result.size()));
        return result;
    }

    static void write_gaussian(std::ostream& output, const math::Gaussian& value) {
        write_value(output, value.mu());
        write_value(output, value.sigma());
    }

    static math::Gaussian read_gaussian(std::istream& input) {
        const double mu = read_value<double>(input);
        const double sigma = read_value<double>(input);
        return math::Gaussian(mu, sigma);
    }

    static void write_doubles(std::ostream& output, const std::vector<double>& values) {
        write_size(output, values.size());
        for (const double value : values) {
            write_value(output, value);
        }
    }

    static std::vector<double> read_doubles(std::istream& input) {
        std::vector<double> result(read_size(input));
        for (double& value : result) {
            value = read_value<double>(input);
        }
        return result;
    }

    static void write_weights(std::ostream& output, const Game::Weights& weights) {
        write_size(output, weights.size());
        for (const auto& team : weights) {
            write_doubles(output, team);
        }
    }

    static Game::Weights read_weights(std::istream& input) {
        Game::Weights result;
        const std::size_t size = read_size(input);
        result.reserve(size);
        for (std::size_t i = 0; i < size; ++i) {
            result.push_back(read_doubles(input));
        }
        return result;
    }

    static void write_composition(
        std::ostream& output,
        const GameComposition& composition) {
        write_size(output, composition.size());
        for (const auto& team : composition) {
            write_size(output, team.size());
            for (const auto& name : team) {
                write_string(output, name);
            }
        }
    }

    static GameComposition read_composition(std::istream& input) {
        GameComposition result;
        const std::size_t team_count = read_size(input);
        result.reserve(team_count);
        for (std::size_t team = 0; team < team_count; ++team) {
            result.emplace_back();
            const std::size_t player_count = read_size(input);
            result.back().reserve(player_count);
            for (std::size_t player = 0; player < player_count; ++player) {
                result.back().push_back(read_string(input));
            }
        }
        return result;
    }

    std::size_t event_count = 0;
    std::size_t revealed_batches = 0;
};

} // namespace trueskill

#endif // HISTORY_HPP
