#include <doctest/doctest.h>

#include <sstream>

#include "../src/history.hpp"

using namespace trueskill;
using namespace trueskill::math;

TEST_CASE("History performs forward filtering") {
    const History::Composition composition = {
        {{"a"}, {"b"}},
        {{"a"}, {"c"}},
        {{"b"}, {"c"}}
    };
    const History::Results results = {{1, 0}, {0, 1}, {1, 0}};
    History::Priors priors;
    for (const auto* name : {"a", "b", "c"}) {
        priors.emplace(name, Player(Gaussian(25, 25.0 / 3.0), 25.0 / 6.0, 1.25));
    }

    History history(composition, results, {1, 2, 3}, priors);
    history.forward_propagation();
    const auto curves = history.learning_curves();
    CHECK(curves.at("a")[0].second.mu() == doctest::Approx(29.205).epsilon(1e-3));
    CHECK(curves.at("a")[0].second.sigma() == doctest::Approx(7.19448).epsilon(1e-3));
}

TEST_CASE("History converges across time") {
    const History::Composition composition = {
        {{"a"}, {"b"}},
        {{"a"}, {"c"}},
        {{"b"}, {"c"}}
    };
    const History::Results results = {{1, 0}, {0, 1}, {1, 0}};
    History history(composition, results, {}, {}, 0, 6, 1, 0.05);
    history.convergence(14, 0);

    const auto curves = history.learning_curves();
    CHECK(curves.at("a")[0].second.mu() == doctest::Approx(0.000548).epsilon(1e-4));
    CHECK(curves.at("a")[0].second.sigma() == doctest::Approx(2.395715).epsilon(1e-4));
    CHECK(curves.at("b")[0].second.mu() == doctest::Approx(-0.001641).epsilon(1e-4));
    CHECK(curves.at("b").back().second.mu() == doctest::Approx(0.001762).epsilon(1e-4));
}

TEST_CASE("History groups simultaneous games") {
    const History::Composition composition = {
        {{"a"}, {"b"}},
        {{"b"}, {"c"}},
        {{"c"}, {"a"}}
    };
    History history(composition, {{1, 0}, {1, 0}, {1, 0}}, {1, 1, 1}, {},
        25, 25.0 / 3.0, 25.0 / 6.0, 1.25);
    history.convergence(10, 1e-4);

    CHECK(history.batch_count() == 1);
    const auto curves = history.learning_curves();
    CHECK(curves.at("a")[0].second.mu() == doctest::Approx(25).epsilon(1e-3));
    CHECK(curves.at("a")[0].second.sigma() == doctest::Approx(5.419).epsilon(1e-3));
}

TEST_CASE("History applies elapsed-time skill drift") {
    const History::Composition composition = {
        {{"a"}, {"b"}},
        {{"a"}, {"b"}}
    };
    History history(composition, {{1, 0}, {1, 0}}, {0, 1}, {}, 0, 6, 1, 10);
    history.forward_propagation();

    CHECK(history.batch_skills[1].at("a").forward.mu()
        == doctest::Approx(3.3390790626));
    CHECK(history.batch_skills[1].at("a").forward.sigma()
        == doctest::Approx(11.1736543));
}

TEST_CASE("History supports incremental online processing") {
    const History::Composition first = {
        {{"a"}, {"b"}},
        {{"b"}, {"c"}}
    };
    const History::Composition second = {
        {{"c"}, {"a"}}
    };
    History history({}, {}, {}, {}, 0, 3, 1, 0.15, 0, true);
    history.add_history(first);
    history.convergence(8, 0);
    CHECK(history.revealed_batch_count() == 2);

    history.add_history(second);
    history.convergence(8, 0);
    CHECK(history.size() == 3);
    CHECK(history.revealed_batch_count() == 3);
    CHECK(history.learning_curves({"a"}, true).at("a").size() == 2);
    CHECK(std::isfinite(history.log_evidence()));
    CHECK(std::isfinite(history.geometric_mean()));
}

TEST_CASE("History supports mixed observation models") {
    History::Priors priors = {
        {"a1", Player(Gaussian(2, 2), 1, 0.1)},
        {"a2", Player(Gaussian(2, 2), 1, 0.1)},
        {"b", Player(Gaussian(4, 2), 1, 0.1)},
        {"c", Player(Gaussian(3, 2), 1, 0.1)}
    };
    const History::GameComposition game = {{"a1", "a2"}, {"c"}, {"b"}};
    const History::Composition composition = {game, game, game};
    const History::Results results = {{4, 0, 2}, {4.0, 0.0, 2.1}, {4, 0, 2}};
    const History::Observations observations = {
        Observation::Ordinal,
        Observation::Continuous,
        Observation::Discrete
    };
    History history(
        composition, results, {0, 1, 2}, priors,
        0, 3, 1, 0.15, 0, false, {}, observations);
    history.forward_propagation();

    const auto curves = history.learning_curves();
    CHECK(curves.at("a1")[0].second.is_approx(Gaussian(3.053757, 1.782038), 1e-5));
    CHECK(curves.at("a1")[1].second.is_approx(Gaussian(2.996293, 1.480994), 1e-5));
    CHECK(curves.at("a1")[2].second.is_approx(Gaussian(2.413075, 1.269865), 1e-5));
    CHECK(curves.at("c")[0].second.is_approx(Gaussian(1.925404, 1.686565), 1e-5));
}

TEST_CASE("History serialization preserves converged incremental state") {
    const History::Composition first = {
        {{"a"}, {"b"}},
        {{"b"}, {"c"}}
    };
    History original(first, {}, {}, {}, 0, 3, 1, 0.15, 0, true);
    original.convergence(8, 0);

    std::stringstream stream(std::ios::in | std::ios::out | std::ios::binary);
    original.save(stream);
    stream.seekg(0);
    History restored = History::load(stream);

    const History::Composition second = {{{"c"}, {"a"}}};
    original.add_history(second);
    restored.add_history(second);
    original.convergence(8, 0);
    restored.convergence(8, 0);

    const auto expected = original.learning_curves();
    const auto actual = restored.learning_curves();
    CHECK(actual.size() == expected.size());
    for (const auto& player : expected) {
        CHECK(actual.at(player.first).size() == player.second.size());
        for (std::size_t i = 0; i < player.second.size(); ++i) {
            CHECK(actual.at(player.first)[i].second.is_approx(player.second[i].second));
        }
    }
    CHECK(restored.log_evidence() == doctest::Approx(original.log_evidence()));
}

TEST_CASE("History exposes online batches and validates event shapes") {
    History history(
        {{{"a"}, {"b"}}, {{"b"}, {"c"}}},
        {}, {}, {}, 0, 3, 1, 0.15, 0, true);
    CHECK(history.revealed_batch_count() == 0);
    history.unveil_batch();
    CHECK(history.revealed_batch_count() == 1);
    history.forward_propagation();
    CHECK(history.learning_curves({"a"}).at("a").size() == 1);

    CHECK_THROWS_AS(
        History(History::Composition{{{"a"}}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        History(History::Composition{{{"a"}, {"b"}}}, {{1, 2, 3}}),
        std::invalid_argument);
    CHECK_THROWS_AS(
        History(
            History::Composition{{{"a"}, {"b"}}}, {}, {}, {},
            0, 3, 1, 0.15, 0, false, {{{1}, {1, 2}}}),
        std::invalid_argument);
}

TEST_CASE("History evidence distinguishes online and smoothed estimates") {
    const History::Composition composition = {
        {{"c"}, {"e"}}, {{"h"}, {"e"}}, {{"a"}, {"f"}}, {{"a"}, {"b"}},
        {{"c"}, {"f"}}, {{"b"}, {"e"}}, {{"c"}, {"b"}}, {{"f"}, {"e"}},
        {{"b"}, {"h"}}, {{"c"}, {"e"}}, {{"f"}, {"h"}}, {{"b"}, {"h"}},
        {{"a"}, {"e"}}, {{"h"}, {"d"}}, {{"d"}, {"h"}}, {{"f"}, {"a"}}
    };
    History smoothed(composition);
    History online(composition, {}, {}, {}, 0, 3, 1, 0.15, 0, true);
    smoothed.convergence(16, 0);
    online.convergence(16, 0);

    CHECK(smoothed.geometric_mean() == doctest::Approx(0.4873978682));
    CHECK(online.geometric_mean() == doctest::Approx(0.5283514026));
}

TEST_CASE("History propagates per-event weights") {
    const History::Composition composition = {
        {{"a"}, {"b"}},
        {{"b"}, {"a"}}
    };
    const History::HistoryWeights weights = {
        {{5}, {4}},
        {{5}, {4}}
    };
    History history(
        composition, {}, {}, {}, 2, 6, 1, 0, 0, false, weights);
    history.forward_propagation();
    const auto curves = history.learning_curves();
    CHECK(curves.at("a")[0].second.is_approx(Gaussian(5.53765944, 4.758722), 1e-4));
    CHECK(curves.at("b")[0].second.is_approx(Gaussian(-0.83012755, 5.2395689), 1e-4));
    CHECK(curves.at("a")[1].second.is_approx(Gaussian(1.7922776, 4.099566689), 1e-4));
    CHECK(curves.at("b")[1].second.is_approx(Gaussian(4.8455331752, 3.7476161), 1e-4));
}

TEST_CASE("History applies one observation model to all games") {
    History history(
        {{{"a"}, {"b"}}, {{"a"}, {"b"}}},
        {{2.0, 1.0}, {3.0, 2.0}},
        {}, {}, 0, 3, 1, 0.15, 0, false, {},
        {Observation::Continuous});
    history.forward_propagation();
    CHECK(history.batch_observations[0][0] == Observation::Continuous);
    CHECK(history.batch_observations[1][0] == Observation::Continuous);

    const auto priors = history.within_priors(1, 0);
    CHECK(priors.size() == 2);
    CHECK(priors[0].size() == 1);
    CHECK_THROWS_AS(history.within_priors(3, 0), std::out_of_range);
}
