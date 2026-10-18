#include <doctest/doctest.h>

#include "../src/match.hpp"

using namespace trueskill;
using namespace trueskill::math;

TEST_CASE("Player defaults and performance") {
    const Player player;

    CHECK(player.prior.mu() == doctest::Approx(0));
    CHECK(player.prior.sigma() == doctest::Approx(6));
    CHECK(player.beta == doctest::Approx(1));
    CHECK(player.gamma == doctest::Approx(0.03));

    const Gaussian performance = player.performance();
    CHECK(performance.mu() == doctest::Approx(0));
    CHECK(performance.sigma() == doctest::Approx(std::sqrt(37.0)));
}

TEST_CASE("Player supports custom and deterministic performance") {
    const Player custom(Gaussian(2, 2), 1, 0.15);
    const Gaussian noisy_performance = custom.performance();

    CHECK(noisy_performance.mu() == doctest::Approx(2));
    CHECK(noisy_performance.sigma() == doctest::Approx(std::sqrt(5.0)));
    CHECK(custom.gamma == doctest::Approx(0.15));

    const Player deterministic(Gaussian(-3, 4), 0, 0);
    const Gaussian deterministic_performance = deterministic.performance();
    CHECK(deterministic_performance.mu() == doctest::Approx(-3));
    CHECK(deterministic_performance.sigma() == doctest::Approx(4));
}

TEST_CASE("Player rejects negative dynamics parameters") {
    CHECK_THROWS_AS(Player(Gaussian(), -1, 0), std::invalid_argument);
    CHECK_THROWS_AS(Player(Gaussian(), 1, -0.1), std::invalid_argument);
}

TEST_CASE("MatchPlayer preserves identified rating holder") {
    MatchPlayer player(42, Rating(25, 8));
    CHECK(player.id == 42);
    CHECK(player.getRating().mu() == doctest::Approx(25));

    player.setRating(Rating(30, 6));
    CHECK(player.getRating().mu() == doctest::Approx(30));
    CHECK(player.getRating().sigma() == doctest::Approx(6));
}
