#include <doctest/doctest.h>

#include "../src/game.hpp"

using namespace trueskill;
using namespace trueskill::math;

TEST_CASE("Ordinal Game computes 1v1 likelihoods and posteriors") {
    const double sigma = 25.0 / 3.0;
    const double beta = 25.0 / 6.0;
    const Game game(
        {{Player(Gaussian(25, sigma), beta, 0)},
         {Player(Gaussian(25, sigma), beta, 0)}},
        {0, 1});

    const auto posteriors = game.posteriors();
    CHECK(posteriors[0][0].mu() == doctest::Approx(20.7947792561));
    CHECK(posteriors[0][0].sigma() == doctest::Approx(7.1944814226));
    CHECK(posteriors[1][0].mu() == doctest::Approx(29.2052207439));
    CHECK(posteriors[1][0].sigma() == doctest::Approx(7.1944814226));
    CHECK(game.evidence == doctest::Approx(0.5));

    CHECK((game.teams[0][0].prior * game.likelihoods[0][0]).mu()
        == doctest::Approx(posteriors[0][0].mu()));
}

TEST_CASE("Ordinal Game supports draws") {
    const double sigma = 25.0 / 3.0;
    const double beta = 25.0 / 6.0;
    const Game game(
        {{Player(Gaussian(25, sigma), beta, 0)},
         {Player(Gaussian(25, sigma), beta, 0)}},
        {0, 0},
        0.25);

    const auto posteriors = game.posteriors();
    CHECK(posteriors[0][0].mu() == doctest::Approx(25));
    CHECK(posteriors[0][0].sigma() == doctest::Approx(6.469).epsilon(0.001));
    CHECK(posteriors[1][0].mu() == doctest::Approx(25));
    CHECK(posteriors[1][0].sigma() == doctest::Approx(6.469).epsilon(0.001));
    CHECK(game.evidence == doctest::Approx(0.113315).epsilon(1e-5));

    const Game known_skills(
        {{Player(Gaussian(25, 1e-7), beta, 0)},
         {Player(Gaussian(25, 1e-7), beta, 0)}},
        {0, 0},
        0.25);
    CHECK(known_skills.evidence == doctest::Approx(0.25));
}

TEST_CASE("Ordinal Game supports multiple teams") {
    const double sigma = 25.0 / 3.0;
    const double beta = 25.0 / 6.0;
    const Player player(Gaussian(25, sigma), beta, 0);
    const Game game({{player}, {player}, {player}}, {1, 2, 0});

    const auto posteriors = game.posteriors();
    CHECK(posteriors[0][0].mu() == doctest::Approx(25.0).epsilon(1e-5));
    CHECK(posteriors[0][0].sigma() == doctest::Approx(6.238469796).epsilon(1e-5));
    CHECK(posteriors[1][0].mu() == doctest::Approx(31.3113582213).epsilon(1e-5));
    CHECK(posteriors[1][0].sigma() == doctest::Approx(6.69881865).epsilon(1e-5));
    CHECK(posteriors[2][0].mu() == doctest::Approx(18.6886417787).epsilon(1e-5));
}

TEST_CASE("Ordinal Game leaves zero-weight players unchanged") {
    const Game game(
        {{Player(Gaussian(2, 6), 1, 0)},
         {Player(Gaussian(2, 6), 1, 0)}},
        {},
        0,
        {{1}, {0}});

    const auto posteriors = game.posteriors();
    CHECK(posteriors[0][0].mu() == doctest::Approx(5.557176746).epsilon(1e-4));
    CHECK(posteriors[0][0].sigma() == doctest::Approx(4.0527906913).epsilon(1e-4));
    CHECK(posteriors[1][0].mu() == doctest::Approx(2));
    CHECK(posteriors[1][0].sigma() == doctest::Approx(6));
    CHECK(game.likelihoods[1][0].sigma() == std::numeric_limits<double>::infinity());
}

TEST_CASE("Ordinal Game supports negative weights") {
    const Game game(
        {{Player(Gaussian(2, 6), 1, 0)},
         {Player(Gaussian(2, 6), 1, 0)}},
        {},
        0,
        {{1}, {-1}});

    const auto posteriors = game.posteriors();
    CHECK(posteriors[0][0].mu() == doctest::Approx(posteriors[1][0].mu()));
    CHECK(posteriors[0][0].sigma() == doctest::Approx(posteriors[1][0].sigma()));
}

TEST_CASE("Ordinal Game returns uninformative likelihood for an expected result") {
    const Game game(
        {{Player(Gaussian(1.139, 0.531), 1, 0.2125)},
         {Player(Gaussian(15.568, 0.51), 1, 0.2125)}},
        {0, 1});

    CHECK(game.likelihoods[0][0].sigma() == std::numeric_limits<double>::infinity());
    CHECK(game.likelihoods[1][0].sigma() == std::numeric_limits<double>::infinity());
}

TEST_CASE("Ordinal Game validates its shape") {
    CHECK_THROWS_AS(Game(Game::Teams{}), std::invalid_argument);
    CHECK_THROWS_AS(Game(Game::Teams{{Player()}, {Player()}}, {1}), std::invalid_argument);
    CHECK_THROWS_AS(
        Game(Game::Teams{{Player()}, {Player()}}, {}, 0, {{1}, {}}),
        std::invalid_argument);
    CHECK_THROWS_AS(Game(Game::Teams{{Player()}, {Player()}}, {}, 1.1), std::invalid_argument);
}

TEST_CASE("Continuous Game observes performance differences") {
    const Game game(
        {{Player(Gaussian(2, 2), 1, 0)},
         {Player(Gaussian(1, 2), 1, 0)}},
        {45.24, 44.24},
        0,
        {},
        Observation::Continuous);

    const auto posteriors = game.posteriors();
    CHECK(posteriors[0][0].is_approx(Gaussian(2, 1.549193), 1e-5));
    CHECK(posteriors[1][0].is_approx(Gaussian(1, 1.549193), 1e-5));
    CHECK(game.evidence == doctest::Approx(
        pdf(1, 1, std::sqrt(10.0))));

    const Game weighted(
        {{Player(Gaussian(2, 2), 1, 0)},
         {Player(Gaussian(1, 2), 1, 0)}},
        {45.24, 44.24},
        0,
        {{1}, {2}},
        Observation::Continuous);
    const auto weighted_posteriors = weighted.posteriors();
    CHECK(weighted_posteriors[0][0].is_approx(Gaussian(2.16, 1.833030), 1e-5));
    CHECK(weighted_posteriors[1][0].is_approx(Gaussian(0.68, 1.2), 1e-5));
}

TEST_CASE("Continuous Game supports multiple teams") {
    const Game game(
        {{Player(Gaussian(2, 2), 1, 0), Player(Gaussian(2, 2), 1, 0)},
         {Player(Gaussian(3, 2), 1, 0)},
         {Player(Gaussian(4, 2), 1, 0)}},
        {4.2, 0.2, 2.1},
        0,
        {},
        Observation::Continuous);
    const auto posteriors = game.posteriors();
    CHECK(posteriors[0][0].is_approx(Gaussian(2.816, 1.649242), 1e-5));
    CHECK(posteriors[0][1].is_approx(Gaussian(2.816, 1.649242), 1e-5));
    CHECK(posteriors[1][0].is_approx(Gaussian(2.232, 1.442221), 1e-5));
    CHECK(posteriors[2][0].is_approx(Gaussian(3.952, 1.442221), 1e-5));
}

TEST_CASE("Discrete Game uses Poisson fixed-point approximation") {
    CHECK(Game::fixed_point_approx(2, 2, 2).is_approx(
        Gaussian(0.6551929425, 0.6218258872), 1e-6));

    const Game game(
        {{Player(Gaussian(4, 6), 1, 0)},
         {Player(Gaussian(0, 6), 1, 0)}},
        {0, 54},
        0,
        {{1}, {1}},
        Observation::Discrete);
    const auto posteriors = game.posteriors();
    CHECK(posteriors[0][0].is_approx(Gaussian(0.118952, 4.300102), 1e-5));
    CHECK(posteriors[1][0].is_approx(Gaussian(3.881048, 4.300102), 1e-5));
    CHECK(game.evidence > 0);

    const Game calibrated(
        {{Player(Gaussian(4, 1), 1, 0)},
         {Player(Gaussian(0, 1), 1, 0)}},
        {std::exp(4.0), 0},
        0,
        {{1}, {1}},
        Observation::Discrete);
    const auto calibrated_posteriors = calibrated.posteriors();
    CHECK(calibrated_posteriors[0][0].is_approx(Gaussian(3.997732, 0.866683), 1e-5));
    CHECK(calibrated_posteriors[1][0].is_approx(Gaussian(0.002268, 0.866683), 1e-5));
}

TEST_CASE("Ordinal Game computes many-team evidence") {
    Game::Teams teams;
    std::vector<double> result;
    for (int i = 0; i < 12; ++i) {
        teams.push_back({Player(Gaussian(0, 0), 1, 0)});
        result.push_back(i);
    }
    const Game game(teams, result);
    CHECK(game.evidence == doctest::Approx(1.2180818904e-8));
}
