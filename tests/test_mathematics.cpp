#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "../src/mathematics.hpp"

using namespace trueskill::math;

TEST_CASE("Gaussian class tests") {
    SUBCASE("Constructor and basic getters") {
        Gaussian g(10, 2);
        CHECK(g.mu() == doctest::Approx(10));
        CHECK(g.sigma() == doctest::Approx(2));
        CHECK(g.variance() == doctest::Approx(4));
        CHECK(g.precision() == doctest::Approx(0.25));
        CHECK(g.precisionMean() == doctest::Approx(2.5));
    }

    SUBCASE("Operators") {
        Gaussian g1(10, 2);
        Gaussian g2(20, 3);

        Gaussian g3 = g1 * g2;
        CHECK(g3.mu() == doctest::Approx(13.0769));
        CHECK(g3.sigma() == doctest::Approx(1.6641));

        Gaussian g4 = g1 / g2;
        CHECK(g4.mu() == doctest::Approx(2));
        CHECK(g4.sigma() == doctest::Approx(2.68328));
    }
}

TEST_CASE("Math function tests") {
    CHECK(cdf(0) == doctest::Approx(0.5));
    CHECK(pdf(0) == doctest::Approx(0.398942));
    CHECK(ppf(0.5) == doctest::Approx(0).epsilon(0.01));
}