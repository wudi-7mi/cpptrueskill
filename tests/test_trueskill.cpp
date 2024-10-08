#include <doctest/doctest.h>
#include "../src/trueskill.hpp"

using namespace trueskill;
using namespace math;


TEST_CASE("vw function tests") {
    double diff = 1.0;
    double draw_margin = 0.5;

    double v_w = v_win(diff, draw_margin);
    double v_d = v_draw(diff, draw_margin);
    double w_w = w_win(diff, draw_margin);
    double w_d = w_draw(diff, draw_margin);

    CHECK(v_w == doctest::Approx(0.509160));
    CHECK(v_d == doctest::Approx(-0.920644));
    CHECK(w_w == doctest::Approx(0.513825));
    CHECK(w_d == doctest::Approx(0.923058));
}

TEST_CASE("Quality and rate 1vs1 test") {
    Rating alice;
    Rating bob;

    double quality = quality_1vs1(alice, bob);
    CHECK(quality == doctest::Approx(0.447214));

    auto result = rate_1vs1(alice, bob);
    auto new_alice = result[0][0];
    auto new_bob = result[1][0];

    double alice_mu = new_alice.mu();
    double alice_sigma = new_alice.sigma();
    double bob_mu = new_bob.mu();
    double bob_sigma = new_bob.sigma();

    CHECK(alice_mu == doctest::Approx(29.396));
    CHECK(alice_sigma == doctest::Approx(7.17148));
    CHECK(bob_mu == doctest::Approx(20.604));
    CHECK(bob_sigma == doctest::Approx(7.17148));
}