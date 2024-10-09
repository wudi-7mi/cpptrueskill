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

    SUBCASE("1v1 win") {
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
    
    SUBCASE("1v1 draw") {
        Rating alice(21, 8);
        Rating bob;

        double quality = quality_1vs1(alice, bob);
        CHECK(quality == doctest::Approx(0.433285));

        auto result = rate_1vs1(alice, bob, true);
        auto new_alice = result[0][0];
        auto new_bob = result[1][0];

        double alice_mu = new_alice.mu();
        double alice_sigma = new_alice.sigma();
        double bob_mu = new_bob.mu();
        double bob_sigma = new_bob.sigma();

        CHECK(alice_mu == doctest::Approx(22.5207));
        CHECK(alice_sigma == doctest::Approx(6.29868));
        CHECK(bob_mu == doctest::Approx(23.350));
        CHECK(bob_sigma == doctest::Approx(6.38765));
    }

}

TEST_CASE("Quality and rate multiple teams test") {

    SUBCASE("Calculate ratings for different teams") {
        trueskill::TrueSkill ts(25.0, 8.333333, 4.2, 0.3, 0.1);

        Rating alice(34.23, 7.22);
        Rating bob(24.23, 3.22);   
        Rating charlie(12.23, 3.9);
        Rating david(26.23, 6.11);
        Rating eve(41.2, 4.22);
        Rating frank(24.23, 5.08);
        Rating gloria(4.23, 3.08);

        std::vector<std::vector<double>> weights = {{1.0, 1.0}, {0.8, 1.0}, {1.0, 1.0, 1.0}};

        std::vector<std::vector<Rating>> teams = {{alice, bob}, {charlie, david}, {eve, frank, gloria}};

        double quality = ts.quality(teams, weights);
        CHECK(quality == doctest::Approx(0.01439));

        auto newTeams = ts.rate(teams, {2, 1, 3}, weights);
        auto new_alice = newTeams[0][0];
        auto new_bob = newTeams[0][1];
        auto new_charlie = newTeams[1][0];
        auto new_david = newTeams[1][1];
        auto new_eve = newTeams[2][0];
        auto new_frank = newTeams[2][1];
        auto new_gloria = newTeams[2][2];

        CHECK(new_alice.mu() == doctest::Approx(31.0297));
        CHECK(new_bob.mu() == doctest::Approx(23.589));
        CHECK(new_charlie.mu() == doctest::Approx(15.5848));
        CHECK(new_david.mu() == doctest::Approx(36.4867));
        CHECK(new_eve.mu() == doctest::Approx(37.3913));
        CHECK(new_frank.mu() == doctest::Approx(18.7194));
        CHECK(new_gloria.mu() == doctest::Approx(2.1922));
        
        CHECK(new_alice.sigma() == doctest::Approx(5.8226));
        CHECK(new_bob.sigma() == doctest::Approx(3.11829));
        CHECK(new_charlie.sigma() == doctest::Approx(3.77367));
        CHECK(new_david.sigma() == doctest::Approx(5.24585));
        CHECK(new_eve.sigma() == doctest::Approx(4.00044));
        CHECK(new_frank.sigma() == doctest::Approx(4.68291));
        CHECK(new_gloria.sigma() == doctest::Approx(3.00565));
    }
}