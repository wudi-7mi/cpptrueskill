#include <doctest/doctest.h>

#include "../src/skill.hpp"

using namespace trueskill;
using namespace trueskill::math;

TEST_CASE("Skill combines temporal and event messages") {
    Skill skill(2);
    skill.forward = Gaussian(1, 3);
    skill.backward = Gaussian(2, 4);
    const Gaussian event_likelihood(3, 5);

    const auto initial_step = skill.update_likelihood(2, event_likelihood);
    CHECK(initial_step.first == doctest::Approx(3));
    CHECK(initial_step.second == std::numeric_limits<double>::infinity());
    CHECK(skill.posterior().is_approx(
        skill.forward * skill.backward * event_likelihood));
    CHECK(skill.forward_posterior().is_approx(skill.forward * event_likelihood));
    CHECK(skill.backward_posterior().is_approx(skill.backward * event_likelihood));
}

TEST_CASE("Skill tracks multiple events in one time batch") {
    Skill skill(0);
    skill.add_event(3);
    skill.update_likelihood(3, Gaussian(-2, 6));

    CHECK(skill.events.size() == 2);
    CHECK(skill.likelihood(3).mu() == doctest::Approx(-2));
    CHECK_THROWS_AS(skill.likelihood(1), std::out_of_range);
}
