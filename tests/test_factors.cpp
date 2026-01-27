#include <doctest/doctest.h>
#include "../src/factors.hpp"

using namespace trueskill;
using namespace math;

// Test for Variable class
TEST_CASE("Variable class set and delta") {
    Variable var;
    Gaussian g1(2.0, 3.0);  // Create a Gaussian distribution
    var.set(g1);  // Set the Gaussian distribution value

    Gaussian g2(2.5, 3.5);  // Create another Gaussian distribution
    double delta_value = var.delta(g2);

    CHECK(delta_value == doctest::Approx(0.171693));  // Assume the difference should be 0.5
}

// Test for PriorFactor class
TEST_CASE("PriorFactor down") {
    Variable var;
    Gaussian val(1.0, 2.0);
    PriorFactor factor(&var, val, 0.5);  // Add a dynamic factor
    double down_result = factor.down();

    CHECK(down_result == doctest::Approx(0.485071));
    CHECK(var.messages[&factor].pi() == doctest::Approx(0.235294));
    CHECK(var.messages[&factor].tau() == doctest::Approx(0.235294));
}