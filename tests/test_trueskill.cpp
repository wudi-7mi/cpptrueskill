#include <doctest/doctest.h>
#include "../src/trueskill.hpp"

using namespace trueskill;
using namespace math;


TEST_CASE("Trueskill function tests") {
    TrueSkill ts;
    double diff = 1.0;
    double draw_margin = 0.5;

    double v_win = ts.v_win(diff, draw_margin);
    double v_draw = ts.v_draw(diff, draw_margin);
    double w_win = ts.w_win(diff, draw_margin);
    double w_draw = ts.w_draw(diff, draw_margin);

    CHECK(v_win == doctest::Approx(0.509160));
    CHECK(v_draw == doctest::Approx(-0.920644));
    CHECK(w_win == doctest::Approx(0.513825));
    CHECK(w_draw == doctest::Approx(0.923058));
}
