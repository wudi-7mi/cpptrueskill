#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN
#include <doctest/doctest.h>
#include "../src/mathematics.hpp"

using namespace trueskill::math;


TEST_CASE("Math function tests") {
    CHECK(cdf(0) == doctest::Approx(0.5));
    CHECK(pdf(0) == doctest::Approx(0.398942));
    CHECK(ppf(0.5) == doctest::Approx(0).epsilon(0.01));
}

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

TEST_CASE("Matrix class tests") {
    SUBCASE("Matrix Initialization") {
        // 测试通过二维向量初始化
        std::vector<std::vector<double>> data = {
            {1, 2, 3},
            {4, 5, 6},
            {7, 8, 9}
        };
        Matrix m1(data);
        CHECK(m1.height == 3);
        CHECK(m1.width == 3);
        CHECK(m1.data[0][0] == 1);
        CHECK(m1.data[1][1] == 5);
        CHECK(m1.data[2][2] == 9);

        // 测试通过 map 初始化
        std::map<std::pair<int, int>, double> sparse_data = {
            {{0, 0}, 1.0},
            {{1, 1}, 2.0},
            {{2, 2}, 3.0}
        };
        Matrix m2(sparse_data, 3, 3);
        CHECK(m2.data[0][0] == 1.0);
        CHECK(m2.data[1][1] == 2.0);
        CHECK(m2.data[2][2] == 3.0);
        CHECK(m2.data[0][1] == 0.0);
    }

    SUBCASE("Matrix Transpose") {
        std::vector<std::vector<double>> data = {
            {1, 2, 3},
            {4, 5, 6}
        };
        Matrix m1(data);
        Matrix t = m1.transpose();
        CHECK(t.height == 3);
        CHECK(t.width == 2);
        CHECK(t.data[0][0] == 1);
        CHECK(t.data[1][0] == 2);
        CHECK(t.data[0][1] == 4);
    }

    SUBCASE("Matrix Addition") {
        std::vector<std::vector<double>> data1 = {
            {1, 2},
            {3, 4}
        };
        std::vector<std::vector<double>> data2 = {
            {5, 6},
            {7, 8}
        };
        Matrix m1(data1);
        Matrix m2(data2);
        Matrix sum = m1 + m2;
        CHECK(sum.data[0][0] == 6);
        CHECK(sum.data[0][1] == 8);
        CHECK(sum.data[1][0] == 10);
        CHECK(sum.data[1][1] == 12);
    }

    SUBCASE("Matrix Multiplication") {
        std::vector<std::vector<double>> data1 = {
            {1, 2},
            {3, 4}
        };
        std::vector<std::vector<double>> data2 = {
            {2, 0},
            {1, 2}
        };
        Matrix m1(data1);
        Matrix m2(data2);
        Matrix product = m1 * m2;
        CHECK(product.data[0][0] == 4);
        CHECK(product.data[0][1] == 4);
        CHECK(product.data[1][0] == 10);
        CHECK(product.data[1][1] == 8);
    }

    SUBCASE("Matrix Determinant") {
        std::vector<std::vector<double>> data = {
            {1, 2},
            {3, 4}
        };
        Matrix m(data);
        CHECK(m.determinant() == -2);
        
        std::vector<std::vector<double>> data2 = {
            {6, 1, 1},
            {4, -2, 5},
            {2, 8, 7}
        };
        Matrix m2(data2);
        CHECK(m2.determinant() == -306);
    }

    SUBCASE("Matrix Inverse") {
        std::vector<std::vector<double>> data = {
            {4, 7},
            {2, 6}
        };
        Matrix m(data);
        Matrix inv = m.inverse();
        CHECK(inv.data[0][0] == doctest::Approx(0.6));
        CHECK(inv.data[0][1] == doctest::Approx(-0.7));
        CHECK(inv.data[1][0] == doctest::Approx(-0.2));
        CHECK(inv.data[1][1] == doctest::Approx(0.4));
    }

    SUBCASE("Matrix Adjugate") {
        std::vector<std::vector<double>> data = {
            {1, 2, 3},
            {0, 1, 4},
            {5, 6, 0}
        };
        Matrix m(data);
        Matrix adj = m.adjugate();
        CHECK(adj.data[0][0] == -24);
        CHECK(adj.data[0][1] == 18);
        CHECK(adj.data[0][2] == 5);
        CHECK(adj.data[1][0] == 20);
        CHECK(adj.data[1][1] == -15);
        CHECK(adj.data[1][2] == -4);
        CHECK(adj.data[2][0] == -5);
        CHECK(adj.data[2][1] == 4);
        CHECK(adj.data[2][2] == 1);
    }
}
