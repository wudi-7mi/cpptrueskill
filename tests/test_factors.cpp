#include <doctest/doctest.h>
#include "../src/factors.hpp"

using namespace trueskill;
using namespace math;

// 创建一个辅助函数用于比较浮点数
bool double_compare(double a, double b, double epsilon = 1e-6) {
    return std::fabs(a - b) < epsilon;
}

// Test for Variable class
TEST_CASE("Variable class set and delta") {
    Variable var;
    Gaussian g1(2.0, 3.0);  // 创建一个高斯分布
    var.set(g1);  // 设置高斯分布值

    Gaussian g2(2.5, 3.5);  // 创建另一个高斯分布
    double delta_value = var.delta(g2);

    CHECK(double_compare(delta_value, 0.5));  // 假设差值应为 0.5
}

TEST_CASE("Variable class update_message") {
    Variable var;
    Node factor;
    Gaussian message(1.0, 2.0);
    var.update_message(&factor, message._pi, message._tau, message);

    CHECK(double_compare(var.messages[&factor]._pi, 1.0));
    CHECK(double_compare(var.messages[&factor]._tau, 2.0));
}

// Test for PriorFactor class
TEST_CASE("PriorFactor down") {
    Variable var;
    Gaussian val(1.0, 2.0);
    PriorFactor factor(&var, val, 0.5);  // 添加一个动态因素
    double down_result = factor.down();

    CHECK(double_compare(down_result, 0.0));  // 假设输出结果
    CHECK(double_compare(var.messages[&factor]._pi, val._pi));
    CHECK(double_compare(var.messages[&factor]._tau, val._tau));
}

// Test for LikelihoodFactor class
TEST_CASE("LikelihoodFactor down") {
    Variable mean_var, value_var;
    LikelihoodFactor factor(&mean_var, &value_var, 0.1);

    // 先更新 mean_var 的消息，再检查 value_var 的更新
    factor.down();
    CHECK(double_compare(value_var.messages[&factor]._pi, 0.9090909));  // 假设预期输出
    CHECK(double_compare(value_var.messages[&factor]._tau, 0.9090909));  // 假设预期输出
}

TEST_CASE("LikelihoodFactor up") {
    Variable mean_var, value_var;
    LikelihoodFactor factor(&mean_var, &value_var, 0.1);

    // 先更新 value_var 的消息，再检查 mean_var 的更新
    factor.up();
    CHECK(double_compare(mean_var.messages[&factor]._pi, 0.9090909));  // 假设预期输出
    CHECK(double_compare(mean_var.messages[&factor]._tau, 0.9090909));  // 假设预期输出
}

// Test for SumFactor class
TEST_CASE("SumFactor down") {
    Variable sum_var;
    Variable term_var1, term_var2;
    std::vector<Variable*> terms = { &term_var1, &term_var2 };
    std::vector<double> coeffs = { 0.5, 0.5 };
    SumFactor sum_factor(&sum_var, terms, coeffs);

    sum_factor.down();

    CHECK(double_compare(sum_var.messages[&sum_factor]._pi, 1.0));  // 假设预期输出
    CHECK(double_compare(sum_var.messages[&sum_factor]._tau, 1.0));  // 假设预期输出
}

// Test for TruncateFactor class
double v_func(double x, double margin) {
    return std::exp(-x * margin);  // 示例 v 函数
}

double w_func(double x, double margin) {
    return std::exp(-x * margin);  // 示例 w 函数
}

TEST_CASE("TruncateFactor up") {
    Variable var;
    TruncateFactor truncate_factor(&var, v_func, w_func, 0.5);

    truncate_factor.up();

    CHECK(double_compare(var.messages[&truncate_factor]._pi, 1.0));  // 假设预期输出
    CHECK(double_compare(var.messages[&truncate_factor]._tau, 1.0));  // 假设预期输出
}
