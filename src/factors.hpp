#ifndef FACTORS_HPP
#define FACTORS_HPP

#include <vector>
#include <memory>
#include <map>
#include <cmath>
#include <cassert>
#include <iostream>

#include "mathematics.hpp"

namespace trueskill {

constexpr double inf = std::numeric_limits<double>::infinity();

class Node {
public:
    virtual ~Node() = default;
};

class Variable : public Node, public math::Gaussian {
public:
    std::map<Node*, math::Gaussian> messages;

    Variable() : math::Gaussian() {}

    double set(const Gaussian& val) {
        double delta = this->delta(val);
        this->set_pi(val.pi());
        this->set_tau(val.tau());
        return delta;
    }

    double delta(const Gaussian& other) const {
        double pi_delta = std::abs(this->pi() - other.pi());
        if (pi_delta == inf) {
            return 0.0;
        }
        return std::max(std::abs(this->tau() - other.tau()), std::sqrt(pi_delta));
    }

    double update_message(Node* factor, double pi = 0, double tau = 0) {
        auto old_message = this->messages[factor];
        auto message = Gaussian::fromPiTau(pi, tau);
        this->messages[factor] = message;
        return this->set(*this / old_message * message);
    }

    double update_message(Node* factor, Gaussian message) {
        auto old_message = this->messages[factor];
        this->messages[factor] = message;
        return this->set(*this / old_message * message);
    }

    double update_value(Node* factor, double pi = 0, double tau = 0) {
        auto old_message = this->messages[factor];
        auto value = Gaussian::fromPiTau(pi, tau);
        this->messages[factor] = value * old_message / *this;
        return this->set(value);
    }

    double update_value(Node* factor, Gaussian value) {
        auto old_message = this->messages[factor];
        this->messages[factor] = value * old_message / *this;
        return this->set(value);
    }

    Gaussian operator[](Node* factor) const {
        return this->messages.at(factor);
    }

    Gaussian& operator[](Node* factor) {
        return this->messages[factor];
    }
};

class Factor : public Node {
public:
    std::vector<Variable*> vars;

    Factor(const std::vector<Variable*>& variables) : vars(variables) {
        for (Variable* var : variables) {
            var->messages[this] = math::Gaussian();
        }
    }

    virtual double down() {
        return 0;
    }

    virtual double up() {
        return 0;
    }

    Variable* var() const {
        assert(vars.size() == 1);
        return vars[0];
    }
};

class PriorFactor : public Factor {
public:
    math::Gaussian val;
    double dynamic;

    PriorFactor(Variable* var, const math::Gaussian& val, double dynamic = 0)
        : Factor({ var }), val(val), dynamic(dynamic) {}

    double down() override {
        double sigma = std::sqrt(std::pow(val.sigma(), 2) + std::pow(dynamic, 2));
        math::Gaussian value(val.mu(), sigma);
        return vars[0]->update_value(this, value);
    }
};

class LikelihoodFactor : public Factor {
public:
    double variance;

    LikelihoodFactor(Variable* mean_var, Variable* value_var, double variance)
        : Factor({ mean_var, value_var }), variance(variance) {}

    double calc_a(const math::Gaussian& var) const {
        return 1.0 / (1.0 + variance * var.pi());
    }

    double down() override {
        math::Gaussian msg = *vars[0] / (*vars[0])[this];
        double a = calc_a(msg);
        return vars[1]->update_message(this, a * msg.pi(), a * msg.tau());
    }

    double up() override {
        math::Gaussian msg = *vars[1] / (*vars[1])[this];
        double a = calc_a(msg);
        return vars[0]->update_message(this, a * msg.pi(), a * msg.tau());
    }
};

class SumFactor : public Factor {
public:
    std::vector<double> coeffs;

    SumFactor(Variable* sum_var, const std::vector<Variable*>& term_vars, const std::vector<double>& coeffs)
        : Factor({ sum_var }), coeffs(coeffs) {
        vars.insert(vars.end(), term_vars.begin(), term_vars.end());
    }

    double down() override {
        return update(vars[0], std::vector<Variable*>(vars.begin() + 1, vars.end()),
                      coeffs);
    }

    double up(size_t index) {
        double coeff = coeffs[index];
        std::vector<double> new_coeffs;
        for (size_t i = 0; i < coeffs.size(); ++i) {
            if (i == index) {
                new_coeffs.push_back(1.0 / coeff);
            } else {
                new_coeffs.push_back(-coeffs[i] / coeff);
            }
        }
        std::vector<Variable*> vals(vars.begin() + 1, vars.end());
        vals[index] = vars[0];
        return update(vars[index + 1], vals, new_coeffs);
    }

    double update(Variable* var, const std::vector<Variable*>& vals, const std::vector<double>& coeffs) {
        double pi_inv = 0;
        double mu = 0;
        for (size_t i = 0; i < vals.size(); ++i) {
            math::Gaussian div = *vals[i] / (*vals[i])[this];
            mu += coeffs[i] * div.mu();
            if (pi_inv == inf) continue;
            pi_inv += std::pow(coeffs[i], 2) / div.pi();
        }
        double pi = 1.0 / pi_inv;
        double tau = pi * mu;
        return var->update_message(this, pi, tau);
    }
};

class TruncateFactor : public Factor {
public:
    double(*v_func)(double, double);
    double(*w_func)(double, double);
    double draw_margin;

    TruncateFactor(Variable* var, double(*v_func)(double, double), double(*w_func)(double, double), double draw_margin)
        : Factor({ var }), v_func(v_func), w_func(w_func), draw_margin(draw_margin) {}

    double up() override {
        Variable* val = vars[0];
        math::Gaussian msg = (*val) / (*val)[this];
        double sqrt_pi = std::sqrt(msg.pi());
        double v = v_func(msg.tau() / sqrt_pi, draw_margin * sqrt_pi);
        double w = w_func(msg.tau() / sqrt_pi, draw_margin * sqrt_pi);
        double denom = 1.0 - w;
        double pi = msg.pi() / denom;
        double tau = (msg.tau() + sqrt_pi * v) / denom;
        return val->update_value(this, pi, tau);
    }
};

} // namespace trueskill

#endif // FACTORS_HPP