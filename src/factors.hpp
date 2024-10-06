#ifndef FACTORS_HPP
#define FACTORS_HPP

#include <vector>
#include <memory>
#include <map>
#include <cmath>
#include "mathematics.hpp"

class Factor;

class Variable {
public:
    Variable(double mu = 0.0, double sigma = 1.0) : value_(mu, sigma) {}
    void updateValue(const trueskill::math::Gaussian& newValue) { value_ = newValue; }
    trueskill::math::Gaussian getValue() const { return value_; }
    void updateMessage(const Factor* sender, const trueskill::math::Gaussian& msg) {
        messages_[sender] = msg;
    }
    trueskill::math::Gaussian getMessageFrom(const Factor* sender) const {
        auto it = messages_.find(sender);
        return it != messages_.end() ? it->second : trueskill::math::Gaussian();
    }

private:
    trueskill::math::Gaussian value_;
    std::map<const Factor*, trueskill::math::Gaussian> messages_;
};

class Factor {
public:
    virtual void updateMessage(const Variable* receiver) = 0;
    virtual ~Factor() = default;
protected:
    std::vector<std::shared_ptr<Variable>> variables_;
};

class PriorFactor : public Factor {
public:
    PriorFactor(std::shared_ptr<Variable> variable, double mu, double sigma)
        : prior_(mu, sigma) {
        variables_.push_back(variable);
    }

    void updateMessage(const Variable* receiver) override {
        receiver->updateValue(prior_);
    }

private:
    trueskill::math::Gaussian prior_;
};

class PerformanceFactor : public Factor {
public:
    PerformanceFactor(std::shared_ptr<Variable> team_perf, const std::vector<std::shared_ptr<Variable>>& player_skills)
        : beta_(4.166) {
        variables_.push_back(team_perf);
        variables_.insert(variables_.end(), player_skills.begin(), player_skills.end());
    }

    void updateMessage(const Variable* receiver) override {
        trueskill::math::Gaussian teamPerf;
        for (const auto& var : variables_) {
            if (var.get() != receiver) {
                teamPerf = teamPerf * var->getValue();
            }
        }
        trueskill::math::Gaussian msg = trueskill::math::Gaussian(teamPerf.mu(), std::sqrt(teamPerf.variance() + beta_ * beta_));
        receiver->updateMessage(this, msg);
    }

private:
    double beta_;
};

class ComparisonFactor : public Factor {
public:
    ComparisonFactor(std::shared_ptr<Variable> team1_perf, std::shared_ptr<Variable> team2_perf, bool draw)
        : draw_(draw), epsilon_(0.1) {
        variables_.push_back(team1_perf);
        variables_.push_back(team2_perf);
    }

    void updateMessage(const Variable* receiver) override {
        const Variable* other = (receiver == variables_[0].get()) ? variables_[1].get() : variables_[0].get();
        trueskill::math::Gaussian msgFromOther = other->getMessageFrom(this);
        trueskill::math::Gaussian otherValue = other->getValue();
        trueskill::math::Gaussian likelihoodTimes = otherValue / msgFromOther;

        double c = std::sqrt(1 + epsilon_ * epsilon_);
        double v = trueskill::math::V((likelihoodTimes.mu() - receiver->getValue().mu()) / c, draw_ ? epsilon_ : 0);
        double w = trueskill::math::W((likelihoodTimes.mu() - receiver->getValue().mu()) / c, draw_ ? epsilon_ : 0);
        double newMu = receiver->getValue().mu() + receiver->getValue().sigma() * receiver->getValue().sigma() * v / c;
        double newSigma = receiver->getValue().sigma() * std::sqrt(1 - receiver->getValue().sigma() * receiver->getValue().sigma() * w / (c * c));

        trueskill::math::Gaussian newMsg(newMu, newSigma);
        receiver->updateMessage(this, newMsg);
    }

private:
    bool draw_;
    double epsilon_;
};

class FactorGraph {
public:
    std::shared_ptr<Variable> createVariable(double mu = 0.0, double sigma = 1.0) {
        auto var = std::make_shared<Variable>(mu, sigma);
        variables_.push_back(var);
        return var;
    }

    void addFactor(std::shared_ptr<Factor> factor) {
        factors_.push_back(factor);
    }

    void runInference(int iterations = 5) {
        for (int i = 0; i < iterations; ++i) {
            for (const auto& factor : factors_) {
                for (const auto& var : variables_) {
                    factor->updateMessage(var.get());
                }
            }
        }
    }

    const std::vector<std::shared_ptr<Variable>>& getVariables() const {
        return variables_;
    }

private:
    std::vector<std::shared_ptr<Variable>> variables_;
    std::vector<std::shared_ptr<Factor>> factors_;
};

#endif // FACTORS_HPP