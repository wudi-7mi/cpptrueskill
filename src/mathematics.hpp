#ifndef MATHEMATICS_HPP
#define MATHEMATICS_HPP

#include <cmath>
#include <limits>

namespace trueskill {
namespace math {

const double PI = 3.14159265358979323846;
const double SQRT_2PI = std::sqrt(2 * PI);
const double SQRT_PI = std::sqrt(PI);

inline double erf_approx(double x) {
    const double a1 = 0.278393;
    const double a2 = 0.230389;
    const double a3 = 0.000972;
    const double a4 = 0.078108;
    double abs_x = std::abs(x);
    double t = 1.0 / (1.0 + a1 * abs_x + a2 * abs_x * abs_x + a3 * abs_x * abs_x * abs_x + a4 * abs_x * abs_x * abs_x * abs_x);
    double erf_val = 1.0 - t * t * t * t;
    return x >= 0 ? erf_val : -erf_val;
}

inline double erf_inv(double x) {
    const double a = 0.147;  // 常数项
    double sign_x = x >= 0 ? 1 : -1;
    double ln_term = std::log(1 - x * x);
    double part1 = 2 / (PI * a) + ln_term / 2;
    double part2 = std::pow(part1, 2) - ln_term / a;
    return sign_x * std::sqrt(std::sqrt(part2) - part1);
}

// 正态分布的累积分布函数 (CDF)
inline double cdf(double x) {
    return 0.5 * (1 + erf_approx(x / std::sqrt(2)));
}

// 正态分布的概率密度函数 (PDF)
inline double pdf(double x) {
    return std::exp(-std::pow(x, 2) / 2) / SQRT_2PI;
}

// 正态分布的逆累积分布函数
inline double ppf(double x) {
    return erf_inv(2 * x - 1) * std::sqrt(2);
}

// 计算两个高斯分布的KL散度
inline double kl_divergence(double mu1, double sigma1, double mu2, double sigma2) {
    return std::log(sigma2 / sigma1) + (std::pow(sigma1, 2) + std::pow(mu1 - mu2, 2)) / (2 * std::pow(sigma2, 2)) - 0.5;
}

// 计算两个高斯分布的Hellinger距离
inline double hellinger_distance(double mu1, double sigma1, double mu2, double sigma2) {
    double term1 = std::pow(mu1 - mu2, 2) / (4 * (std::pow(sigma1, 2) + std::pow(sigma2, 2)));
    double term2 = 0.5 * std::log((std::pow(sigma1, 2) + std::pow(sigma2, 2)) / (2 * sigma1 * sigma2));
    return std::sqrt(1 - std::exp(-term1 - term2));
}

// 计算两个高斯分布的Wasserstein距离
inline double wasserstein_distance(double mu1, double sigma1, double mu2, double sigma2) {
    return std::sqrt(std::pow(mu1 - mu2, 2) + std::pow(sigma1 - sigma2, 2));
}

class Gaussian {
public:

    Gaussian(double mu, double sigma) {
        if (sigma <= 0) {
            throw std::invalid_argument("sigma should be greater than 0");
        }

        double pi = pow(sigma, -2);
        double tau = mu * pi;
        
        this->_mu = mu;
        this->_sigma = sigma;
        this->_pi = pi;
        this->_tau = tau;
    }

    double mu() const { return _mu; }
    double sigma() const { return _sigma; }
    double variance() const { return _sigma * _sigma; }
    double precision() const { return 1.0 / variance(); }
    double precisionMean() const { return precision() * _mu; }

    static Gaussian fromPrecisionMean(double precisionMean, double precision) {
        double sigma = std::sqrt(1.0 / precision);
        return Gaussian(precisionMean / precision, sigma);
    }

    Gaussian operator*(const Gaussian& other) const {
        double newPrecision = precision() + other.precision();
        double newPrecisionMean = precisionMean() + other.precisionMean();
        return fromPrecisionMean(newPrecisionMean, newPrecision);
    }

    Gaussian operator/(const Gaussian& other) const {
        double newPrecision = precision() - other.precision();
        double newPrecisionMean = precisionMean() - other.precisionMean();
        return fromPrecisionMean(newPrecisionMean, newPrecision);
    }

    double logProductNormalization(const Gaussian& other) const {
        if (_sigma == 0 || other._sigma == 0) {
            return std::numeric_limits<double>::infinity();
        }
        double varianceSum = variance() + other.variance();
        double muDiff = _mu - other._mu;
        double logSqrt2Pi = std::log(std::sqrt(2 * PI));
        return -logSqrt2Pi - std::log(std::sqrt(varianceSum)) - 
               std::pow(muDiff, 2) / (2 * varianceSum);
    }

    bool operator==(const Gaussian& other) const {
        return (_mu == other._mu) && (_sigma == other._sigma);
    }

    bool operator<(const Gaussian& other) const {
        return _mu < other._mu;
    }

    bool operator<=(const Gaussian& other) const {
        return _mu <= other._mu;
    }

    bool operator>(const Gaussian& other) const {
        return _mu > other._mu;
    }

    bool operator>=(const Gaussian& other) const {
        return _mu >= other._mu;
    }

    // Print function
    friend std::ostream& operator<<(std::ostream& os, const Gaussian& g) {
        os << "Gaussian(mu=" << g.mu() << ", sigma=" << g.sigma() << ")";
        return os;
    }

private:
    double _mu = 25;         // Mean
    double _sigma = 8.333;   // Standard deviation
    double _pi;
    double _tau;
};

} // namespace math
} // namespace trueskill

#endif // MATHEMATICS_HPP