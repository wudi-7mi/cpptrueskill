#ifndef MATHEMATICS_HPP
#define MATHEMATICS_HPP

#include <cmath>
#include <limits>
#include <vector>
#include <map>
#include <iomanip>


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

class Matrix {
public:
    std::vector<std::vector<double>> data;
    int height, width;

    Matrix(const std::vector<std::vector<double>>& src) : height(src.size()), width(src[0].size()), data(src) {
        for (const auto& row : src) {
            if (row.size() != width) {
                throw std::invalid_argument("src must be a rectangular array of numbers");
            }
        }
    }

    Matrix(const std::map<std::pair<int, int>, double>& src, int height, int width) : height(height), width(width) {
        data.resize(height, std::vector<double>(width, 0));
        for (const auto& item : src) {
            int r = item.first.first, c = item.first.second;
            data[r][c] = item.second;
        }
    }

    Matrix transpose() const {
        std::vector<std::vector<double>> transposed(width, std::vector<double>(height));
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                transposed[c][r] = data[r][c];
            }
        }
        return Matrix(transposed);
    }

    Matrix minor(int row_n, int col_n) const {
        if (row_n < 0 || row_n >= height || col_n < 0 || col_n >= width) {
            throw std::invalid_argument("Invalid row or column number");
        }
        std::vector<std::vector<double>> minorMatrix;
        for (int r = 0; r < height; r++) {
            if (r == row_n) continue;
            std::vector<double> row;
            for (int c = 0; c < width; c++) {
                if (c == col_n) continue;
                row.push_back(data[r][c]);
            }
            minorMatrix.push_back(row);
        }
        return Matrix(minorMatrix);
    }
    
    Matrix adjugate() const {
        if (height != width) {
            throw std::invalid_argument("Only square matrix can be adjugated");
        }

        std::vector<std::vector<double>> cofactors(height, std::vector<double>(width));
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                double sign = ((r + c) % 2 == 0) ? 1.0 : -1.0;
                cofactors[r][c] = sign * minor(r, c).determinant();
            }
        }
        return Matrix(cofactors).transpose();
    }

    Matrix inverse() const {
        double det = determinant();
        if (det == 0) {
            throw std::invalid_argument("Matrix is singular, can't find its inverse");
        }
        return adjugate() * (1.0 / det);
    }
    
    Matrix operator+(const Matrix& other) const {
        if (height != other.height || width != other.width) {
            throw std::invalid_argument("Matrices must be the same size for addition");
        }
        std::vector<std::vector<double>> result(height, std::vector<double>(width));
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                result[r][c] = data[r][c] + other.data[r][c];
            }
        }
        return Matrix(result);
    }

    Matrix operator*(const Matrix& other) const {
        if (width != other.height) {
            throw std::invalid_argument("Bad size for multiplication");
        }
        std::vector<std::vector<double>> result(height, std::vector<double>(other.width));
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < other.width; c++) {
                for (int k = 0; k < width; k++) {
                    result[r][c] += data[r][k] * other.data[k][c];
                }
            }
        }
        return Matrix(result);
    }

    Matrix operator*(double scalar) const {
        std::vector<std::vector<double>> result(height, std::vector<double>(width));
        for (int r = 0; r < height; r++) {
            for (int c = 0; c < width; c++) {
                result[r][c] = data[r][c] * scalar;
            }
        }
        return Matrix(result);
    }

    double determinant() const {
        if (height != width) {
            throw std::invalid_argument("Only square matrix can calculate a determinant");
        }
        if (height == 1) return data[0][0];

        if (height == 2) {
            return data[0][0] * data[1][1] - data[0][1] * data[1][0];
        }

        double det = 0;
        for (int c = 0; c < width; c++) {
            det += ((c % 2 == 0) ? 1 : -1) * data[0][c] * minor(0, c).determinant();
        }
        return det;
    }

    friend std::ostream& operator<<(std::ostream& os, const Matrix& matrix) {
        for (const auto& row : matrix.data) {
            for (double val : row) {
                os << std::setw(10) << val << " ";
            }
            os << std::endl;
        }
        return os;
    }
};

} // namespace math
} // namespace trueskill

#endif // MATHEMATICS_HPP