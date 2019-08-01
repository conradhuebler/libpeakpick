/*
 * <Math containing Header file.>
 * Copyright (C) 2017 - 2019 Conrad Hübler <Conrad.Huebler@gmx.net>
 * 
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 * 
 */

#pragma once

#include <Eigen/Dense>

#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

typedef Eigen::VectorXd Vector;
static double pi = 3.14159265;

namespace PeakPick {

struct LinearRegression {
    Vector y, x, y_head;
    double m = 0;
    double n = 0;
    double R = 0;
    double sum_err = 0;
};

inline double SimpsonIntegrate(double lower, double upper, std::function<double(double, const Vector&)> function, const Vector& parameter, double prec = 1e-3)
{
    double integ = 0;
    double delta = prec;

    int increments = (upper - lower) / prec;
    std::cout << "prepare for numerical integration " << integ << parameter.transpose() << std::endl;
#pragma omp parallel for reduction(+ \
                                   : integ)
    for (int i = 0; i < increments; ++i) {
        double x = lower + i / double(increments);
        double b = x + delta;
        integ += (b - x) / 6 * (function(x, parameter) + 4 * function((x + b) / 2, parameter) + function(b, parameter));
    }
    std::cout << "result from integration over part of function " << integ << parameter.transpose() << std::endl;

    return integ;
}

inline double mean(const Vector& vector, unsigned int* min = NULL, unsigned int* max = NULL)
{
    if (vector.size() == 0) {
        std::cout << "PeakPick::mean() vector is empty" << std::endl;
        return 0;
    }
    double sum = 0;
    double temp_min = vector(0);
    unsigned int pos_min = 0;
    double temp_max = vector(0);
    unsigned int pos_max = 0;
    // #pragma omp parallel for shared(sum, temp_min, pos_min, temp_max, pos_max) reduction(+:sum, &:temp_min, &:pos_min, &:temp_max, &:pos_max)
    for (unsigned int i = 0; i < vector.size(); ++i) {
        sum += vector(i);
        if (vector(i) > temp_max) {
            temp_max = vector(i);
            pos_max = i;
        }

        if (vector(i) < temp_min) {
            temp_min = vector(i);
            pos_min = i;
        }

        if (max != NULL)
            *max = pos_max;

        if (min != NULL)
            *min = pos_min;
    }
    return sum / double(vector.size());
}

inline double meanThreshold(const Vector& vector, double threshold, unsigned int* min = NULL, unsigned int* max = NULL)
{
    if (vector.size() == 0)
        return 0;

    double sum = 0;
    double temp_min = vector(0);
    unsigned int pos_min = 0;
    double temp_max = vector(0);
    unsigned int pos_max = 0;
    for (int i = 0; i < vector.size(); ++i) {
        if (std::abs(vector(i)) >= threshold)
            continue;
        sum += vector(i);
        if (vector(i) > temp_max) {
            temp_max = vector(i);
            pos_max = i;
        }
        if (vector(i) < temp_min) {
            temp_min = vector(i);
            pos_min = i;
        }

        if (max != NULL)
            *max = pos_max;

        if (min != NULL)
            *min = pos_min;
    }

    return sum / double(vector.size());
}

inline double stddev(const Vector& vector, double mean)
{
    if (vector.size() == 0)
        return 0;

    double sum = 0;
#pragma omp parallel for reduction(+ \
                                   : sum)
    for (unsigned int i = 0; i < vector.size(); ++i) {
        sum += (vector(i) - mean) * (vector(i) - mean);
    }
    return sqrt(sum / double(vector.size()));
}

inline double stddevThreshold(const Vector& vector, double mean, double threshold)
{
    if (vector.size() == 0)
        return 0;

    double sum = 0;

#pragma omp parallel for reduction(+ \
                                   : sum)
    for (unsigned int i = 0; i < vector.size(); ++i) {
        if (std::abs(vector(i)) >= threshold)
            continue;
        sum += (vector(i) - mean) * (vector(i) - mean);
    }
    return sqrt(sum / double(vector.size()));
}

inline double Gaussian(double x, double a, double x_0, double c)
{
    return a * exp(-pow((x - x_0), 2) / (2 * c * c));
}

inline double Lorentzian(double x, double x_0, double gamma)
{
    return 1 / pi * (0.5 * gamma) / (pow(x - x_0, 2) + pow(0.5 * gamma, 2));
}

inline double Signal(double x, const Vector& parameter)
{
    double signal = 0;
    for (unsigned int i = 0; i < parameter.size() / 6; ++i) {
        double gaussian = Gaussian(x, parameter(1 + i * 6), parameter(0 + i * 6), parameter(2 + i * 6));
        double lorentzian = Lorentzian(x, parameter(0 + i * 6), parameter(3 + i * 6));
        signal += ((1 - parameter(5 + i * 6)) * gaussian + parameter(5 + i * 6) * lorentzian) * parameter(4 + i * 6);
    }
    return signal;
}

inline double SignalSingle(double x, const Vector& parameter, int function)
{
    double signal = 0;
    if (function >= parameter.size() / 6)
        return 0;

    double gaussian = Gaussian(x, parameter(1 + function * 6), parameter(0 + function * 6), parameter(2 + function * 6));
    double lorentzian = Lorentzian(x, parameter(0 + function * 6), parameter(3 + function * 6));
    signal += ((1 - parameter(5 + function * 6)) * gaussian + parameter(5 + function * 6) * lorentzian) * parameter(4 + function * 6);

    return signal;
}

inline double IntegrateGLFunction(const Vector& parameter)
{
    double integ = 0;
    for (unsigned int i = 0; i < parameter.size() / 6; ++i) {
        integ += (((1 - parameter(5 + i * 6)) * parameter(1 + i * 6) * parameter(2 + i * 6)) + parameter(5 + i * 6)) * parameter(4 + i * 6);
    }
    std::cout << "result from integration over whole function " << integ << parameter.transpose() << std::endl;
    return integ;
}

inline double IntegrateGLSignal(const Vector& parameter, double start, double end)
{
    std::function<double(double, const Vector&)> function = Signal;
    return SimpsonIntegrate(start, end, function, parameter);
}

/*! \brief Evaluates a polynomial function at point x with a vector of coeff
 * y = x^0*coeff[0] + x^1*coeff[1] + x^2*coeff[2] ...
 */
inline double Polynomial(double x, const Vector& coeff)
{
    double y = 0;
    for (unsigned int i = 0; i < coeff.size(); ++i)
        y += pow(x, i) * coeff(i);

    return y;
}

inline LinearRegression LeastSquares(const Vector& x, const Vector& y)
{
    LinearRegression regression;

    if (x.size() != y.size())
        return regression;
    // http://www.bragitoff.com/2015/09/c-program-to-linear-fit-the-data-using-least-squares-method/ //
    double xsum = 0, x2sum = 0, ysum = 0, xysum = 0; //variables for sums/sigma of xi,yi,xi^2,xiyi etc
    int n = x.size();
    for (int i = 0; i < n; ++i) {
        xsum += x[i]; //calculate sigma(xi)
        ysum += y[i]; //calculate sigma(yi)
        x2sum += (x[i] * x[i]); //calculate sigma(x^2i)
        xysum += x[i] * y[i]; //calculate sigma(xi*yi)
    }
    regression.m = (n * xysum - xsum * ysum) / (n * x2sum - xsum * xsum); //calculate slope
    regression.n = (x2sum * ysum - xsum * xysum) / (x2sum * n - xsum * xsum); //calculate intercept
    double mean_x = xsum / double(n);
    double mean_y = ysum / double(n);
    double x_ = 0;
    double y_ = 0;
    double xy_ = 0;
    regression.x = x;
    regression.y = y;
    std::vector<double> v_y_head;
    for (int i = 0; i < n; ++i) {
        double y_head = regression.m * x[i] + regression.n;
        v_y_head.push_back(y_head); // regression.y_head << y_head;
        regression.sum_err += (y_head - y[i]) * (y_head - y[i]);
        x_ += (x[i] - mean_x) * (x[i] - mean_x);
        y_ += (y[i] - mean_y) * (y[i] - mean_y);
        xy_ += (x[i] - mean_x) * (y[i] - mean_y);
    }
    regression.y_head = Vector::Map(&v_y_head[0], n);
    regression.R = (xy_ / sqrt(x_ * y_)) * (xy_ / sqrt(x_ * y_));

    return regression;
}

inline LinearRegression LeastSquares(const std::vector<double>& x, const std::vector<double>& y)
{
    Vector _x, _y;
    _x = Vector::Map(&x[0], x.size());
    _y = Vector::Map(&y[0], y.size());
    return LeastSquares(_x, _y);
}
}
