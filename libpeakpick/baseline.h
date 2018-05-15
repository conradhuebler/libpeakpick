/*
 * <BaseLine Header file.>
 * Copyright (C) 2018  Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include <Eigen/Sparse>
#include <unsupported/Eigen/NonLinearOptimization>

#include <cmath>
#include <iostream>
#include <vector>

#include "analyse.h"
#include "mathhelper.h"
#include "spectrum.h"

namespace PeakPick {

template <typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>

struct BaseLineFit {
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    int m_inputs, m_values;

    inline BaseLineFit(int inputs, int values)
        : m_inputs(inputs)
        , m_values(values)
    {
    }

    int inputs() const { return m_inputs; }
    int values() const { return m_values; }
};

struct BaseLineFitFunction : BaseLineFit<double> {
    inline BaseLineFitFunction(const Vector& x, const Vector& y, int size)
        : m_x(x)
        , m_y(y)
        , BaseLineFit(x.size(), size)
        , no_parameter(size)
        , no_points(x.size())
    {
    }

    inline ~BaseLineFitFunction() {}

    inline int operator()(const Eigen::VectorXd parameter, Eigen::VectorXd& fvec) const
    {
        for (int i = 0; i < m_x.size(); ++i)
            fvec(i) = m_y(i) - Polynomial(m_x(i), parameter);
        return 0;
    }

    int no_parameter, no_points, start, end;
    inline int inputs() const { return no_parameter; }
    inline int values() const { return no_points; }
    const Vector &m_x, &m_y;
};

struct BaseLineFitFunctionDiff : Eigen::NumericalDiff<BaseLineFitFunction> {
};

inline Vector FitBaseLine(const Vector& x, const Vector& y, int size, double mean)
{
    BaseLineFitFunction fit(x, y, size);
    Vector parameter(size);
    parameter(0) = mean;
    for (int i = 1; i < size; ++i)
        parameter(i) = pow(10, -2 * (size));

    Eigen::NumericalDiff<BaseLineFitFunction> numDiff(fit);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<BaseLineFitFunction>> lm(numDiff);
    Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(parameter);

    qreal diff = 1;

    for (int iter = 0; iter < 100 && diff > 1e-5; ++iter) {

        Vector param = parameter;
        status = lm.minimizeOneStep(parameter);
        for (int i = 0; i < size; ++i)
            diff += (parameter(i) - param(i)) * (parameter(i) - param(i));
        diff = sqrt(diff);
    }

    return parameter;
}

class BaseLine {
public:
    inline BaseLine(const spectrum* spec)
        : m_spec(spec)
        , m_no_coeff(2)
        , m_lower(0)
        , m_upper(0)
    {
    }

    inline void setNoCoeffs(int no_coeff) { m_no_coeff = no_coeff; }
    inline void setLower(double lower) { m_lower = lower; }
    inline void setUpper(double upper) { m_upper = upper; }

    inline Vector Fit()
    {
        Vector vector(m_no_coeff);
        Vector x, y;
        if (!m_peak_list)
            ApplyFilter(x, y);
        else
            FromPeaks(x, y);

        if (m_no_coeff >= 3) {
            vector = FitBaseLine(x, y, m_no_coeff, m_spec->Mean());
        } else if (m_no_coeff == 2) {

            LinearRegression regression = LeastSquares(x, y);
            vector(0) = regression.n;
            vector(1) = regression.m;
        }

        return vector;
    }
    inline void setPeaks(std::vector<Peak>* peaks)
    {
        m_peaks = peaks;
        m_peak_list = true;
    }
    inline void clearPeaks()
    {
        m_peaks->clear();
        m_peak_list = false;
    }

private:
    void ApplyFilter(Vector& x, Vector& y)
    {
        std::vector<double> v_x, v_y;
        for (int i = 0; i < m_spec->size(); ++i) {
            if (ok(m_spec->Y(i))) {
                v_x.push_back(m_spec->X(i));
                v_y.push_back(m_spec->Y(i));
            }
        }
        x = Vector::Map(&v_x[0], v_x.size());
        y = Vector::Map(&v_y[0], v_y.size());
    }

    void FromPeaks(Vector& x, Vector& y)
    {
        std::vector<double> v_x, v_y;
        for (int i = 1; i < m_peaks->size() - 1; ++i) {
            {
                v_x.push_back(m_spec->X((*m_peaks)[i].start));
                v_x.push_back(m_spec->X((*m_peaks)[i].end));
                v_y.push_back(m_spec->Y((*m_peaks)[i].start));
                v_y.push_back(m_spec->Y((*m_peaks)[i].end));
            }
        }
        x = Vector::Map(&v_x[0], v_x.size());
        y = Vector::Map(&v_y[0], v_y.size());
    }

    inline bool ok(double y)
    {
        return (m_lower <= y && y <= m_upper) || (m_lower == m_upper);
    }

    const spectrum* m_spec;
    int m_no_coeff;
    double m_lower, m_upper;
    bool m_peak_list = false;
    std::vector<Peak>* m_peaks;
};
}
