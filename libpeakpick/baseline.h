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

struct BaseLineResult {
    std::vector<Vector> baselines;
    std::vector<Vector> x_grid_points;
    std::vector<Vector> y_grid_points;
};

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
    inline BaseLineFitFunction(const Vector& x, const Vector& y, unsigned int size)
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
        for (unsigned int i = 0; i < m_x.size(); ++i)
            fvec(i) = m_y(i) - Polynomial(m_x(i), parameter);
        return 0;
    }

    int no_parameter, no_points, start, end;
    inline unsigned int inputs() const { return no_parameter; }
    inline unsigned int values() const { return no_points; }
    const Vector &m_x, &m_y;
};

struct BaseLineFitFunctionDiff : Eigen::NumericalDiff<BaseLineFitFunction> {
};

inline Vector FitBaseLine(const Vector& x, const Vector& y, unsigned int size, double mean, Vector initial = Vector(0))
{
    BaseLineFitFunction fit(x, y, size);
    Vector parameter(size);
    if (initial.size() != parameter.size()) {
        parameter(0) = mean;
        for (unsigned int i = 1; i < size; ++i)
            parameter(i) = pow(10, -2 * int(size));
    } else
        parameter = initial;

    Eigen::NumericalDiff<BaseLineFitFunction> numDiff(fit);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<BaseLineFitFunction>> lm(numDiff);
    Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(parameter);

    qreal diff = 1;

    for (unsigned int iter = 0; iter < 100 && diff > 1e-5; ++iter) {

        Vector param = parameter;
        status = lm.minimizeOneStep(parameter);
        for (unsigned int i = 0; i < size; ++i)
            diff += (parameter(i) - param(i)) * (parameter(i) - param(i));
        diff = sqrt(diff);
    }

    return parameter;
}

inline Vector FitBaseLineMLR(const Vector& x, const Vector& y, unsigned int size)
{
    Vector parameter = Vector::Zero(size);
    if (size > x.size())
        return parameter;

    double max = x[x.size() - 1];

    //std::cout << x << std::endl;
    Eigen::MatrixXd X = Eigen::MatrixXd::Ones(x.size(), size);

    for (int j = 0; j < x.size(); ++j) {
        const double quotient = x(j); ///max;
        for (int i = size - 1; i >= 0; --i)
            X(j, i) *= pow(quotient, i);
    }
    //std::cout << X << std::endl;

    parameter = (X.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(y));
    std::cout << parameter.transpose() << std::endl;
    //for(int i = size - 1; i >= 1; --i)
    //    parameter(i) /= max;
    std::cout << parameter.transpose() << std::endl;
    return parameter;
}

inline Vector FitBaseLineIterative(const Vector& x, const Vector& y, unsigned int size, double mean, Vector initial = Vector(0))
{
    //    std::cout << "Intitial guess: " << initial.transpose() << std::endl;

    Vector increment = initial;
    if (initial.size() > size) {
        increment = Vector(size - 1);
        for (unsigned int i = 0; i < size - 1; ++i)
            increment(i) = initial(i);
        initial = increment;
    }

    if (initial.size() < 1) {
        LinearRegression regression = LeastSquares(x, y);

        initial = Vector(2);
        initial(0) = regression.n;
        initial(1) = regression.m;
    }
    increment = initial;
    if (x.size() == 0 || y.size() == 0)
        return initial;

    Vector old = initial;
    //    std::cout << "Intitial guess created: " << initial.transpose() << std::endl;
    for (unsigned int i = initial.size(); i < size; ++i) {
        if (i < size) {
            initial = Vector(i + 1);
            for (unsigned int j = 0; j < i; ++j)
                initial(j) = increment(j);
            initial(i) = 1e-3;
            //std::cout << "Step: "<< i << ": Guess " << initial.transpose() << std::endl << std::endl;
        }

        increment = FitBaseLine(x, y, i + 1, mean, initial);

        //std::cout << "Step: "<< i << ": Result " << increment.transpose() << std::endl;

        if (increment == initial) {
            //  std::cout << "something converged although should not";
            increment = old;
            break;
        }
        old = increment;
    }

    return increment;
}

class BaseLine {

public:
    enum Type {
        Automatic = 0,
        StartEnd = 1
    };

    enum BLR {
        FullSpectrum = 0,
        PeakWise = 1
    };

    enum Polynom {
        Fast = 0,
        Slow = 1,
        MLR = 2
    };

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

    inline void setBaseLineRange(BLR blr) { m_blr = blr; }

    inline void setPolynomFit(Polynom polynom) { m_polynom = polynom; }

    inline void setInitialBaseLine(const Vector& vector) { m_baseline = vector; }

    inline BaseLineResult Fit(Type type = Automatic)
    {
        m_baselineresult.x_grid_points.clear();
        m_baselineresult.y_grid_points.clear();

        if (m_blr == BLR::FullSpectrum) {
            m_baselineresult.baselines.clear();
            m_baselineresult.baselines.push_back(FitFullSpectrum(type));

        } else if (m_blr == BLR::PeakWise) {

            m_baselineresult.baselines = FitPeakWise();
        }

        return m_baselineresult;
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

    inline spectrum Corrected() const
    {
        std::vector<double> y;

        if (m_baselineresult.baselines.size() == 1) {
            for (unsigned int i = 1; i <= m_spec->size(); ++i) {
                y.push_back(m_spec->Y(i) - Polynomial(m_spec->X(i), m_baseline));
            }
        } else if (m_baselineresult.baselines.size() == m_peaks->size()) {
            for (unsigned int i = 0; i < m_peaks->size(); ++i) {
                for (unsigned int k = m_peaks->at(i).start; k <= m_peaks->at(i).end; ++k)
                    y.push_back(m_spec->Y(k) - Polynomial(m_spec->X(k), m_baselineresult.baselines[i]));
            }
        }

        Vector vec = Vector::Map(&y[0], y.size());
        spectrum spec(vec, m_spec->X(1), m_spec->X(m_spec->size()));
        return spec;
    }

private:
    void ApplyFilter(Vector& x, Vector& y)
    {
        std::vector<double> v_x, v_y;
        for (unsigned int i = 1; i <= m_spec->size(); ++i) {
            if (ok(m_spec->Y(i))) {
                v_x.push_back(m_spec->X(i));
                v_y.push_back(m_spec->Y(i));
            }
        }
        x = Vector::Map(&v_x[0], v_x.size());
        y = Vector::Map(&v_y[0], v_y.size());

        m_baselineresult.x_grid_points.push_back(x);
        m_baselineresult.y_grid_points.push_back(y);
    }

    void FromPeaks(Vector& x, Vector& y)
    {
        std::vector<double> v_x, v_y;
        for (unsigned int i = 0; i < m_peaks->size(); ++i) {
            {
                v_x.push_back(m_spec->X((*m_peaks)[i].start));
                v_x.push_back(m_spec->X((*m_peaks)[i].end));
                v_y.push_back(m_spec->Y((*m_peaks)[i].start));
                v_y.push_back(m_spec->Y((*m_peaks)[i].end));
            }
        }
        x = Vector::Map(&v_x[0], v_x.size());
        y = Vector::Map(&v_y[0], v_y.size());
        m_baselineresult.x_grid_points.push_back(x);
        m_baselineresult.y_grid_points.push_back(y);
    }

    inline bool ok(double y)
    {
        return (m_lower <= y && y <= m_upper) || (m_lower == m_upper);
    }

    Vector FitFullSpectrum(Type type)
    {
        Vector vector(m_no_coeff);
        Vector x, y;
        if (type == Automatic) {
            if (!m_peak_list)
                ApplyFilter(x, y);
            else
                FromPeaks(x, y);
        } else if (type == StartEnd) {
            std::vector<double> t_x, t_y;
            t_x.push_back(0);
            t_y.push_back(m_spec->Y(1));
            t_x.push_back(m_spec->size());
            t_y.push_back(m_spec->LastY());
            x = Vector::Map(&t_x[0], 2);
            y = Vector::Map(&t_y[0], 2);

            m_baselineresult.x_grid_points.push_back(x);
            m_baselineresult.y_grid_points.push_back(y);

            m_no_coeff = 2;
        } else
            throw - 1;

        if (m_no_coeff >= 3) {
            if (m_polynom == Polynom::MLR) {
                vector = FitBaseLineMLR(x, y, m_no_coeff);
                m_baseline = vector;
            } else if (m_polynom == Polynom::Fast) {
                vector = FitBaseLine(x, y, m_no_coeff, m_spec->Mean());
                m_baseline = vector;
            } else {
                vector = FitBaseLineIterative(x, y, m_no_coeff, m_spec->Mean(), m_baseline);
                m_baseline = vector;
            }

        } else if (m_no_coeff == 2) {
            LinearRegression regression = LeastSquares(x, y);

            vector(0) = regression.n;
            vector(1) = regression.m;
        } else if (m_no_coeff == 1)
            vector(0) = mean(y);

        return vector;
    }
    std::vector<Vector> FitPeakWise()
    {
        std::vector<Vector> baseline;
        if (m_no_coeff > 2)
            m_no_coeff = 2;
        for (unsigned int i = 0; i < m_peaks->size(); ++i) {
            {
                Vector vector(m_no_coeff);
                std::vector<double> v_x, v_y;
                v_x.push_back(m_spec->X((*m_peaks)[i].start));
                v_x.push_back(m_spec->X((*m_peaks)[i].end));
                v_y.push_back(m_spec->Y((*m_peaks)[i].start));
                v_y.push_back(m_spec->Y((*m_peaks)[i].end));

                Vector x = Vector::Map(&v_x[0], v_x.size());
                Vector y = Vector::Map(&v_y[0], v_y.size());

                m_baselineresult.x_grid_points.push_back(x);
                m_baselineresult.y_grid_points.push_back(y);
                if (m_no_coeff == 2) {
                    LinearRegression regression = LeastSquares(x, y);
                    vector(0) = regression.n;
                    vector(1) = regression.m;
                } else
                    vector(0) = (y[1] + y[0]) / 2.0;
                baseline.push_back(vector);
            }
        }

        return baseline;
    }

    const spectrum* m_spec;
    unsigned int m_no_coeff;
    double m_lower, m_upper;
    bool m_peak_list = false;
    std::vector<Peak>* m_peaks;
    Vector m_baseline = Vector(0);
    BLR m_blr = BLR::FullSpectrum;
    Polynom m_polynom = Polynom::Slow;
    BaseLineResult m_baselineresult;
};
}
