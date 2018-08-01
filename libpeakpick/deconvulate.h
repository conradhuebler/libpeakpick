/*
 * <Math containing Header file.>
 * Copyright (C) 2017  Conrad Hübler <Conrad.Huebler@gmx.net>
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
#include "glfit.h"
#include "mathhelper.h"
#include "spectrum.h"

typedef Eigen::VectorXd Vector;

namespace PeakPick {

struct FitResult {
    Vector parameter;
    double sum_error = 0;
    double sum_squared = 0;
    double integral = 0;
};

template <typename _Scalar, unsigned int NX = Eigen::Dynamic, unsigned int NY = Eigen::Dynamic>

struct GaussianLorentzian {
    typedef _Scalar Scalar;
    enum {
        InputsAtCompileTime = NX,
        ValuesAtCompileTime = NY
    };
    typedef Eigen::Matrix<Scalar, InputsAtCompileTime, 1> InputType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, 1> ValueType;
    typedef Eigen::Matrix<Scalar, ValuesAtCompileTime, InputsAtCompileTime> JacobianType;

    unsigned int m_inputs, m_values;

    inline GaussianLorentzian(unsigned int inputs, unsigned int values)
        : m_inputs(inputs)
        , m_values(values)
    {
    }

    unsigned int inputs() const { return m_inputs; }
    unsigned int values() const { return m_values; }
};

struct LiberalGLFit : GaussianLorentzian<double> {
    inline LiberalGLFit(GLFit* glfit)
        : m_glfit(glfit)
        , GaussianLorentzian(glfit->ParamSize(), glfit->PointsSize())
        , no_parameter(glfit->ParamSize())
        , no_points(glfit->PointsSize())
        , start(glfit->Start())
        , end(glfit->End())
    {
    }

    inline ~LiberalGLFit() {}

    inline int operator()(const Eigen::VectorXd parameter, Eigen::VectorXd& fvec) const
    {
        unsigned int j = 0;
        m_glfit->UpdateParamater(parameter);
        for (unsigned int i = start; i <= end; ++i) {
            fvec(j) = (*m_glfit)(i);
            ++j;
        }
        return 0;
    }

    unsigned int no_parameter, no_points, start, end;
    inline unsigned int inputs() const { return no_parameter; }
    inline unsigned int values() const { return no_points; }

    PeakPick::GLFit* m_glfit;
};

struct LiberalGLFitNumericalDiff : Eigen::NumericalDiff<LiberalGLFit> {
};

inline FitResult* Deconvulate(GLFit* glfit)
{
    LiberalGLFit fit(glfit);
    Vector parameter = glfit->Parameter();

    Eigen::NumericalDiff<LiberalGLFit> numDiff(fit);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LiberalGLFit>> lm(numDiff);
    Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(parameter);

    lm.minimize(parameter);

    FitResult* result = new FitResult;
    result->parameter = parameter;
    result->sum_error = glfit->SumError();
    result->sum_squared = glfit->SumSquared();
    glfit->Print();
    result->integral = IntegrateGLFunction(parameter);

    return result;
}
}
