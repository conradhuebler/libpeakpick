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
#include "spectrum.h"
#include "mathhelper.h"

typedef Eigen::VectorXd Vector;


struct FitResult
{
    Vector parameter;
    double sum_error = 0;
    double sum_squared = 0;
    double integral = 0;
    
};

namespace PeakPick{

    template<typename _Scalar, int NX = Eigen::Dynamic, int NY = Eigen::Dynamic>
    
    struct GaussianLorentzian
    {
        typedef _Scalar Scalar;
        enum {
            InputsAtCompileTime = NX,
            ValuesAtCompileTime = NY
        };
        typedef Eigen::Matrix<Scalar,InputsAtCompileTime,1> InputType;
        typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,1> ValueType;
        typedef Eigen::Matrix<Scalar,ValuesAtCompileTime,InputsAtCompileTime> JacobianType;
        
        int m_inputs, m_values;
        
        inline GaussianLorentzian(int inputs, int values) : m_inputs(inputs), m_values(values) {}
        
        int inputs() const { return m_inputs; }
        int values() const { return m_values; }
        
    };
    
    struct LiberalGLFit : GaussianLorentzian<double>
    {
        inline LiberalGLFit(int inputs, int values) : GaussianLorentzian(inputs, values), no_parameter(inputs),  no_points(values)
        {

        }
        inline ~LiberalGLFit() { }
        inline int operator()(Eigen::VectorXd parameter, Eigen::VectorXd &fvec) const
        {
            int j = 0;
            double err = 0;
            for(int i = 0; i < parameter.size(); ++i)
            {
                if(m_lock(i) == 1)
                    parameter(i) = m_original(i);

            }
            for(int i = start; i <= end; ++i)
            {
                double x = spec->X(i);
                double Y = spec->Y(i);
                double Y_=  Signal(x, parameter);
                fvec(j) =  Y - Y_  + err;
                ++j;
            }
            return 0;
        }

        inline int operator()(const Eigen::VectorXd &parameter)
        {
            m_sum_error = 0;
            m_sum_squared = 0;
            for(int i = start; i <= end; ++i)
            {
                double x = spec->X(i);
                double Y = spec->Y(i);
                double Y_=  Signal(x, parameter);
                m_sum_error += Y-Y_;
                m_sum_squared += (Y-Y_)*(Y-Y_);
            }
            return 0;
        }

        int no_parameter, no_points, start, end;
        double m_sum_error, m_sum_squared;
        const spectrum *spec;
        Vector m_original, m_lock;
        inline int inputs() const { return no_parameter; }
        inline int values() const { return no_points; }
    };

    struct LiberalGLFitNumericalDiff : Eigen::NumericalDiff<LiberalGLFit> {};
    
    FitResult* LiberalDeconvulate(const spectrum *spec, double start, double end, double ratio, const Vector &guess, int fittype)
    {
        std::cout << "Having guess size: " << guess.size() << std::endl;
        LiberalGLFit functor(6*guess.size(), end-start+1);
        functor.start = start;
        functor.end = end;
        functor.spec = spec;
        Vector parameter(6*guess.size());
        Vector lock(6*guess.size());
        for(int i = 0; i < guess.size(); ++i)
        {
            parameter(0+i*6) = guess(i);
            if(fittype == Conservative)
                lock(0+i*6) = 1;
            else
                lock(0+i*6) = 0;

            parameter(1+i*6) = 1;
            lock(1+i*6) = 0;

            parameter(2+i*6) = 10;
            lock(2+i*6) = 0;

            parameter(3+i*6) = 1/double(50);
            lock(3+i*6) = 0;

            parameter(4+i*6) = 1/double(30);
            lock(4+i*6) = 0;

            parameter(5+i*6) = ratio;
            if(fittype == Innovative)
                lock(5+i*6) = 0;
            else
                lock(5+i*6) = 1;
        }

        functor.m_original = parameter;
        functor.m_lock = lock;
        Eigen::NumericalDiff<LiberalGLFit> numDiff(functor);
        Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LiberalGLFit> > lm(numDiff);
        Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(parameter);
        lm.minimize(parameter);
        functor(parameter);
        FitResult *result = new FitResult;
        result->parameter = parameter;
        result->sum_error = functor.m_sum_error;
        result->sum_squared = functor.m_sum_squared;
        result->integral = IntegrateGLFunction(parameter);
        std::cout << parameter << std::endl;
        return result;
    }
}

