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

namespace PeakPick{
    
    struct FitResult
    {
        Vector parameter;
        double sum_error = 0;
        double sum_squared = 0;
        double integral = 0;
    };
    
    class GLFit;

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

                if((m_lock(i) == 2) && (parameter(i) < 0))
                {
                    err += -1*parameter(i)*1e10;
                    parameter(i) = parameter(i)*-1;
                    std::cout << i << " adding error for scaling " << parameter(i) << std::endl;
                }

                if((m_lock(i) == 3) && (parameter(i) < 0))
                {
                    err += 0.1;
                    std::cout << i << " adding error" << std::endl;
                }
            }
            
            err /= double(fvec.size());
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
        
        PeakPick::GLFit *m_glfit;
    };

    struct LiberalGLFitNumericalDiff : Eigen::NumericalDiff<LiberalGLFit> {};
    
    class GLFit 
    {
        
    public:
        inline GLFit(const spectrum *spec, double start, double end, double ratio = 0.9, int fittype = 1)  : m_spec(spec), m_start(start), m_end(end), m_ratio(ratio), m_fittype(fittype) 
        { 
            
        }
        inline ~GLFit() { };
        
        inline void setFitType(int fittype) { m_fittype = fittype; }
        inline void setGLRatio(double ratio) { m_ratio = ratio; }
        
        inline void setParameter(const Vector &parameter) { m_parameter = parameter; }
        
        inline void setGuess(const Vector &guess) 
        { 
            Vector parameter(6*guess.size());
            Vector lock(6*guess.size());
            for(int i = 0; i < guess.size(); ++i)
            {
                parameter(0+i*6) = guess(i);
                if(m_fittype == Conservative)
                    lock(0+i*6) = 1;
                else
                    lock(0+i*6) = 0;

                parameter(1+i*6) = 1;
                lock(1+i*6) = 3;

                parameter(2+i*6) = 10;
                lock(2+i*6) = 3;

                parameter(3+i*6) = 1/double(50);
                lock(3+i*6) = 3;

                parameter(4+i*6) = 1/double(30);
                lock(4+i*6) = 2;

                parameter(5+i*6) = m_ratio;
                if(m_fittype == Innovative)
                    lock(5+i*6) = 0;
                else
                    lock(5+i*6) = 1;
            }
            m_parameter = parameter;
            m_guess = guess;
            m_lock = lock;
        }
        
        inline void UpdateParamater(const Vector &parameter)
        {
            for(int i = 0; i < m_parameter.size(); ++i)
            {
                    if(m_lock(i) == 0)
                        m_parameter(i) = parameter(i);
            }
        }
        
        inline FitResult* Deconvulate()
        {
            LiberalGLFit glfit(6*m_guess.size(), m_end-m_start+1);
            
            glfit.start = m_start;
            glfit.end = m_end;
            glfit.spec = m_spec;
            glfit.m_glfit = this;
            
            glfit.m_original = m_parameter;
            glfit.m_lock = m_lock;
            
            Eigen::NumericalDiff<LiberalGLFit> numDiff(glfit);
            Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LiberalGLFit> > lm(numDiff);
            Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(m_parameter);
            
            lm.minimize(m_parameter);
            //glfit(m_parameter);
            
            FitResult *result = new FitResult;
            result->parameter = m_parameter;
            result->sum_error = glfit.m_sum_error;
            result->sum_squared = glfit.m_sum_squared;
            result->integral = IntegrateGLFunction(m_parameter);
            // std::cout << parameter << std::endl;
            return result;
        }
        
    private:
        inline void releaseLock()
        {
            Vector lock(m_parameter.size());
            for(int i = 0; i < m_parameter.size(); ++i)
                lock(i) = 0;
            m_lock = lock;
        }
        
        inline void createLock()
        {
            Vector lock(m_parameter.size());
            for(int i = 0; i < m_parameter.size(); ++i)
                lock(i) = 1;
            m_lock = lock;
        }
        
        inline void Lock()
        {
            Vector lock(m_parameter.size());
            for(int i = 0; i < m_parameter.size()/6; ++i)
            {
                if(m_fittype == Conservative)
                    lock(0+i*6) = 1;
                else
                    lock(0+i*6) = 0;

                lock(1+i*6) = 3;

                lock(2+i*6) = 3;

                lock(3+i*6) = 3;

                lock(4+i*6) = 2;

                if(m_fittype == Innovative)
                    lock(5+i*6) = 0;
                else
                    lock(5+i*6) = 1;
            }
            m_lock = lock;
        }
        
        Vector m_parameter, m_guess, m_lock;
        int m_fittype;
        const spectrum *m_spec;
        double m_start, m_end, m_ratio;
    };
}

