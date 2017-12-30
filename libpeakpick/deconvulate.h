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
    
    class GLFit
    {
        
    public:
        inline GLFit(const spectrum *spec, double start, double end, double ratio = 0.9, int fittype = 1)  : m_spec(spec), m_start(start), m_end(end), m_ratio(ratio), m_fittype(fittype) 
        { 
            
        }

        inline ~GLFit() { }
        
        inline void setFitType(int fittype) { m_fittype = fittype; }
        inline void setGLRatio(double ratio) { m_ratio = ratio; }
        
        inline void setParameter(const Vector &parameter) 
        { 
            m_parameter = parameter; 
            Lock();
        }
        
        inline Vector Parameter() const { return m_parameter; }
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
            m_sum_error = 0;
            m_sum_squared = 0;
            m_err = 0;
            for(int i = 0; i < parameter.size(); ++i)
            {
                if(m_lock(i) == 0 || m_lock(i) == 2 || m_lock(i) == 3)
                    m_parameter(i) = parameter(i);

                if((m_lock(i) == 2) && (parameter(i) < 0))
                {
                    m_err += -1*parameter(i)*1e10;
                    m_parameter(i) = parameter(i)*-1;
                    std::cout << i << " adding error for scaling " << parameter(i) << std::endl;
                }

                if((m_lock(i) == 3) && (parameter(i) < 0))
                {
                    m_err += 0.1;
                    std::cout << i << " adding error" << std::endl;
                }
            }
            m_err /= double(m_end-m_start+1); 
        }
        
        inline void setLock(const Vector &lock)
        {
            m_lock = lock; 
        }
        
        inline Vector Lock() const { return m_lock; }
    
        inline double operator()(int i)
        {
            double x = m_spec->X(i);
            double Y = m_spec->Y(i);
            double Y_=  Signal(x, m_parameter);
            
            m_sum_error += Y-Y_;
            m_sum_squared += (Y-Y_)*(Y-Y_);
            
            return (Y - Y_  + m_err);
        }
        
        inline int PointsSize() const { return m_end-m_start+1; }
        inline int ParamSize() const { return m_parameter.size(); }
        inline double Start() const { return m_start; }
        inline double End() const { return m_end; }
        inline double SumError() const { return m_sum_error; }
        inline double SumSquared() const { return m_sum_squared; }
        const spectrum  *Spec() const { return m_spec; }
        
        void Print() const
        {
         std::cout << "Gauss/Lorentz Parameter:" << std::endl;
         for(int i = 0; i < m_parameter.size()/6; ++i)
         {
             double a = m_parameter(1+i*6);
             double x_0 = m_parameter(0+i*6);
             double c = m_parameter(2+i*6);
             double gamma = m_parameter(3+i*6);
             std::cout << i + 1 << ": Function as = ( (" << 1-m_parameter(5+i*6) << "*Gaussian) + (" << m_parameter(5+i*6) << "*Lorentzian) )*" << m_parameter(4+i*6) << std::endl;
             std::cout << i + 1 << ": Gaussian type function = 1/(" << a << "*sqrt(2*pi))*exp(-pow((x-" << x_0 << "),2)/(2*pow(" <<c << ",2)))" << std::endl;
             std::cout << i + 1 << ": Lorentzian type function = 1/pi*(0.5*" << gamma << ")/(pow(x-"<< x_0 << ",2)+pow(0.5*1" <<gamma<<",2))" << std::endl;
         }
            
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
        double m_start, m_end, m_ratio, m_err, m_sum_error, m_sum_squared;

    };
    
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
        inline LiberalGLFit(GLFit *glfit) : m_glfit(glfit), GaussianLorentzian(glfit->ParamSize(), glfit->PointsSize()), no_parameter(glfit->ParamSize()),  no_points(glfit->PointsSize()), start(glfit->Start()), end(glfit->End())
        {
            
        }
        
        inline ~LiberalGLFit() { }
        
        inline int operator()(const Eigen::VectorXd parameter, Eigen::VectorXd &fvec) const
        {
            int j = 0;
            m_glfit->UpdateParamater(parameter);
            for(int i = start; i <= end; ++i)
            {
                 fvec(j) =  (*m_glfit)(i); 
                ++j;
            }
            return 0;
        }

        int no_parameter, no_points, start, end;
        inline int inputs() const { return no_parameter; }
        inline int values() const { return no_points; }
        
        PeakPick::GLFit *m_glfit;
    };

    struct LiberalGLFitNumericalDiff : Eigen::NumericalDiff<LiberalGLFit> {};

    
    inline FitResult* Deconvulate(GLFit *glfit)
        {
            LiberalGLFit fit(glfit);
            Vector parameter = glfit->Parameter();
            
            Eigen::NumericalDiff<LiberalGLFit> numDiff(fit);
            Eigen::LevenbergMarquardt<Eigen::NumericalDiff<LiberalGLFit> > lm(numDiff);
            Eigen::LevenbergMarquardtSpace::Status status = lm.minimizeInit(parameter);
            
            lm.minimize(parameter);
            
            FitResult *result = new FitResult;
            result->parameter = parameter;
            result->sum_error = glfit->SumError();
            result->sum_squared = glfit->SumSquared();
            glfit->Print();
            result->integral = IntegrateGLFunction(parameter);

            return result;
        }
        
}

