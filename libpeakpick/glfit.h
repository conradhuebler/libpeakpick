/*
 * <Math containing Header file.>
 * Copyright (C) 2017 - 2020 Conrad Hübler <Conrad.Huebler@gmx.net>
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

typedef Eigen::VectorXd Vector;

namespace PeakPick {

class GLFit {

public:
    inline GLFit(const spectrum* spec, double start, double end, double ratio = 0.9, unsigned int fittype = 1)
        : m_spec(spec)
        , m_fittype(fittype)
        , m_start(start)
        , m_end(end)
        , m_ratio(ratio)

    {
    }

    inline ~GLFit() {}

    inline void setFitType(int fittype) { m_fittype = fittype; }
    inline void setGLRatio(double ratio) { m_ratio = ratio; }

    inline void setParameter(const Vector& parameter)
    {
        m_parameter = parameter;
        Lock();
    }

    inline Vector Parameter() const { return m_parameter; }

    inline void UpdateParamater(const Vector& parameter)
    {
        m_sum_error = 0;
        m_sum_squared = 0;
        m_err = 0;
        for (int i = 0; i < parameter.size(); ++i) {
            if (m_lock(i) == 0 || m_lock(i) == 2 || m_lock(i) == 3)
                m_parameter(i) = parameter(i);

            if ((m_lock(i) == 2) && (parameter(i) < 0)) {
                m_err += -1 * parameter(i) * 1e10;
                m_parameter(i) = parameter(i) * -1;
                std::cout << i << " adding error for scaling " << parameter(i) << std::endl;
            }

            if ((m_lock(i) == 3) && (parameter(i) < 0)) {
                m_err += 0.1;
                std::cout << i << " adding error" << std::endl;
            }
        }
        m_err /= double(m_end - m_start + 1);
    }

    inline void setLock(const Vector& lock)
    {
        m_lock = lock;
    }

    inline Vector Lock() const { return m_lock; }

    inline Vector StdLock(int fittype) const
    {
        Vector lock(6);
        if (fittype == Conservative)
            lock(0) = 1;
        else
            lock(0) = 0;
        lock(1) = 3;
        lock(2) = 3;
        lock(3) = 3;
        lock(4) = 2;
        if (fittype == Innovative)
            lock(5) = 0;
        else
            lock(5) = 1;
        return lock;
    }

    inline Vector StdGuess(double x) const
    {
        Vector parameter(6);

        parameter(0) = x;
        parameter(1) = 1;
        parameter(2) = 10;
        parameter(3) = 1 / double(50);
        parameter(4) = 1 / double(30);
        parameter(5) = m_ratio;

        return parameter;
    }

    inline void addFunction(const Vector& guess, const Vector& lock)
    {
        int size = guess.size() + m_parameter.size();
        Vector parameter(size);
        Vector _lock(size);

        for (int i = 0; i < m_parameter.size(); ++i) {
            parameter(i) = m_parameter(i);
            _lock(i) = m_lock(i);
        }
        for (int i = 0; i < lock.size(); ++i) {
            parameter(i + m_parameter.size()) = guess(i);
            _lock(i + m_parameter.size()) = lock(i);
        }
        m_parameter = parameter;
        m_lock = _lock;
    }

    inline bool removeFunction(int function)
    {
        if (function >= m_parameter.size() / 6)
            return false;
// #warning not yet implemented
    }

    inline Vector Function(int function) const
    {
        Vector parameter(6);

        if (function >= m_parameter.size() / 6)
            return parameter;
        for (int i = 0; i < 6; ++i)
            parameter(i) = m_parameter(i + function * 6);
        return parameter;
    }

    inline double X_0(int function) const
    {
        if (function >= m_parameter.size() / 6)
            return 0;
        else
            return m_parameter(0 + function * 6);
    }

    inline int Functions() const { return m_parameter.size() / 6; }

    inline double operator()(int i)
    {
// #warning maybe 'i - 1' now
        double x = m_spec->X(i);
        double Y = m_spec->Y(i);
        double Y_ = Signal(x, m_parameter);

        m_sum_error += Y - Y_;
        m_sum_squared += (Y - Y_) * (Y - Y_);

        return (Y - Y_ + m_err);
    }

    inline void setGuess(const Vector& guess)
    {
        for (int i = 0; i < guess.size(); ++i) {
            addFunction(StdGuess(guess(i)), StdLock(m_fittype));
        }
    }

    inline int PointsSize() const { return m_end - m_start + 1; }
    inline unsigned int ParamSize() const { return m_parameter.size(); }
    inline double Start() const { return m_start; }
    inline double End() const { return m_end; }
    inline double SumError() const { return m_sum_error; }
    inline double SumSquared() const { return m_sum_squared; }
    const spectrum* Spec() const { return m_spec; }

    inline void releaseLock()
    {
        Vector lock(m_parameter.size());
        for (int i = 0; i < m_parameter.size(); ++i)
            lock(i) = 0;
        m_lock = lock;
    }

    inline void createLock()
    {
        Vector lock(m_parameter.size());
        for (int i = 0; i < m_parameter.size(); ++i)
            lock(i) = 1;
        m_lock = lock;
    }

    inline void Lock()
    {
        Vector lock(m_parameter.size());
        for (int i = 0; i < m_parameter.size() / 6; ++i) {
            if (m_fittype == Conservative)
                lock(0 + i * 6) = 1;
            else
                lock(0 + i * 6) = 0;

            lock(1 + i * 6) = 3;
            lock(2 + i * 6) = 3;
            lock(3 + i * 6) = 3;
            lock(4 + i * 6) = 2;

            if (m_fittype == Innovative)
                lock(5 + i * 6) = 0;
            else
                lock(5 + i * 6) = 1;
        }
        m_lock = lock;
    }

    void Print() const
    {
        std::cout << "Gauss/Lorentz Parameter:" << std::endl;
        for (int i = 0; i < m_parameter.size() / 6; ++i) {
            double a = m_parameter(1 + i * 6);
            double x_0 = m_parameter(0 + i * 6);
            double c = m_parameter(2 + i * 6);
            double gamma = m_parameter(3 + i * 6);
            std::cout << i + 1 << ": Function as = ( (" << 1 - m_parameter(5 + i * 6) << "*Gaussian) + (" << m_parameter(5 + i * 6) << "*Lorentzian) )*" << m_parameter(4 + i * 6) << std::endl;
            std::cout << i + 1 << ": Gaussian type function = 1/(" << a << "*sqrt(2*pi))*exp(-pow((x-" << x_0 << "),2)/(2*pow(" << c << ",2)))" << std::endl;
            std::cout << i + 1 << ": Lorentzian type function = 1/pi*(0.5*" << gamma << ")/(pow(x-" << x_0 << ",2)+pow(0.5*" << gamma << ",2))" << std::endl;
        }
    }

private:
    const spectrum* m_spec;

    Vector m_parameter, m_guess, m_lock;
    unsigned int m_fittype;
    double m_start, m_end, m_ratio, m_err, m_sum_error, m_sum_squared;
};
}
