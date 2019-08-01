/*
 * <Spectrum Header file.>
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

#include <cmath>
#include <iostream>
#include <vector>

#include "mathhelper.h"

typedef Eigen::VectorXd Vector;

namespace PeakPick {

enum { Innovative = 1,
    Liberal = 2,
    Conservative = 3 };

class spectrum {

public:
    inline spectrum(const Vector& x, const Vector& y)
        : m_x(x)
        , m_y(y)
    {
        if (x.size() != y.size())
            throw 2;

        Analyse();
    }

    inline spectrum() {}

    inline spectrum(const spectrum* other)
    {

        m_x = other->m_x;
        m_y = other->m_y;
        Analyse();
    }
    inline spectrum(const spectrum& other)
    {
        m_x = other.m_x;
        m_y = other.m_y;
        Analyse();
    }

    inline ~spectrum() {}

    inline void setSpectrum(const Vector& x, const Vector& y)
    {
        if (x.size() != y.size())
            throw 2;

        m_x = x;
        m_y = y;
        Analyse();
    }

    inline void setSpectrum(const Vector& y)
    {
        if (y.size() != m_y.size())
            throw 3;

        m_y = y;
        Analyse();
    }

    inline void Analyse()
    {
        m_mean = mean(m_y, &m_pos_min, &m_pos_max);
        m_stddev = stddev(m_y, m_mean);
    }

    inline Vector getRangedSpectrum(unsigned int start, unsigned int end)
    {
        Vector vector(end - start);
        if (start >= m_y.size() || end >= m_y.size())
            return vector;
        for (unsigned int i = start; i < end; ++i)
            vector(i) = Y(i);
        return vector;
    }

    inline Vector getRangedSpectrum(double start, double end)
    {
        std::vector<double> entries;

        if (start >= m_y.size() || end >= m_y.size())
            return Vector(0);
        int number = 0;
        for (unsigned int i = 0; i < m_y.size(); ++i) {
            double val = X(i);
            if (val <= end && val >= start) {
                entries.push_back(Y(i));
                number++;
            }
        }
        Vector vector = Vector::Map(&entries[0], number);
        return vector;
    }

    inline double Mean() const { return m_mean; }
    inline double Max() const { return m_y(m_pos_max); }
    inline int IndexMax() const { return m_pos_max; }
    inline double PosMax() const { return X(m_pos_max); }
    inline unsigned int IndexMin() const { return m_pos_min; }
    inline double Min() const { return m_y(m_pos_min); }
    inline double PosMin() const { return X(m_pos_min); }
    inline double StdDev() const { return m_stddev; }
    inline double Threshold() const { return stddevThreshold(m_y, m_mean, m_stddev); }

    inline double X(unsigned int i) const
    {
        if (i > m_x.size())
            return 0;
        else
            return m_x(i);
    }

    inline double X(int i) const
    {
        if (i > m_x.size() && i < 0)
            return 0;
        else
            return m_x(i);
    }

    inline double Y(unsigned int i) const
    {
        if (i > m_y.size())
            return 0;
        else
            return m_y(i);
    }

    inline double Y(int i) const
    {
        if (i > m_y.size() && i < 0)
            return 0;
        else
            return m_y(i);
    }

    inline double Y(double x) const
    {
        unsigned int i = XtoIndex(x);
        return Y(i);
    }

    inline int XtoIndex(double x) const
    {
        double step = Step();
        double diff = (x - XMin()) / step;
        int val = diff;
        double m_diff = abs(x - m_x[val]);
        for (int i = diff - 4; i < diff + 4 && i < m_x.size(); ++i) {
            // std::cout << i << " " << m_x(i) << std::endl;
            if (abs(x - m_x(i)) < m_diff) {
                val = i;
                m_diff = abs(x - m_x(i));
            }
        }
        return val;
    }

    inline double Step() const
    {
        return m_x[1] - m_x[0];
    }

    inline double XMin() const { return m_x[0]; }
    inline double XMax() const { return m_x[m_x.size() - 1]; }

    inline void setY(unsigned int i, double value) { m_y(i - 1) = value; }

    inline unsigned int size() const { return m_y.size(); }

    inline void setZero(unsigned int start, unsigned int end)
    {
        for (unsigned int i = start; i < end; ++i)
            setY(i, 0);
    }

    inline spectrum* operator=(const spectrum* other)
    {
        m_x = other->m_x;
        m_y = other->m_y;
        Analyse();
        return this;
    }
    inline spectrum& operator=(const spectrum& other)
    {
        m_x = other.m_x;
        m_y = other.m_y;
        Analyse();
        return *this;
    }

    inline void print() const
    {
        //std::cout << m_y << std::endl;
        unsigned int i = 1;
        double step = Step();
        std::cout << "Step size " << step << " starting from " << XMin() << " to " << XMax() << " in " << size() << " steps." << std::endl;
        for (double x = XMin(); x <= XMax(); x += step) {
            std::cout << i << " " << Y(i) << " (" << x << "," << Y(x) << ") " << XtoIndex(x) << std::endl;
            //            std::cout << "l" << x << " " << XtoIndex(x) << std::endl;
            ++i;
        }
    }

    void center()
    {
        for (unsigned int i = 0; i < m_y.size(); ++i)
            m_y[i] -= m_mean;
    }

    void InvertSgn()
    {
        for (unsigned int i = 0; i < m_y.size(); ++i)
            m_y[i] *= -1;
    }

    inline double LastY() const
    {
        return m_y[m_y.size() - 1];
    }

    inline Vector x() const { return m_x; }
    inline Vector y() const { return m_y; }

private:
    Vector m_y;
    Vector m_x;

    double m_mean, m_stddev;
    unsigned int m_pos_min, m_pos_max;
};
}
