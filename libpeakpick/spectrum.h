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
    inline spectrum(const Vector& y, double min, double max)
        : m_y(y)
        , m_min(min)
        , m_max(max)
    {
        Analyse();
    }

    inline spectrum() {}

    inline spectrum(const spectrum* other)
    {
        m_y = other->m_y;
        m_min = other->m_min;
        m_max = other->m_max;
        Analyse();
    }
    inline spectrum(const spectrum& other)
    {
        m_y = other.m_y;
        m_min = other.m_min;
        m_max = other.m_max;
        Analyse();
    }

    inline ~spectrum() {}

    inline void setSpectrum(const Vector& y, double min, double max)
    {
        m_y = y;
        m_min = min;
        m_max = max;
        Analyse();
    }

    inline void setSpectrum(const Vector& y)
    {
        if (y.size() == m_y.size())
            m_y = y;
        else {
            std::cout << "not updating data" << std::endl;
            return;
        }
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
#ifdef X
        if (i >= m_x.size() || i < 0)
            return 0;
        else
            return m_x(i);
#else
        return XFromIndex(i);
#endif
    }
    inline double Y(unsigned int i) const
    {
        if (i > m_y.size() || i == 0)
            return 0;
        else
            return m_y(i - 1);
    }

    inline double Y(int i) const
    {
        if (i > m_y.size() || i <= 0)
            return 0;
        else
            return m_y(i - 1);
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
        double val = diff + 1;
        //         std::cout  << x << " "<< step << " " << x - XMin() << " " << diff  << " hier:"  <<val << std::endl;
        //         std::cout << "dort: " <<val << std::endl;
        return val;
    }

    inline double XFromIndex(unsigned int index) const
    {
        //         std::cout << "index" <<index << " XMin" << XMin() << " iStep" << index*Step() << std::endl;
        return XMin() + (index - 1) * Step();
    }

    inline double Step() const
    {
        return (XMax() - XMin()) / double(m_y.size() - 1);
    }

    inline double XMin() const { return m_min; }
    inline double XMax() const { return m_max; }

    inline void setY(unsigned int i, double value) { m_y(i - 1) = value; }

    inline unsigned int size() const { return m_y.size(); }

    inline void setZero(unsigned int start, unsigned int end)
    {
        for (unsigned int i = start; i < end; ++i)
            setY(i, 0);
    }

    inline spectrum* operator=(const spectrum* other)
    {
        m_y = other->m_y;
        m_min = other->m_min;
        m_max = other->m_max;
        Analyse();
        return this;
    }
    inline spectrum& operator=(const spectrum& other)
    {
        m_y = other.m_y;
        m_min = other.m_min;
        m_max = other.m_max;
        Analyse();
        return *this;
    }

    inline void print() const
    {
        std::cout << m_y << std::endl;
        unsigned int i = 0;
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

private:
    Vector m_y;
#ifdef X
    Vector m_x;
#endif
    double m_mean, m_stddev;
    double m_min, m_max;
    unsigned int m_pos_min, m_pos_max;
};
}
