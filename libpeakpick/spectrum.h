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

namespace PeakPick{

    enum { Innovative = 1, Liberal = 2, Conservative = 3 };


class spectrum
{
    
public:
    inline spectrum(const Vector &y, double min, double max) : m_y(y), m_min(min), m_max(max)
    { 
        Analyse();
    }

    inline spectrum() { }
    
    inline spectrum(const spectrum *other) { m_y = other->m_y; m_min = other->m_min;  m_max = other->m_max; Analyse();}
    inline spectrum(const spectrum &other) { m_y = other.m_y; m_min = other.m_min;  m_max = other.m_max; Analyse();  }

    inline ~spectrum() { }
    
    inline void setSpectrum(const Vector &y, double min, double max)
    { 
        m_y = y;
        m_min = min;
        m_max = max;
        Analyse();
    }
    
    inline void setSpectrum(const Vector &y)
    {
        if(y.size() == m_y.size())
            m_y = y;
        else
        {
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
    
    inline Vector getRangedSpectrum(int start, int end)
    {
        Vector vector(end-start);
        if(start >= m_y.size() || end >= m_y.size())
            return vector;
        for(int i = start; i < end; ++i)
            vector(i) = Y(i);
        return vector;
    }
    
    inline Vector getRangedSpectrum(double start, double end)
    {
        std::vector<double> entries;
        
        if(start >= m_y.size() || end >= m_y.size())
            return Vector(0);
        int number = 0;
        for(int i = 0; i <  m_y.size(); ++i)
        {
            double val = X(i);
            if(val <= end && val >=start)
            {
                entries.push_back(Y(i));
                number++;
            }
        }
        Vector vector = Vector::Map(&entries[0], number); 
        return vector;
    }
    
    inline double Mean() const { return m_mean; }
    inline double Max() const { return m_y(m_pos_max); }
    inline double PosMax() const { return X(m_pos_max);       } // m_pos_max*double(m_y.size())/m_max+m_min; }
    inline double Min() const { return m_y(m_pos_min); }
    inline double PosMin() const { return m_pos_min*double(m_y.size())/m_max+m_min; }
    inline double StdDev() const { return m_stddev; }
    inline double Threshold() const { return stddevThreshold(m_y, m_mean, m_stddev); }
    
    inline double X(int i) const { return m_min+(i)*(m_max-m_min)/double(m_y.size()); }
    inline double Y(int i) const 
        { 
            if(i >= m_y.size() || i < 0)
                return 0;
            else 
                return m_y(i); 
        }

    inline double Y(double x) const
    {
        int i = (XMin()+x)*double(m_y.size())/(XMax()-XMin());
        return Y(i);
    }
        
    inline int XtoIndex(double x) const
    {
        return (XMin()+x)*double(m_y.size())/(XMax()-XMin());
    }

    inline double XMin() const { return m_min; }
    inline double XMax() const { return m_max; }
    
    inline double PosOfPoint(double X) const { return (X-m_min)*m_y.size()/(m_max-m_min)/2; }
        
    inline void setY(int i, double value) { m_y(i) = value; }
    
    inline int size() const { return m_y.size(); }
    
    inline void setZero(int start, int end) 
    { 
        for(int i = start; i < end; ++i)
            setY(i, 0);
    }
    
    inline spectrum *operator=(const spectrum *other) { m_y = other->m_y; m_min = other->m_min;  m_max = other->m_max; Analyse(); return this; }
    inline spectrum &operator=(const spectrum &other) { m_y = other.m_y; m_min = other.m_min;  m_max = other.m_max; Analyse();  return *this; }
    inline void print() const
    {
        int i = 0;
        for(double min = XMin(); min < XMax(); min += (XMax()-XMin())/double(m_y.size()-1))
            std::cout << i++ << " " << Y(i) << " " << min << " " << Y(min) << " " << XtoIndex(min) << std::endl;
    }
private:
    Vector m_y;
    double m_mean, m_stddev;
    double m_min, m_max;
    int m_pos_min, m_pos_max;
};


}
