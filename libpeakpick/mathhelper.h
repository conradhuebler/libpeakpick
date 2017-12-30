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


#include <cmath>
#include <functional>
#include <iostream>
#include <vector>

typedef Eigen::VectorXd Vector;
static double pi = 3.14159265;

namespace PeakPick{

    inline double SimpsonIntegrate(double lower, double upper, std::function<double(double, const Vector &)> function, const Vector &parameter, double prec = 1e-3)
    {
        double integ = 0;
        double delta = prec;

        int increments = (upper - lower)/prec;
        std::cout << "prepare for numerical integration " << integ << parameter.transpose() << std::endl;
#pragma omp parallel for reduction(+:integ)
         for(int i = 0; i < increments; ++i)
        {
            double x = lower+i/double(increments);
            double b = x + delta;
            integ += (b-x)/6*(function(x, parameter)+4*function((x+b)/2,parameter)+function(b,parameter));
        }
        std::cout << "result from integration over part of function " << integ << parameter.transpose() << std::endl;

        return integ;
    }


    inline double mean(const Vector &vector, int *min = NULL, int *max = NULL)
    {
        if(vector.size() == 0)
        {
            std::cout << "PeakPick::mean() vector is empty" << std::endl;
            return 0;
        }
        double sum = 0;
        double temp_min = vector(0);
        int pos_min = 0;
        double temp_max = vector(0);
        int pos_max = 0;
// #pragma omp parallel for shared(sum, temp_min, pos_min, temp_max, pos_max) reduction(+:sum, &:temp_min, &:pos_min, &:temp_max, &:pos_max)
        for(int i = 0; i < vector.size(); ++i)
        {
            sum += vector(i);
            if(vector(i) > temp_max)
            {
                temp_max = vector(i);
                pos_max = i;
            }
            if(vector(i) < temp_min)
            {
                temp_min = vector(i);
                pos_min = i;
            }
            
            if(max != NULL)
                *max = pos_max;
            
            if(min != NULL)
                *min = pos_min;
        }
        return sum/double(vector.size());
    }
    
    inline double meanThreshold(const Vector &vector, double threshold, int *min = NULL, int *max = NULL)
    {
        if(vector.size() == 0)
            return 0;
        
        double sum = 0;
        double temp_min = vector(0);
        int pos_min = 0;
        double temp_max = vector(0);
        int pos_max = 0;
        for(int i = 0; i < vector.size(); ++i)
        {
            if(std::abs(vector(i)) >= threshold)
                continue;
            sum += vector(i);
            if(vector(i) > temp_max)
            {
                temp_max = vector(i);
                pos_max = i;
            }
            if(vector(i) < temp_min)
            {
                temp_min = vector(i);
                pos_min = i;
            }
            
            if(max != NULL)
                *max = pos_max;
            
            if(min != NULL)
                *min = pos_min;
        }
        
        return sum/double(vector.size());
    }
    
    
    inline double stddev(const Vector &vector, double mean)
    {
        if(vector.size() == 0)
            return 0;
        
        double sum = 0;
#pragma omp parallel for reduction(+:sum)
        for(int i = 0; i < vector.size(); ++i)
        {
            sum += (vector(i) - mean)*(vector(i) - mean);
        }
        return sqrt(sum/double(vector.size()));
    }
    
    inline double stddevThreshold(const Vector &vector, double mean, double threshold)
    {
        if(vector.size() == 0)
            return 0;
        
        double sum = 0;
        
#pragma omp parallel for reduction(+:sum)
        for(int i = 0; i < vector.size(); ++i)
        {
            if(std::abs(vector(i)) >= threshold)
                continue;
            sum += (vector(i) - mean)* (vector(i) - mean);
        }
        return sqrt(sum/double(vector.size()));
    }

        
    inline double Gaussian(double x, double a, double x_0, double c)
    {
        return a*exp(-pow((x-x_0),2)/(2*c*c));
    }
    
    inline double Lorentzian(double x, double x_0, double gamma)
    {
        return 1/pi*(0.5*gamma)/(pow(x-x_0,2)+pow(0.5*gamma,2));
    }
    
    inline double Signal(double x, const Vector &parameter)
    {
        double signal = 0;
        for(int i = 0; i < parameter.size()/6; ++i)
        {
            double gaussian = Gaussian(x, parameter(1+i*6), parameter(0+i*6), parameter(2+i*6));
            double lorentzian = Lorentzian(x, parameter(0+i*6), parameter(3+i*6));
            signal += ((1-parameter(5+i*6))*gaussian + parameter(5+i*6)*lorentzian)*parameter(4+i*6);
        }
        return signal;
    }

    
    inline double SignalSingle(double x, const Vector &parameter, int function)
    {
        double signal = 0;
        if(function >= parameter.size()/6)
            return 0;
        
        double gaussian = Gaussian(x, parameter(1+function*6), parameter(0+function*6), parameter(2+function*6));
        double lorentzian = Lorentzian(x, parameter(0+function*6), parameter(3+function*6));
        signal += ((1-parameter(5+function*6))*gaussian + parameter(5+function*6)*lorentzian)*parameter(4+function*6);
        
        return signal;
    }
    
    
    inline double IntegrateGLFunction(const Vector &parameter)
    {
        double integ = 0;
        for(int i = 0; i < parameter.size()/6; ++i)
        {
            integ += (((1-parameter(5+i*6))*parameter(1+i*6)*parameter(2+i*6))+parameter(5+i*6))*parameter(4+i*6);
        }
        std::cout << "result from integration over whole function " << integ << parameter.transpose() << std::endl;
        return integ;
    }

    inline double IntegrateGLSignal(const Vector &parameter, double start, double end)
    {
        std::function<double(double, const Vector &)> function = Signal;
        return SimpsonIntegrate(start, end, function, parameter);
    }
}

