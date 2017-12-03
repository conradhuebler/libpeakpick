/*
 * <Main header file for PeakPick Library.>
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

#include <iostream>
#include <fstream>

#include "math.h"
#include "analyse.h"
#include "spectrum.h"


inline PeakPick::spectrum loadFormFile(const std::string &filename, double min = 0, double max = 0)
{
    Vector y;
    std::ifstream myfile;
    std::string line;
    myfile.open (filename); 
    std::string ignore = "#";
    if (myfile.is_open())
    {
        std::vector<double> entries;
        int rows(0);
        while ( getline (myfile,line) )
        {
            rows++;
            if(line.find(ignore) != std::string::npos)
                continue;
            std::cout << line << " " << std::stod(line) << std::endl;
            entries.push_back(std::stod(line));
        }
        y = Vector::Map(&entries[0], rows); //[0], rows, 1);
        myfile.close();
    }
    else 
        std::cout << "Unable to open file" << std::endl;; 
    return PeakPick::spectrum(y,min,max); 
}

