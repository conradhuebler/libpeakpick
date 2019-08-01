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

#include <fstream>
#include <iostream>

#include "analyse.h"
#include "math.h"
#include "spectrum.h"

inline PeakPick::spectrum loadFromFile(const std::string& filename, double min = 0, double max = 0)
{
    Vector x, y;
    std::ifstream myfile;
    std::string line;
    myfile.open(filename);
    std::string ignore = "#";
    std::string start = "#start = ";
    std::string end = "#end = ";
    if (myfile.is_open()) {
        std::vector<double> entries_x, entries_y;
        int rows(0);
        while (getline(myfile, line)) {

            if (line.find(ignore) != std::string::npos)
                continue;

            entries_y.push_back(std::stod(line));
            entries_x.push_back(min + rows);
            rows++;
        }
        y = Vector::Map(&entries_y[0], rows); //[0], rows, 1);
        x = Vector::Map(&entries_x[0], rows);
        myfile.close();
    } else
        std::cout << "Unable to open file" << std::endl;
    return PeakPick::spectrum(x, y);
}
