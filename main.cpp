/*
 * <PeakPick Example file - How to use it quick 'n' dirty.>
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

#include "libpeakpick/peakpick.h"

#include <iostream>

#include "src/testLorentzian.h"

int FirstStart() 
{
    
    Vector y(10);
    
    y << 1,1,2,3,4,4,5,4,2,1;

    PeakPick::spectrum spec(y, 1, 10);

    std::cout << "Spectrum with " << spec.Mean() << " as mean. The maximal value is ("<< spec.PosMax() + 1 << "," << spec.Max() << ") and the minimal is ("<< spec.PosMin() + 1 << "," << spec.Min() << "). The stddev " << spec.StdDev() << ". Fine" << std::endl;
    std::cout << y << std::endl;
        
    y = Eigen::VectorXd::Random(10);

    spec = PeakPick::spectrum(y,1,10);

    std::cout << "Spectrum with " << spec.Mean() << " as mean. The maximal value is ("<< spec.PosMax() + 1 << "," << spec.Max() << ") and the minimal is ("<< spec.PosMin() + 1 << "," << spec.Min() << "). The stddev " << spec.StdDev() << ". Fine" << std::endl;
    std::cout << y << std::endl;
    
    y = Eigen::VectorXd::Random(10);
    spec = PeakPick::spectrum(y, -5, 5);

    std::cout << "Spectrum with " << spec.Mean() << " as mean. The maximal value is ("<< spec.PosMax() + 1 << "," << spec.Max() << ") and the minimal is ("<< spec.PosMin() + 1 << "," << spec.Min() << "). The stddev " << spec.StdDev() << ". Fine" << std::endl;
    std::cout << y << std::endl;
    
    
    return 0;
}

int main()
{

    if(FirstStart() == 0)
        std::cout << "libpeak basis working " << std::endl;
    
    testLorentzien();
    
    return 0;
}
