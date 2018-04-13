/*
 * <Savitzky-Golay-Filter Coeff containing Header file.>
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

namespace PeakPick {

inline double SavitzkyGolayCoefficient(int points, int position)
{
    double coeff = 0;
    if (position > points)
        return 0;

    switch (points) {
    case 2:
        if (position == 0)
            coeff = 2;
        else if (position == 1)
            coeff = 12;
        else
            coeff = -3;
        break;

    case 3:
        if (position == 0)
            coeff = 7;
        else if (position == 1)
            coeff = 6;
        else if (position == 2)
            coeff = 3;
        else
            coeff = -2;
        break;

    case 4:
        if (position == 0)
            coeff = 59;
        else if (position == 1)
            coeff = 54;
        else if (position == 2)
            coeff = 39;
        else if (position == 3)
            coeff = 14;
        else
            coeff = -21;
        break;

    case 5:
        if (position == 0)
            coeff = 89;
        else if (position == 1)
            coeff = 84;
        else if (position == 2)
            coeff = 69;
        else if (position == 3)
            coeff = 44;
        else if (position == 4)
            coeff = 9;
        else
            coeff = -36;
        break;

    case 6:
        if (position == 0)
            coeff = 25;
        else if (position == 1)
            coeff = 24;
        else if (position == 2)
            coeff = 21;
        else if (position == 3)
            coeff = 16;
        else if (position == 4)
            coeff = 9;
        else if (position == 5)
            coeff = 0;
        else
            coeff = -11;
        break;

    case 7:
        if (position == 0)
            coeff = 167;
        else if (position == 1)
            coeff = 162;
        else if (position == 2)
            coeff = 147;
        else if (position == 3)
            coeff = 122;
        else if (position == 4)
            coeff = 87;
        else if (position == 5)
            coeff = 42;
        else if (position == 6)
            coeff = -13;
        else
            coeff = -78;
        break;

    case 8:
        if (position == 0)
            coeff = 43;
        else if (position == 1)
            coeff = 42;
        else if (position == 2)
            coeff = 39;
        else if (position == 3)
            coeff = 34;
        else if (position == 4)
            coeff = 27;
        else if (position == 5)
            coeff = 18;
        else if (position == 6)
            coeff = 7;
        else if (position == 7)
            coeff = -6;
        else
            coeff = -21;
        break;

    case 9:
        if (position == 0)
            coeff = 269;
        else if (position == 1)
            coeff = 264;
        else if (position == 2)
            coeff = 249;
        else if (position == 3)
            coeff = 224;
        else if (position == 4)
            coeff = 189;
        else if (position == 5)
            coeff = 144;
        else if (position == 6)
            coeff = 89;
        else if (position == 7)
            coeff = 24;
        else if (position == 8)
            coeff = -51;
        else
            coeff = -136;

    case 10:
        if (position == 0)
            coeff = 329;
        else if (position == 1)
            coeff = 324;
        else if (position == 2)
            coeff = 309;
        else if (position == 3)
            coeff = 284;
        else if (position == 4)
            coeff = 249;
        else if (position == 5)
            coeff = 204;
        else if (position == 6)
            coeff = 149;
        else if (position == 7)
            coeff = 84;
        else if (position == 8)
            coeff = 9;
        else if (position == 9)
            coeff = -76;
        else
            coeff = -171;
        break;

    case 11:
        if (position == 0)
            coeff = 79;
        else if (position == 1)
            coeff = 78;
        else if (position == 2)
            coeff = 75;
        else if (position == 3)
            coeff = 70;
        else if (position == 4)
            coeff = 63;
        else if (position == 5)
            coeff = 54;
        else if (position == 6)
            coeff = 43;
        else if (position == 7)
            coeff = 30;
        else if (position == 8)
            coeff = 15;
        else if (position == 9)
            coeff = -2;
        else if (position == 10)
            coeff = -21;
        else
            coeff = -42;
        break;

    case 12:
        if (position == 0)
            coeff = 467;
        else if (position == 1)
            coeff = 462;
        else if (position == 2)
            coeff = 447;
        else if (position == 3)
            coeff = 422;
        else if (position == 4)
            coeff = 387;
        else if (position == 5)
            coeff = 343;
        else if (position == 6)
            coeff = 287;
        else if (position == 7)
            coeff = 222;
        else if (position == 8)
            coeff = 147;
        else if (position == 9)
            coeff = 62;
        else if (position == 10)
            coeff = -33;
        else if (position == 11)
            coeff = -138;
        else
            coeff = -253;
        break;
    }
    return coeff;
}

inline double SavitzkyGolayNorm(int points)
{
    double norm;

    if (points == 2)
        norm = 35;
    else if (points == 3)
        norm = 21;
    else if (points == 4)
        norm = 231;
    else if (points == 5)
        norm = 429;
    else if (points == 6)
        norm = 143;
    else if (points == 7)
        norm = 1105;
    else if (points == 8)
        norm = 323;
    else if (points == 9)
        norm = 2261;
    else if (points == 10)
        norm = 3059;
    else if (points == 11)
        norm = 805;
    else if (points == 12)
        norm = 5175;
    return norm;
}
}
