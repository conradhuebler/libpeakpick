/*
 * <Unit tests for spectrum class>
 * Copyright (C) 2024  Conrad Hübler <Conrad.Huebler@gmx.net>
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

#include "libpeakpick/spectrum.h"
#include <iostream>
#include <cmath>

#define TEST_ASSERT(condition, message) \
    if (!(condition)) { \
        std::cerr << "FAILED: " << message << std::endl; \
        return 1; \
    }

#define TEST_ASSERT_NEAR(val1, val2, epsilon, message) \
    if (std::abs((val1) - (val2)) > (epsilon)) { \
        std::cerr << "FAILED: " << message << " (" << val1 << " != " << val2 << ")" << std::endl; \
        return 1; \
    }

int main() {
    std::cout << "Running spectrum tests..." << std::endl;

    // Test data
    Vector test_y(10);
    test_y << 1, 1, 2, 3, 4, 5, 4, 3, 2, 1;

    Vector test_x(10);
    test_x << 0, 1, 2, 3, 4, 5, 6, 7, 8, 9;

    // Test 1: Constructor with X and Y
    {
        PeakPick::spectrum spec(test_x, test_y);
        TEST_ASSERT(spec.size() == 10, "Constructor with X,Y: size");
        TEST_ASSERT_NEAR(spec.X(0), 0.0, 0.0001, "Constructor with X,Y: X(0)");
        TEST_ASSERT_NEAR(spec.Y(0), 1.0, 0.0001, "Constructor with X,Y: Y(0)");
        std::cout << "  ✓ Constructor with X and Y" << std::endl;
    }

    // Test 2: Constructor with Y and range
    {
        PeakPick::spectrum spec(test_y, 0, 10);
        TEST_ASSERT(spec.size() == 10, "Constructor with Y and range: size");
        TEST_ASSERT_NEAR(spec.XMin(), 0.0, 0.0001, "Constructor with Y and range: XMin");
        TEST_ASSERT(spec.XMax() > 9.0, "Constructor with Y and range: XMax");
        std::cout << "  ✓ Constructor with Y and range" << std::endl;
    }

    // Test 3: Constructor throws on size mismatch
    {
        Vector wrong_x(5);
        wrong_x << 1, 2, 3, 4, 5;
        bool threw_exception = false;
        try {
            PeakPick::spectrum spec(wrong_x, test_y);
        } catch (...) {
            threw_exception = true;
        }
        TEST_ASSERT(threw_exception, "Constructor should throw on size mismatch");
        std::cout << "  ✓ Constructor throws on size mismatch" << std::endl;
    }

    // Test 4: Copy constructor
    {
        PeakPick::spectrum spec1(test_x, test_y);
        PeakPick::spectrum spec2(spec1);
        TEST_ASSERT(spec2.size() == spec1.size(), "Copy constructor: size");
        TEST_ASSERT_NEAR(spec2.Mean(), spec1.Mean(), 0.0001, "Copy constructor: Mean");
        std::cout << "  ✓ Copy constructor" << std::endl;
    }

    // Test 5: Mean calculation
    {
        PeakPick::spectrum spec(test_x, test_y);
        // Mean of [1,1,2,3,4,5,4,3,2,1] = 26/10 = 2.6
        TEST_ASSERT_NEAR(spec.Mean(), 2.6, 0.0001, "Mean calculation");
        std::cout << "  ✓ Mean calculation" << std::endl;
    }

    // Test 6: Max and Min
    {
        PeakPick::spectrum spec(test_x, test_y);
        TEST_ASSERT_NEAR(spec.Max(), 5.0, 0.0001, "Max value");
        TEST_ASSERT_NEAR(spec.Min(), 1.0, 0.0001, "Min value");
        TEST_ASSERT(spec.IndexMax() == 5, "Max index");
        TEST_ASSERT(spec.IndexMin() == 0, "Min index");
        std::cout << "  ✓ Max and Min calculations" << std::endl;
    }

    // Test 7: Standard deviation
    {
        PeakPick::spectrum spec(test_x, test_y);
        TEST_ASSERT(spec.StdDev() > 0.0, "StdDev should be > 0");
        std::cout << "  ✓ Standard deviation calculation" << std::endl;
    }

    // Test 8: X and Y accessors
    {
        PeakPick::spectrum spec(test_x, test_y);
        for (unsigned int i = 0; i < 10; ++i) {
            TEST_ASSERT_NEAR(spec.X(i), test_x(i), 0.0001, "X accessor");
            TEST_ASSERT_NEAR(spec.Y(i), test_y(i), 0.0001, "Y accessor");
        }
        std::cout << "  ✓ X and Y accessors" << std::endl;
    }

    // Test 9: Out of bounds access
    {
        PeakPick::spectrum spec(test_x, test_y);
        TEST_ASSERT_NEAR(spec.Y(100), 0.0, 0.0001, "Out of bounds Y");
        TEST_ASSERT_NEAR(spec.X(100), 0.0, 0.0001, "Out of bounds X");
        std::cout << "  ✓ Out of bounds access" << std::endl;
    }

    // Test 10: setY
    {
        PeakPick::spectrum spec(test_x, test_y);
        spec.setY(5, 10.0);
        TEST_ASSERT_NEAR(spec.Y(5), 10.0, 0.0001, "setY");
        std::cout << "  ✓ setY method" << std::endl;
    }

    // Test 11: Step size
    {
        PeakPick::spectrum spec(test_x, test_y);
        TEST_ASSERT_NEAR(spec.Step(), 1.0, 0.0001, "Step size");
        std::cout << "  ✓ Step size calculation" << std::endl;
    }

    // Test 12: XtoIndex
    {
        PeakPick::spectrum spec(test_x, test_y);
        TEST_ASSERT(spec.XtoIndex(0.0) == 0, "XtoIndex(0)");
        TEST_ASSERT(spec.XtoIndex(5.0) == 5, "XtoIndex(5)");
        std::cout << "  ✓ XtoIndex method" << std::endl;
    }

    // Test 13: Center
    {
        PeakPick::spectrum spec(test_x, test_y);
        spec.center();
        TEST_ASSERT_NEAR(spec.Mean(), 0.0, 0.0001, "Center method");
        std::cout << "  ✓ Center method" << std::endl;
    }

    // Test 14: InvertSgn
    {
        PeakPick::spectrum spec(test_x, test_y);
        double original_max = spec.Max();
        spec.InvertSgn();
        TEST_ASSERT_NEAR(spec.Min(), -original_max, 0.0001, "InvertSgn");
        std::cout << "  ✓ InvertSgn method" << std::endl;
    }

    // Test 15: LastY
    {
        PeakPick::spectrum spec(test_x, test_y);
        TEST_ASSERT_NEAR(spec.LastY(), 1.0, 0.0001, "LastY");
        std::cout << "  ✓ LastY method" << std::endl;
    }

    std::cout << "\nAll spectrum tests passed!" << std::endl;
    return 0;
}
