/*
 * <Unit tests for mathhelper functions>
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

#include "libpeakpick/mathhelper.h"
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
    std::cout << "Running mathhelper tests..." << std::endl;

    // Test 1: Mean calculation
    {
        Vector v(5);
        v << 1, 2, 3, 4, 5;
        double m = PeakPick::mean(v);
        TEST_ASSERT_NEAR(m, 3.0, 0.0001, "Mean of 1,2,3,4,5");
        std::cout << "  ✓ Mean calculation" << std::endl;
    }

    // Test 2: Mean with min/max tracking
    {
        Vector v(5);
        v << 1, 5, 3, 2, 4;
        unsigned int pos_min = 0, pos_max = 0;
        double m = PeakPick::mean(v, &pos_min, &pos_max);
        TEST_ASSERT_NEAR(m, 3.0, 0.0001, "Mean with min/max");
        TEST_ASSERT(pos_min == 0, "Min position");
        TEST_ASSERT(pos_max == 1, "Max position");
        std::cout << "  ✓ Mean with min/max tracking" << std::endl;
    }

    // Test 3: Standard deviation
    {
        Vector v(5);
        v << 2, 4, 4, 4, 5;
        double m = PeakPick::mean(v);
        double sd = PeakPick::stddev(v, m);
        TEST_ASSERT(sd > 0.0, "StdDev should be > 0");
        std::cout << "  ✓ Standard deviation" << std::endl;
    }

    // Test 4: Gaussian function
    {
        double val = PeakPick::Gaussian(0.0, 1.0, 0.0, 1.0);
        TEST_ASSERT_NEAR(val, 1.0, 0.0001, "Gaussian at peak");
        std::cout << "  ✓ Gaussian function" << std::endl;
    }

    // Test 5: Lorentzian function
    {
        double val = PeakPick::Lorentzian(0.0, 0.0, 1.0);
        TEST_ASSERT(val > 0.0, "Lorentzian at peak");
        std::cout << "  ✓ Lorentzian function" << std::endl;
    }

    // Test 6: Polynomial evaluation
    {
        Vector coeff(3);
        coeff << 1, 2, 3;  // 1 + 2x + 3x^2
        double val = PeakPick::Polynomial(2.0, coeff);
        // 1 + 2*2 + 3*4 = 1 + 4 + 12 = 17
        TEST_ASSERT_NEAR(val, 17.0, 0.0001, "Polynomial evaluation");
        std::cout << "  ✓ Polynomial evaluation" << std::endl;
    }

    // Test 7: Linear regression
    {
        Vector x(5), y(5);
        x << 1, 2, 3, 4, 5;
        y << 2, 4, 6, 8, 10;  // Perfect line: y = 2x

        PeakPick::LinearRegression reg = PeakPick::LeastSquares(x, y);
        TEST_ASSERT_NEAR(reg.m, 2.0, 0.01, "Linear regression slope");
        TEST_ASSERT_NEAR(reg.n, 0.0, 0.01, "Linear regression intercept");
        TEST_ASSERT(reg.R > 0.99, "R² should be close to 1.0");
        std::cout << "  ✓ Linear regression" << std::endl;
    }

    // Test 8: Linear regression with vector input
    {
        std::vector<double> x_vec = {1, 2, 3, 4, 5};
        std::vector<double> y_vec = {1, 3, 5, 7, 9};  // y = 2x - 1

        PeakPick::LinearRegression reg = PeakPick::LeastSquares(x_vec, y_vec);
        TEST_ASSERT_NEAR(reg.m, 2.0, 0.01, "Linear regression slope (std::vector)");
        std::cout << "  ✓ Linear regression with std::vector" << std::endl;
    }

    // Test 9: Mean of empty vector
    {
        Vector v(0);
        double m = PeakPick::mean(v);
        TEST_ASSERT_NEAR(m, 0.0, 0.0001, "Mean of empty vector");
        std::cout << "  ✓ Mean of empty vector" << std::endl;
    }

    // Test 10: StdDev of empty vector
    {
        Vector v(0);
        double sd = PeakPick::stddev(v, 0.0);
        TEST_ASSERT_NEAR(sd, 0.0, 0.0001, "StdDev of empty vector");
        std::cout << "  ✓ StdDev of empty vector" << std::endl;
    }

    std::cout << "\nAll mathhelper tests passed!" << std::endl;
    return 0;
}
