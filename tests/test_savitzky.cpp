/*
 * <Unit tests for Savitzky-Golay filter>
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

#include "libpeakpick/savitzky.h"
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
    std::cout << "Running Savitzky-Golay filter tests..." << std::endl;

    // Test 1: SavitzkyGolayCoefficient with 2 points
    {
        double coeff0 = PeakPick::SavitzkyGolayCoefficient(2, 0);
        double coeff1 = PeakPick::SavitzkyGolayCoefficient(2, 1);

        TEST_ASSERT(coeff0 != 0.0, "Coefficient for 2 points, position 0");
        TEST_ASSERT(coeff1 != 0.0, "Coefficient for 2 points, position 1");
        std::cout << "  ✓ SavitzkyGolayCoefficient with 2 points" << std::endl;
    }

    // Test 2: SavitzkyGolayCoefficient with 3 points
    {
        double coeff0 = PeakPick::SavitzkyGolayCoefficient(3, 0);
        double coeff1 = PeakPick::SavitzkyGolayCoefficient(3, 1);
        double coeff2 = PeakPick::SavitzkyGolayCoefficient(3, 2);

        TEST_ASSERT(coeff0 != 0.0, "Coefficient for 3 points, position 0");
        TEST_ASSERT(coeff1 != 0.0, "Coefficient for 3 points, position 1");
        TEST_ASSERT(coeff2 != 0.0, "Coefficient for 3 points, position 2");
        std::cout << "  ✓ SavitzkyGolayCoefficient with 3 points" << std::endl;
    }

    // Test 3: SavitzkyGolayCoefficient with 4 points
    {
        double coeff0 = PeakPick::SavitzkyGolayCoefficient(4, 0);
        double coeff3 = PeakPick::SavitzkyGolayCoefficient(4, 3);

        TEST_ASSERT(coeff0 != 0.0, "Coefficient for 4 points, position 0");
        TEST_ASSERT(coeff3 != 0.0, "Coefficient for 4 points, position 3");
        std::cout << "  ✓ SavitzkyGolayCoefficient with 4 points" << std::endl;
    }

    // Test 4: SavitzkyGolayCoefficient with 5 points
    {
        double coeff0 = PeakPick::SavitzkyGolayCoefficient(5, 0);
        double coeff4 = PeakPick::SavitzkyGolayCoefficient(5, 4);

        TEST_ASSERT(coeff0 != 0.0, "Coefficient for 5 points, position 0");
        TEST_ASSERT(coeff4 != 0.0, "Coefficient for 5 points, position 4");
        std::cout << "  ✓ SavitzkyGolayCoefficient with 5 points" << std::endl;
    }

    // Test 5: SavitzkyGolayCoefficient out of range
    {
        double coeff = PeakPick::SavitzkyGolayCoefficient(3, 10);
        TEST_ASSERT_NEAR(coeff, 0.0, 0.0001, "Out of range position should return 0");
        std::cout << "  ✓ SavitzkyGolayCoefficient out of range" << std::endl;
    }

    // Test 6: SavitzkyGolayNorm
    {
        double norm2 = PeakPick::SavitzkyGolayNorm(2);
        double norm3 = PeakPick::SavitzkyGolayNorm(3);
        double norm4 = PeakPick::SavitzkyGolayNorm(4);

        TEST_ASSERT(norm2 > 0.0, "Norm for 2 points");
        TEST_ASSERT(norm3 > 0.0, "Norm for 3 points");
        TEST_ASSERT(norm4 > 0.0, "Norm for 4 points");
        std::cout << "  ✓ SavitzkyGolayNorm" << std::endl;
    }

    // Test 7: SavitzkyGolayNorm returns 0 for out of range
    {
        double norm = PeakPick::SavitzkyGolayNorm(100);
        TEST_ASSERT_NEAR(norm, 0.0, 0.0001, "Norm for unsupported points");
        std::cout << "  ✓ SavitzkyGolayNorm for unsupported points" << std::endl;
    }

    std::cout << "\nAll Savitzky-Golay filter tests passed!" << std::endl;
    return 0;
}
