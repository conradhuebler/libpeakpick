/*
 * <Unit tests for baseline correction>
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

#include "libpeakpick/baseline.h"
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
    std::cout << "Running baseline tests..." << std::endl;

    // Test 1: BaseLineResult structure
    {
        PeakPick::BaseLineResult result;
        Vector baseline(10);
        for (int i = 0; i < 10; ++i) {
            baseline(i) = i * 0.1;
        }

        result.baselines.push_back(baseline);
        TEST_ASSERT(result.baselines.size() == 1, "BaseLineResult should contain one baseline");
        std::cout << "  ✓ BaseLineResult structure" << std::endl;
    }

    // Test 2: Simple baseline test with spectrum
    {
        Vector test_x(20), test_y(20);
        // Create spectrum with linear baseline + peaks
        for (int i = 0; i < 20; ++i) {
            test_x(i) = i;
            test_y(i) = i * 0.1;  // Linear baseline
            // Add a peak in the middle
            if (i >= 8 && i <= 12) {
                test_y(i) += 2.0 * exp(-0.5 * pow(i - 10, 2));
            }
        }
        PeakPick::spectrum spec(test_x, test_y);

        TEST_ASSERT(spec.size() == 20, "Spectrum size for baseline test");
        TEST_ASSERT(spec.Max() > 1.0, "Spectrum should have peak above baseline");
        std::cout << "  ✓ Baseline test spectrum creation" << std::endl;
    }

    // Test 3: BaseLineResult with multiple baselines
    {
        PeakPick::BaseLineResult result;

        Vector baseline1(10);
        Vector baseline2(10);
        for (int i = 0; i < 10; ++i) {
            baseline1(i) = i;
            baseline2(i) = i * 2;
        }

        result.baselines.push_back(baseline1);
        result.baselines.push_back(baseline2);

        TEST_ASSERT(result.baselines.size() == 2, "Should have two baselines");
        std::cout << "  ✓ Multiple baselines in result" << std::endl;
    }

    // Test 4: Grid points storage
    {
        PeakPick::BaseLineResult result;

        Vector x_grid(5);
        Vector y_grid(5);
        for (int i = 0; i < 5; ++i) {
            x_grid(i) = i * 2.0;
            y_grid(i) = i * 0.5;
        }

        result.x_grid_points.push_back(x_grid);
        result.y_grid_points.push_back(y_grid);

        TEST_ASSERT(result.x_grid_points.size() == 1, "Should have x grid points");
        TEST_ASSERT(result.y_grid_points.size() == 1, "Should have y grid points");
        std::cout << "  ✓ Grid points storage" << std::endl;
    }

    std::cout << "\nAll baseline tests passed!" << std::endl;
    return 0;
}
