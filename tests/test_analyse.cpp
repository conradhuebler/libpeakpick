/*
 * <Unit tests for analyse functions>
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

#include "libpeakpick/analyse.h"
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
    std::cout << "Running analyse tests..." << std::endl;

    // Test 1: Normalise spectrum
    {
        Vector test_x(10), test_y(10);
        for (int i = 0; i < 10; ++i) {
            test_x(i) = i;
            test_y(i) = i + 1;
        }
        PeakPick::spectrum spec(test_x, test_y);
        PeakPick::Normalise(&spec);

        TEST_ASSERT_NEAR(spec.Max(), 1.0, 0.0001, "Normalised max should be 1.0");
        std::cout << "  ✓ Normalise spectrum" << std::endl;
    }

    // Test 2: Find maximum in peak
    {
        Vector test_x(10), test_y(10);
        test_y << 1, 2, 5, 8, 10, 7, 4, 2, 1, 0;
        for (int i = 0; i < 10; ++i) {
            test_x(i) = i;
        }
        PeakPick::spectrum spec(test_x, test_y);

        PeakPick::Peak peak;
        peak.start = 0;
        peak.end = 10;

        int max_pos = PeakPick::FindMaximum(&spec, peak);
        TEST_ASSERT(max_pos == 4, "Maximum should be at index 4");
        std::cout << "  ✓ Find maximum in peak" << std::endl;
    }

    // Test 3: Find minimum in peak
    {
        Vector test_x(10), test_y(10);
        test_y << 10, 8, 5, 3, 1, 2, 4, 6, 8, 10;
        for (int i = 0; i < 10; ++i) {
            test_x(i) = i;
        }
        PeakPick::spectrum spec(test_x, test_y);

        PeakPick::Peak peak;
        peak.start = 0;
        peak.end = 10;

        int min_pos = PeakPick::FindMinimum(&spec, peak);
        TEST_ASSERT(min_pos == 4, "Minimum should be at index 4");
        std::cout << "  ✓ Find minimum in peak" << std::endl;
    }

    // Test 4: Peak picking
    {
        Vector test_x(20), test_y(20);
        // Create a spectrum with two peaks
        for (int i = 0; i < 20; ++i) {
            test_x(i) = i;
            if (i < 5)
                test_y(i) = 0;
            else if (i >= 5 && i < 8)
                test_y(i) = i - 4;  // First peak
            else if (i >= 8 && i < 12)
                test_y(i) = 12 - i;
            else if (i >= 12 && i < 15)
                test_y(i) = i - 11;  // Second peak
            else
                test_y(i) = 0;
        }
        PeakPick::spectrum spec(test_x, test_y);

        std::vector<PeakPick::Peak> peaks = PeakPick::PickPeaks(&spec, 0.5);
        TEST_ASSERT(peaks.size() > 0, "Should find at least one peak");
        std::cout << "  ✓ Peak picking" << std::endl;
    }

    // Test 5: Numerical integration
    {
        Vector test_x(11), test_y(11);
        // Simple integration: integrate y=2 from 0 to 10, should be ~20
        for (int i = 0; i < 11; ++i) {
            test_x(i) = i;
            test_y(i) = 2.0;
        }
        PeakPick::spectrum spec(test_x, test_y);

        double integ = PeakPick::IntegrateNumerical(&spec, 0, 10);
        TEST_ASSERT_NEAR(integ, 20.0, 0.5, "Integration of constant function");
        std::cout << "  ✓ Numerical integration" << std::endl;
    }

    // Test 6: Integration with offset
    {
        Vector test_x(11), test_y(11);
        for (int i = 0; i < 11; ++i) {
            test_x(i) = i;
            test_y(i) = 5.0;
        }
        PeakPick::spectrum spec(test_x, test_y);

        // Integrate y=5 with offset=5, should be ~0
        double integ = PeakPick::IntegrateNumerical(&spec, 0, 10, 5.0);
        TEST_ASSERT_NEAR(integ, 0.0, 0.5, "Integration with offset");
        std::cout << "  ✓ Integration with offset" << std::endl;
    }

    // Test 7: Integration using Peak structure
    {
        Vector test_x(11), test_y(11);
        for (int i = 0; i < 11; ++i) {
            test_x(i) = i;
            test_y(i) = 3.0;
        }
        PeakPick::spectrum spec(test_x, test_y);

        PeakPick::Peak peak;
        peak.int_start = 0;
        peak.int_end = 10;

        double integ = PeakPick::IntegrateNumerical(&spec, peak);
        TEST_ASSERT_NEAR(integ, 30.0, 0.5, "Integration using Peak");
        TEST_ASSERT_NEAR(peak.integ_num, 30.0, 0.5, "Peak integ_num should be set");
        std::cout << "  ✓ Integration using Peak structure" << std::endl;
    }

    // Test 8: Integration with std::vector
    {
        std::vector<double> x = {0, 1, 2, 3, 4, 5};
        std::vector<double> y = {1, 1, 1, 1, 1, 1};

        double integ = PeakPick::IntegrateNumerical(x, y);
        TEST_ASSERT_NEAR(integ, 5.0, 0.5, "Integration with std::vector");
        std::cout << "  ✓ Integration with std::vector" << std::endl;
    }

    // Test 9: Peak structure setters
    {
        PeakPick::Peak peak;
        peak.setPeakStart(10);
        peak.setPeakEnd(20);

        TEST_ASSERT(peak.start == 10, "Peak start");
        TEST_ASSERT(peak.end == 20, "Peak end");
        TEST_ASSERT(peak.int_start == 10, "Peak int_start");
        TEST_ASSERT(peak.int_end == 20, "Peak int_end");
        std::cout << "  ✓ Peak structure setters" << std::endl;
    }

    std::cout << "\nAll analyse tests passed!" << std::endl;
    return 0;
}
