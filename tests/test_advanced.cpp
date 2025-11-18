/*
 * <Unit tests for advanced algorithms>
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

#include "libpeakpick/advanced.h"
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
    std::cout << "Running advanced algorithm tests..." << std::endl;

    // Test 1: Peak prominence calculation
    {
        std::cout << "  Testing peak prominence..." << std::endl;
        Vector test_x(20), test_y(20);
        for (int i = 0; i < 20; ++i) {
            test_x(i) = i;
            test_y(i) = 1.0; // Baseline

            // Add two peaks with different prominence
            if (i == 5)
                test_y(i) = 5.0; // High peak
            if (i == 15)
                test_y(i) = 2.5; // Lower peak
        }

        PeakPick::spectrum spec(test_x, test_y);

        double prom1 = PeakPick::CalculatePeakProminence(&spec, 5);
        double prom2 = PeakPick::CalculatePeakProminence(&spec, 15);

        TEST_ASSERT(prom1 > prom2, "Higher peak should have higher prominence");
        TEST_ASSERT(prom1 > 3.0, "Prominence of high peak");
        std::cout << "    ✓ Peak prominence calculation (prom1=" << prom1 << ", prom2=" << prom2 << ")" << std::endl;
    }

    // Test 2: Advanced peak picking with prominence filter
    {
        std::cout << "  Testing advanced peak picking..." << std::endl;
        Vector test_x(100), test_y(100);
        for (int i = 0; i < 100; ++i) {
            test_x(i) = i;
            test_y(i) = 0.5 + 0.1 * (rand() % 100) / 100.0; // Noisy baseline

            // Add prominent peaks
            test_y(i) += 5.0 * exp(-0.5 * pow((i - 25.0) / 3.0, 2));  // Peak 1
            test_y(i) += 3.0 * exp(-0.5 * pow((i - 50.0) / 4.0, 2));  // Peak 2
            test_y(i) += 4.0 * exp(-0.5 * pow((i - 75.0) / 3.5, 2)); // Peak 3
        }

        PeakPick::spectrum spec(test_x, test_y);

        // Pick peaks with prominence threshold
        std::vector<PeakPick::Peak> peaks = PeakPick::PickPeaksAdvanced(&spec, 1.0, 1.5, 10);

        TEST_ASSERT(peaks.size() == 3, "Should find 3 prominent peaks");
        std::cout << "    ✓ Advanced peak picking (found " << peaks.size() << " peaks)" << std::endl;
    }

    // Test 3: Noise estimation
    {
        std::cout << "  Testing noise estimation..." << std::endl;
        Vector test_x(100), test_y(100);

        // Create signal with known noise level
        double true_noise = 0.5;
        for (int i = 0; i < 100; ++i) {
            test_x(i) = i;
            test_y(i) = 10.0 + true_noise * (2.0 * (rand() % 100) / 100.0 - 1.0); // Constant + noise
        }

        PeakPick::spectrum spec(test_x, test_y);
        double estimated_noise = PeakPick::EstimateNoise(spec);

        TEST_ASSERT(estimated_noise > 0.0, "Noise estimate should be positive");
        std::cout << "    ✓ Noise estimation (estimated=" << estimated_noise << ")" << std::endl;
    }

    // Test 4: AsLS baseline correction
    {
        std::cout << "  Testing AsLS baseline correction..." << std::endl;
        Vector test_x(100), test_y(100);

        // Create spectrum with polynomial baseline + peak
        for (int i = 0; i < 100; ++i) {
            test_x(i) = i;
            // Polynomial baseline
            test_y(i) = 1.0 + 0.05 * i + 0.001 * i * i;
            // Add peak
            test_y(i) += 5.0 * exp(-0.5 * pow((i - 50.0) / 5.0, 2));
        }

        PeakPick::spectrum spec(test_x, test_y);
        PeakPick::spectrum corrected = PeakPick::BaselineAsLS(spec, 1e6, 0.01, 10);

        TEST_ASSERT(corrected.size() == spec.size(), "Corrected spectrum size");

        // Check that baseline is removed (mean should be close to 0 away from peak)
        double mean_baseline = 0.0;
        int count = 0;
        for (unsigned int i = 0; i < 20; ++i) { // First 20 points, away from peak
            mean_baseline += corrected.Y(i);
            count++;
        }
        mean_baseline /= count;

        TEST_ASSERT(std::abs(mean_baseline) < 1.0, "Baseline should be removed");
        std::cout << "    ✓ AsLS baseline correction (mean baseline=" << mean_baseline << ")" << std::endl;
    }

    // Test 5: Cubic spline interpolation
    {
        std::cout << "  Testing cubic spline interpolation..." << std::endl;

        // Create coarse grid
        Vector test_x_coarse(11), test_y_coarse(11);
        for (int i = 0; i < 11; ++i) {
            test_x_coarse(i) = i * 10.0;
            test_y_coarse(i) = sin(i * 10.0 * 0.1);
        }

        PeakPick::spectrum spec_coarse(test_x_coarse, test_y_coarse);

        // Create fine grid
        Vector test_x_fine(101);
        for (int i = 0; i < 101; ++i) {
            test_x_fine(i) = i;
        }

        // Resample using cubic splines
        PeakPick::spectrum spec_fine = PeakPick::resampleCubic(spec_coarse, test_x_fine);

        TEST_ASSERT(spec_fine.size() == 101, "Resampled spectrum size");

        // Check that interpolated values match original at grid points
        for (int i = 0; i < 11; ++i) {
            int fine_idx = i * 10;
            TEST_ASSERT_NEAR(spec_fine.Y(fine_idx), spec_coarse.Y(i), 0.1, "Interpolation at grid point");
        }

        std::cout << "    ✓ Cubic spline interpolation" << std::endl;
    }

    // Test 6: Cubic spline evaluation
    {
        std::cout << "  Testing cubic spline direct evaluation..." << std::endl;

        Vector x(5), y(5);
        x << 0, 1, 2, 3, 4;
        y << 0, 1, 0, 1, 0; // Oscillating function

        PeakPick::CubicSpline spline(x, y);

        // Evaluate at grid points
        for (int i = 0; i < 5; ++i) {
            double val = spline.evaluate(x(i));
            TEST_ASSERT_NEAR(val, y(i), 0.01, "Spline at grid point");
        }

        // Evaluate between grid points
        double mid_val = spline.evaluate(1.5);
        TEST_ASSERT(mid_val >= 0.0 && mid_val <= 1.0, "Spline interpolated value in range");

        std::cout << "    ✓ Cubic spline direct evaluation" << std::endl;
    }

    // Test 7: Edge cases - empty spectrum
    {
        std::cout << "  Testing edge cases..." << std::endl;

        Vector empty_x(0), empty_y(0);
        PeakPick::spectrum empty_spec;

        // Should not crash
        double noise = PeakPick::EstimateNoise(empty_spec);
        TEST_ASSERT(noise == 0.0, "Noise estimation on empty spectrum");

        std::cout << "    ✓ Edge cases handled" << std::endl;
    }

    std::cout << "\nAll advanced algorithm tests passed!" << std::endl;
    return 0;
}
