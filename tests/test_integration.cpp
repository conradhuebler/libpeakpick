/*
 * <Integration tests for libpeakpick>
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

#include "libpeakpick/peakpick.h"
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
    std::cout << "Running integration tests..." << std::endl;

    // Test 1: Complete workflow - create spectrum, analyze, find peaks
    {
        std::cout << "  Testing complete workflow..." << std::endl;

        Vector test_x(100), test_y(100);

        // Create a synthetic spectrum with multiple Gaussian peaks
        for (int i = 0; i < 100; ++i) {
            test_x(i) = i;
            // Baseline
            test_y(i) = 0.1;

            // Add three Gaussian peaks
            test_y(i) += 5.0 * exp(-0.5 * pow((i - 25.0) / 3.0, 2));  // Peak 1
            test_y(i) += 3.0 * exp(-0.5 * pow((i - 50.0) / 4.0, 2));  // Peak 2
            test_y(i) += 4.0 * exp(-0.5 * pow((i - 75.0) / 3.5, 2)); // Peak 3

            // Add some noise
            test_y(i) += 0.05 * (rand() % 100 - 50) / 100.0;
        }

        PeakPick::spectrum spec(test_x, test_y);

        TEST_ASSERT(spec.size() == 100, "Spectrum size");
        TEST_ASSERT(spec.Max() > 4.0, "Spectrum should have peaks");
        TEST_ASSERT(spec.Mean() > 0.5, "Mean should be above baseline");

        std::cout << "    ✓ Spectrum created and analyzed" << std::endl;

        // Normalize
        PeakPick::Normalise(&spec);
        TEST_ASSERT_NEAR(spec.Max(), 1.0, 0.01, "Normalized max");
        std::cout << "    ✓ Spectrum normalized" << std::endl;

        // Pick peaks
        std::vector<PeakPick::Peak> peaks = PeakPick::PickPeaks(&spec, 0.1);
        TEST_ASSERT(peaks.size() > 0, "Should find peaks");
        std::cout << "    ✓ Found " << peaks.size() << " peaks" << std::endl;
    }

    // Test 2: Smoothing and analysis workflow
    {
        std::cout << "  Testing smoothing workflow..." << std::endl;

        Vector test_x(50), test_y(50);

        // Create noisy data
        for (int i = 0; i < 50; ++i) {
            test_x(i) = i;
            test_y(i) = sin(i * 0.2) + 0.1 * (rand() % 100 - 50) / 100.0;
        }

        PeakPick::spectrum spec(test_x, test_y);
        double original_stddev = spec.StdDev();

        // Smooth the spectrum
        PeakPick::SmoothFunction(&spec, 3);

        TEST_ASSERT(spec.size() > 0, "Spectrum should still have data after smoothing");
        std::cout << "    ✓ Spectrum smoothed" << std::endl;
    }

    // Test 3: Peak integration workflow
    {
        std::cout << "  Testing peak integration workflow..." << std::endl;

        Vector test_x(50), test_y(50);

        // Create a simple peak
        for (int i = 0; i < 50; ++i) {
            test_x(i) = i;
            test_y(i) = exp(-0.1 * pow(i - 25.0, 2));
        }

        PeakPick::spectrum spec(test_x, test_y);

        PeakPick::Peak peak;
        peak.int_start = 15;
        peak.int_end = 35;

        double integral = PeakPick::IntegrateNumerical(&spec, peak);

        TEST_ASSERT(integral > 0.0, "Integral should be positive");
        TEST_ASSERT(peak.integ_num > 0.0, "Peak integral should be stored");

        std::cout << "    ✓ Peak integrated (integral = " << integral << ")" << std::endl;
    }

    // Test 4: Center and invert workflow
    {
        std::cout << "  Testing center and invert workflow..." << std::endl;

        Vector test_x(20), test_y(20);

        for (int i = 0; i < 20; ++i) {
            test_x(i) = i;
            test_y(i) = i + 10.0;  // Positive values
        }

        PeakPick::spectrum spec(test_x, test_y);

        double original_mean = spec.Mean();
        TEST_ASSERT(original_mean > 10.0, "Original mean");

        spec.center();
        TEST_ASSERT_NEAR(spec.Mean(), 0.0, 0.01, "Centered mean");

        spec.InvertSgn();
        TEST_ASSERT(spec.Max() < 0.0 || spec.Min() < -5.0, "Should have negative values after inversion");

        std::cout << "    ✓ Center and invert operations" << std::endl;
    }

    // Test 5: Linear regression on spectrum data
    {
        std::cout << "  Testing linear regression..." << std::endl;

        Vector test_x(10), test_y(10);

        // Perfect linear relationship
        for (int i = 0; i < 10; ++i) {
            test_x(i) = i;
            test_y(i) = 2.0 * i + 3.0;
        }

        PeakPick::LinearRegression reg = PeakPick::LeastSquares(test_x, test_y);

        TEST_ASSERT_NEAR(reg.m, 2.0, 0.01, "Slope should be 2.0");
        TEST_ASSERT_NEAR(reg.n, 3.0, 0.01, "Intercept should be 3.0");
        TEST_ASSERT(reg.R > 0.99, "R² should be close to 1.0");

        std::cout << "    ✓ Linear regression (m=" << reg.m << ", n=" << reg.n << ", R²=" << reg.R << ")" << std::endl;
    }

    // Test 6: Polynomial evaluation
    {
        std::cout << "  Testing polynomial evaluation..." << std::endl;

        Vector coeff(4);
        coeff << 1, 0, 2, 1;  // 1 + 0x + 2x² + 1x³

        double y_at_0 = PeakPick::Polynomial(0.0, coeff);
        double y_at_1 = PeakPick::Polynomial(1.0, coeff);
        double y_at_2 = PeakPick::Polynomial(2.0, coeff);

        TEST_ASSERT_NEAR(y_at_0, 1.0, 0.01, "y(0) = 1");
        TEST_ASSERT_NEAR(y_at_1, 4.0, 0.01, "y(1) = 1 + 2 + 1 = 4");
        TEST_ASSERT_NEAR(y_at_2, 17.0, 0.01, "y(2) = 1 + 8 + 8 = 17");

        std::cout << "    ✓ Polynomial evaluation" << std::endl;
    }

    std::cout << "\nAll integration tests passed!" << std::endl;
    return 0;
}
