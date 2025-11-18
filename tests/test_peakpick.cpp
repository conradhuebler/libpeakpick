/*
 * <Unit tests for peakpick file loading>
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
#include <fstream>
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
    std::cout << "Running peakpick file loading tests..." << std::endl;

    // Test 1: Create a test file and load it
    {
        std::ofstream testfile("test_spectrum.dat");
        TEST_ASSERT(testfile.is_open(), "Test file should be created");

        testfile << "# Test spectrum file\n";
        testfile << "1.0\n";
        testfile << "2.0\n";
        testfile << "3.0\n";
        testfile << "4.0\n";
        testfile << "5.0\n";
        testfile.close();

        PeakPick::spectrum spec = loadFromFile("test_spectrum.dat");
        TEST_ASSERT(spec.size() == 5, "Loaded spectrum should have 5 points");
        TEST_ASSERT_NEAR(spec.Y(0), 1.0, 0.0001, "First Y value");
        TEST_ASSERT_NEAR(spec.Y(4), 5.0, 0.0001, "Last Y value");

        std::remove("test_spectrum.dat");
        std::cout << "  ✓ Load spectrum from file" << std::endl;
    }

    // Test 2: Load file with comments
    {
        std::ofstream testfile("test_spectrum_comments.dat");
        testfile << "# This is a comment\n";
        testfile << "# Another comment\n";
        testfile << "10.0\n";
        testfile << "20.0\n";
        testfile << "# Mid-file comment\n";
        testfile << "30.0\n";
        testfile.close();

        PeakPick::spectrum spec = loadFromFile("test_spectrum_comments.dat");
        TEST_ASSERT(spec.size() == 3, "Should have 3 data points (comments ignored)");
        TEST_ASSERT_NEAR(spec.Y(0), 10.0, 0.0001, "First data point");
        TEST_ASSERT_NEAR(spec.Y(2), 30.0, 0.0001, "Third data point");

        std::remove("test_spectrum_comments.dat");
        std::cout << "  ✓ Load file with comments" << std::endl;
    }

    // Test 3: Load non-existent file
    {
        PeakPick::spectrum spec = loadFromFile("nonexistent_file.dat");
        TEST_ASSERT(spec.size() == 0, "Non-existent file should result in empty spectrum");
        std::cout << "  ✓ Handle non-existent file" << std::endl;
    }

    // Test 4: Load file with min/max range
    {
        std::ofstream testfile("test_spectrum_range.dat");
        for (int i = 0; i < 10; ++i) {
            testfile << (double)i << "\n";
        }
        testfile.close();

        PeakPick::spectrum spec = loadFromFile("test_spectrum_range.dat", 100.0, 200.0);
        TEST_ASSERT(spec.size() == 10, "Should load all 10 points");
        TEST_ASSERT_NEAR(spec.XMin(), 100.0, 0.1, "X should start at min");

        std::remove("test_spectrum_range.dat");
        std::cout << "  ✓ Load file with custom X range" << std::endl;
    }

    // Test 5: Load larger file
    {
        std::ofstream testfile("test_spectrum_large.dat");
        for (int i = 0; i < 100; ++i) {
            testfile << sin(i * 0.1) << "\n";
        }
        testfile.close();

        PeakPick::spectrum spec = loadFromFile("test_spectrum_large.dat");
        TEST_ASSERT(spec.size() == 100, "Should load 100 points");
        TEST_ASSERT(spec.Mean() != 0.0, "Mean should be calculated");

        std::remove("test_spectrum_large.dat");
        std::cout << "  ✓ Load larger file" << std::endl;
    }

    std::cout << "\nAll peakpick file loading tests passed!" << std::endl;
    return 0;
}
