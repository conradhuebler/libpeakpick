/*
 * <Utility functions for libpeakpick>
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

#pragma once

#include <Eigen/Dense>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "spectrum.h"

typedef Eigen::VectorXd Vector;

namespace PeakPick {

/*!
 * \brief Save spectrum to file
 * \param spec Spectrum to save
 * \param filename Output filename
 * \param save_x If true, save both X and Y columns; if false, only Y
 * \return true if successful, false otherwise
 */
inline bool saveToFile(const spectrum& spec, const std::string& filename, bool save_x = true)
{
    std::ofstream outfile(filename);
    if (!outfile.is_open()) {
        std::cerr << "Unable to open file for writing: " << filename << std::endl;
        return false;
    }

    outfile << "# Spectrum data exported from libpeakpick" << std::endl;
    outfile << "# Size: " << spec.size() << std::endl;
    outfile << "# Mean: " << spec.Mean() << std::endl;
    outfile << "# Max: " << spec.Max() << std::endl;
    outfile << "# Min: " << spec.Min() << std::endl;

    if (save_x) {
        outfile << "# Format: X Y" << std::endl;
        for (unsigned int i = 0; i < spec.size(); ++i) {
            outfile << spec.X(i) << " " << spec.Y(i) << std::endl;
        }
    } else {
        outfile << "# Format: Y only" << std::endl;
        outfile << "#start = " << spec.XMin() << std::endl;
        outfile << "#end = " << spec.XMax() << std::endl;
        for (unsigned int i = 0; i < spec.size(); ++i) {
            outfile << spec.Y(i) << std::endl;
        }
    }

    outfile.close();
    return true;
}

/*!
 * \brief Calculate signal-to-noise ratio for a spectrum
 * \param spec Input spectrum
 * \param peak_region_start Start index of peak region
 * \param peak_region_end End index of peak region
 * \param noise_region_start Start index of noise region
 * \param noise_region_end End index of noise region
 * \return Signal-to-noise ratio
 */
inline double calculateSNR(const spectrum& spec,
    unsigned int peak_region_start, unsigned int peak_region_end,
    unsigned int noise_region_start, unsigned int noise_region_end)
{
    if (peak_region_end <= peak_region_start || noise_region_end <= noise_region_start)
        return 0.0;

    // Calculate peak signal (max in peak region)
    double peak_signal = spec.Y(peak_region_start);
    for (unsigned int i = peak_region_start; i < peak_region_end && i < spec.size(); ++i) {
        if (spec.Y(i) > peak_signal)
            peak_signal = spec.Y(i);
    }

    // Calculate noise (stddev in noise region)
    double noise_sum = 0.0;
    double noise_sum_sq = 0.0;
    unsigned int count = 0;

    for (unsigned int i = noise_region_start; i < noise_region_end && i < spec.size(); ++i) {
        double val = spec.Y(i);
        noise_sum += val;
        noise_sum_sq += val * val;
        count++;
    }

    if (count == 0)
        return 0.0;

    double noise_mean = noise_sum / count;
    double noise_variance = (noise_sum_sq / count) - (noise_mean * noise_mean);
    double noise = sqrt(std::max(0.0, noise_variance));

    if (noise < 1e-10)
        return peak_signal > 0 ? 1e10 : 0.0;

    return peak_signal / noise;
}

/*!
 * \brief Find baseline points automatically
 * \param spec Input spectrum
 * \param num_points Number of baseline points to find
 * \param percentile Percentile to use for baseline detection (default 0.1 = 10%)
 * \return Vector of baseline points (indices)
 */
inline std::vector<unsigned int> findBaselinePoints(const spectrum& spec,
    unsigned int num_points,
    double percentile = 0.1)
{
    std::vector<unsigned int> baseline_points;

    if (num_points == 0 || spec.size() == 0)
        return baseline_points;

    // Divide spectrum into segments
    unsigned int segment_size = spec.size() / num_points;
    if (segment_size == 0)
        segment_size = 1;

    for (unsigned int seg = 0; seg < num_points; ++seg) {
        unsigned int start = seg * segment_size;
        unsigned int end = std::min(start + segment_size, spec.size());

        if (start >= spec.size())
            break;

        // Find minimum in this segment
        double min_val = spec.Y(start);
        unsigned int min_idx = start;

        for (unsigned int i = start; i < end; ++i) {
            if (spec.Y(i) < min_val) {
                min_val = spec.Y(i);
                min_idx = i;
            }
        }

        baseline_points.push_back(min_idx);
    }

    return baseline_points;
}

/*!
 * \brief Resample spectrum to new X grid
 * \param spec Input spectrum
 * \param new_x New X grid
 * \return Resampled spectrum
 */
inline spectrum resample(const spectrum& spec, const Vector& new_x)
{
    Vector new_y(new_x.size());

    for (int i = 0; i < new_x.size(); ++i) {
        double x_val = new_x(i);

        // Linear interpolation
        if (x_val <= spec.XMin()) {
            new_y(i) = spec.Y(0);
        } else if (x_val >= spec.XMax()) {
            new_y(i) = spec.Y(spec.size() - 1);
        } else {
            // Find surrounding points
            int idx = spec.XtoIndex(x_val);

            if (idx >= spec.size() - 1) {
                new_y(i) = spec.Y(spec.size() - 1);
            } else {
                double x0 = spec.X(idx);
                double x1 = spec.X(idx + 1);
                double y0 = spec.Y(idx);
                double y1 = spec.Y(idx + 1);

                // Linear interpolation
                double t = (x_val - x0) / (x1 - x0);
                new_y(i) = y0 + t * (y1 - y0);
            }
        }
    }

    return spectrum(new_x, new_y);
}

/*!
 * \brief Calculate full width at half maximum (FWHM) for a peak
 * \param spec Input spectrum
 * \param peak Peak structure with max position
 * \return FWHM value in X units
 */
inline double calculateFWHM(const spectrum& spec, const Peak& peak)
{
    if (peak.max >= spec.size())
        return 0.0;

    double peak_height = spec.Y(peak.max);
    double half_height = peak_height / 2.0;

    // Find left half-maximum point
    unsigned int left_idx = peak.max;
    for (int i = peak.max; i >= 0 && i >= (int)peak.start; --i) {
        if (spec.Y(i) <= half_height) {
            left_idx = i;
            break;
        }
    }

    // Find right half-maximum point
    unsigned int right_idx = peak.max;
    for (unsigned int i = peak.max; i < spec.size() && i < peak.end; ++i) {
        if (spec.Y(i) <= half_height) {
            right_idx = i;
            break;
        }
    }

    double fwhm = spec.X(right_idx) - spec.X(left_idx);
    return fwhm;
}

/*!
 * \brief Subtract one spectrum from another
 * \param spec1 First spectrum
 * \param spec2 Second spectrum (to subtract)
 * \return Difference spectrum
 */
inline spectrum subtract(const spectrum& spec1, const spectrum& spec2)
{
    if (spec1.size() != spec2.size()) {
        std::cerr << "Spectra must have same size for subtraction" << std::endl;
        return spec1;
    }

    Vector new_y(spec1.size());
    for (unsigned int i = 0; i < spec1.size(); ++i) {
        new_y(i) = spec1.Y(i) - spec2.Y(i);
    }

    return spectrum(spec1.x(), new_y);
}

/*!
 * \brief Add two spectra
 * \param spec1 First spectrum
 * \param spec2 Second spectrum
 * \return Sum spectrum
 */
inline spectrum add(const spectrum& spec1, const spectrum& spec2)
{
    if (spec1.size() != spec2.size()) {
        std::cerr << "Spectra must have same size for addition" << std::endl;
        return spec1;
    }

    Vector new_y(spec1.size());
    for (unsigned int i = 0; i < spec1.size(); ++i) {
        new_y(i) = spec1.Y(i) + spec2.Y(i);
    }

    return spectrum(spec1.x(), new_y);
}

/*!
 * \brief Multiply spectrum by a scalar
 * \param spec Input spectrum
 * \param factor Multiplication factor
 * \return Scaled spectrum
 */
inline spectrum scale(const spectrum& spec, double factor)
{
    Vector new_y(spec.size());
    for (unsigned int i = 0; i < spec.size(); ++i) {
        new_y(i) = spec.Y(i) * factor;
    }

    return spectrum(spec.x(), new_y);
}

} // namespace PeakPick
