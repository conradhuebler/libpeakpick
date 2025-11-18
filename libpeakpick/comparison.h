/*
 * <Spectrum comparison and metrics>
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

#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

#include "analyse.h"
#include "mathhelper.h"
#include "spectrum.h"
#include "utilities.h"

typedef Eigen::VectorXd Vector;

namespace PeakPick {

/*!
 * \brief Calculate Root Mean Square Deviation between two spectra
 * \param spec1 First spectrum
 * \param spec2 Second spectrum
 * \return RMSD value
 */
inline double RMSD(const spectrum& spec1, const spectrum& spec2)
{
    if (spec1.size() != spec2.size()) {
        std::cerr << "Warning: RMSD - spectra have different sizes, resampling spec2" << std::endl;
        spectrum spec2_resampled = resample(spec2, spec1.x());
        return RMSD(spec1, spec2_resampled);
    }

    double sum = 0.0;
    for (unsigned int i = 0; i < spec1.size(); ++i) {
        double diff = spec1.Y(i) - spec2.Y(i);
        sum += diff * diff;
    }

    return sqrt(sum / spec1.size());
}

/*!
 * \brief Calculate Pearson correlation coefficient between two spectra
 * \param spec1 First spectrum
 * \param spec2 Second spectrum
 * \return Correlation coefficient [-1, 1]
 */
inline double CorrelationCoefficient(const spectrum& spec1, const spectrum& spec2)
{
    if (spec1.size() != spec2.size()) {
        std::cerr << "Warning: CorrelationCoefficient - spectra have different sizes" << std::endl;
        spectrum spec2_resampled = resample(spec2, spec1.x());
        return CorrelationCoefficient(spec1, spec2_resampled);
    }

    double mean1 = spec1.Mean();
    double mean2 = spec2.Mean();

    double numerator = 0.0;
    double denom1 = 0.0;
    double denom2 = 0.0;

    for (unsigned int i = 0; i < spec1.size(); ++i) {
        double diff1 = spec1.Y(i) - mean1;
        double diff2 = spec2.Y(i) - mean2;

        numerator += diff1 * diff2;
        denom1 += diff1 * diff1;
        denom2 += diff2 * diff2;
    }

    if (denom1 == 0.0 || denom2 == 0.0)
        return 0.0;

    return numerator / sqrt(denom1 * denom2);
}

/*!
 * \brief Calculate cosine similarity between two spectra
 * \param spec1 First spectrum
 * \param spec2 Second spectrum
 * \return Cosine similarity [0, 1]
 */
inline double CosineSimilarity(const spectrum& spec1, const spectrum& spec2)
{
    if (spec1.size() != spec2.size()) {
        spectrum spec2_resampled = resample(spec2, spec1.x());
        return CosineSimilarity(spec1, spec2_resampled);
    }

    double dot_product = 0.0;
    double norm1 = 0.0;
    double norm2 = 0.0;

    for (unsigned int i = 0; i < spec1.size(); ++i) {
        dot_product += spec1.Y(i) * spec2.Y(i);
        norm1 += spec1.Y(i) * spec1.Y(i);
        norm2 += spec2.Y(i) * spec2.Y(i);
    }

    if (norm1 == 0.0 || norm2 == 0.0)
        return 0.0;

    return dot_product / (sqrt(norm1) * sqrt(norm2));
}

/*!
 * \brief Calculate Mean Absolute Error between two spectra
 * \param spec1 First spectrum
 * \param spec2 Second spectrum
 * \return MAE value
 */
inline double MAE(const spectrum& spec1, const spectrum& spec2)
{
    if (spec1.size() != spec2.size()) {
        spectrum spec2_resampled = resample(spec2, spec1.x());
        return MAE(spec1, spec2_resampled);
    }

    double sum = 0.0;
    for (unsigned int i = 0; i < spec1.size(); ++i) {
        sum += std::abs(spec1.Y(i) - spec2.Y(i));
    }

    return sum / spec1.size();
}

/*!
 * \brief Align two spectra by finding optimal X-shift
 * \param spec1 Reference spectrum (not modified)
 * \param spec2 Spectrum to align
 * \param search_range Maximum shift to search (in X units)
 * \param step Step size for search
 * \return Aligned spectrum (spec2 shifted)
 */
inline spectrum AlignSpectra(const spectrum& spec1, const spectrum& spec2,
    double search_range = 10.0,
    double step = 0.1)
{
    double best_shift = 0.0;
    double best_correlation = -1.0;

    std::cout << "AlignSpectra: Searching for optimal shift in range [" << -search_range << ", " << search_range << "]" << std::endl;

    for (double shift = -search_range; shift <= search_range; shift += step) {
        // Create shifted X-axis
        Vector shifted_x(spec2.size());
        for (unsigned int i = 0; i < spec2.size(); ++i) {
            shifted_x(i) = spec2.X(i) + shift;
        }

        // Resample to spec1's grid
        spectrum shifted_spec(shifted_x, spec2.y());
        spectrum resampled = resample(shifted_spec, spec1.x());

        // Calculate correlation
        double corr = CorrelationCoefficient(spec1, resampled);

        if (corr > best_correlation) {
            best_correlation = corr;
            best_shift = shift;
        }
    }

    std::cout << "AlignSpectra: Best shift = " << best_shift << " (correlation = " << best_correlation << ")" << std::endl;

    // Apply best shift
    Vector aligned_x(spec2.size());
    for (unsigned int i = 0; i < spec2.size(); ++i) {
        aligned_x(i) = spec2.X(i) + best_shift;
    }

    return spectrum(aligned_x, spec2.y());
}

/*!
 * \brief Calculate spectral angle between two spectra
 * Useful for spectral library matching
 * \param spec1 First spectrum
 * \param spec2 Second spectrum
 * \return Spectral angle in degrees [0, 90]
 */
inline double SpectralAngle(const spectrum& spec1, const spectrum& spec2)
{
    double cos_sim = CosineSimilarity(spec1, spec2);
    cos_sim = std::max(-1.0, std::min(1.0, cos_sim)); // Clamp to [-1, 1]
    double angle_rad = acos(cos_sim);
    return angle_rad * 180.0 / 3.14159265; // Convert to degrees
}

/*!
 * \brief Compare two peak lists
 * \param peaks1 First peak list
 * \param peaks2 Second peak list
 * \param tolerance Position tolerance for matching
 * \return Number of matched peaks
 */
inline int ComparePeaks(const std::vector<Peak>& peaks1,
    const std::vector<Peak>& peaks2,
    double tolerance = 1.0)
{
    int matched = 0;

    for (const auto& p1 : peaks1) {
        for (const auto& p2 : peaks2) {
            if (std::abs(p1.deconv_x - p2.deconv_x) < tolerance) {
                matched++;
                break;
            }
        }
    }

    std::cout << "ComparePeaks: " << matched << " matched out of " << peaks1.size() << " peaks (tolerance=" << tolerance << ")" << std::endl;

    return matched;
}

} // namespace PeakPick
